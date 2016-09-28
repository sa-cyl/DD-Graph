/*
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * The basic GraphChi engine.
 */


#define maxi(a,b) ( ((a)>(b)) ? (a):(b) )
#define mini(a,b) ( ((a)>(b)) ? (b):(a) )

#ifndef DEF_GRAPHCHI_GRAPHCHI_ENGINE
#define DEF_GRAPHCHI_GRAPHCHI_ENGINE

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <sys/time.h>
#include <errno.h>

#include "api/chifilenames.hpp"
#include "api/graph_objects.hpp"
#include "api/ddgraph_context.hpp"
#include "api/ddgraph_program.hpp"
#include "engine/auxdata/degree_data.hpp"
#include "engine/auxdata/vertex_data.hpp"
#include "engine/bitset_scheduler.hpp"
#include "comm/xgworker.hpp"
#include "io/stripedio.hpp"
#include "logger/logger.hpp"
#include "metrics/metrics.hpp"
#include "shards/memoryshard.hpp"
#include "shards/slidingshard.hpp"
#include "util/pthread_tools.hpp"
#include "output/output.hpp"
#include "schedule/task.hpp"
#include "util/dense_bitset.hpp"
#include "util/toplist.hpp"
#include "shards/dynamicdata/dynamicblock.hpp"

/* Modified By YouLi Cheng */
static task *T;
static intervallock *iL;

namespace ddgraph {
    
    template <typename VertexDataType, typename EdgeDataType,   
    typename svertex_t = ddgraph_vertex<VertexDataType, EdgeDataType> >
    
    class ddgraph_engine {
    public:     
        typedef sliding_shard<VertexDataType, EdgeDataType, svertex_t> slidingshard_t;
        typedef memory_shard<VertexDataType, EdgeDataType, svertex_t> memshard_t;
        
    protected:
        std::string base_filename;
        int nshards;
        
        /* IO manager */
        stripedio * iomgr;
        
        /* Shards */
        std::vector<slidingshard_t *> sliding_shards;
        memshard_t * memoryshard;
        std::vector<std::pair<vid_t, vid_t> > intervals;
        
        /* Auxilliary data handlers */
        degree_data * degree_handler;
        vertex_data_store<VertexDataType> * vertex_data_handler;
        
        /* Computational context */
        ddgraph_context chicontext;
        
        /* Scheduler */
        bitset_scheduler * scheduler;
        
        /* Configuration */
        bool modifies_outedges;
        bool modifies_inedges;
        bool disable_outedges;
        bool only_adjacency;
        bool use_selective_scheduling;
        bool enable_deterministic_parallelism;
        bool store_inedges;
        bool disable_vertexdata_storage;
        bool preload_commit; //alow storing of modified edge data on preloaded data into memory
        bool randomization;
        bool initialize_edges_before_run;
        
        size_t blocksize;
        int membudget_mb;
        int load_threads;
        int exec_threads;
        
        /* State */
        vid_t sub_interval_st;
        vid_t sub_interval_en;
        int iter;
        int niters;
        size_t nupdates;
        size_t nedges;
        size_t work; // work is the number of edges processed
        unsigned int maxwindow;
        mutex modification_lock;
        
        bool reset_vertexdata;
        bool save_edgesfiles_after_inmemmode;
        
        /* Outputs */
        std::vector<ioutput<VertexDataType, EdgeDataType> *> outputs;
        
        /* Metrics */
        metrics &m;
        
        void print_config() {
            logstream(LOG_INFO) << "Engine configuration: " << std::endl;
            logstream(LOG_INFO) << " exec_threads = " << exec_threads << std::endl;
            logstream(LOG_INFO) << " load_threads = " << load_threads << std::endl;
            logstream(LOG_INFO) << " membudget_mb = " << membudget_mb << std::endl;
            logstream(LOG_INFO) << " blocksize = " << blocksize << std::endl;
            logstream(LOG_INFO) << " scheduler = " << use_selective_scheduling << std::endl;
        }
        
    public:
        
        /**
         * Initialize GraphChi engine
         * @param base_filename prefix of the graph files
         * @param nshards number of shards
         * @param selective_scheduling if true, uses selective scheduling 
         */
        ddgraph_engine(std::string _base_filename, int _nshards, bool _selective_scheduling, metrics &_m) : base_filename(_base_filename), nshards(_nshards), use_selective_scheduling(_selective_scheduling), m(_m) {
	    /*Modified by YouLi Cheng Start a server for this worker. */
            type_size = sizeof(EdgeDataType);
            inbuf = (char **) calloc(1024, sizeof(char*));
            obuf = (char **) calloc(1024, sizeof(char*));
            for(int i = 0; i<1024; i++){
                olength[i] = -1; 
                inlength[i] = -1;
                inbuf[i] = NULL;
                obuf[i] = NULL;
            }
            pthread_mutex_init(&tlock,NULL);  
            pthread_mutex_init(&olock,NULL);  
	    pthread_t ts1, ts2, ret;
	    ret = pthread_create(&ts1, NULL, &server1, NULL);
            assert(ret>=0);
	    ret = pthread_create(&ts2, NULL, &server2, NULL);
            assert(ret>=0);
            /* Initialize IO */
            m.start_time("iomgr_init");
            iomgr = new stripedio(m);
            if (disable_preloading()) {
                iomgr->set_disable_preloading(true);
            }
            m.stop_time("iomgr_init");
#ifndef DYNAMICEDATA
            logstream(LOG_INFO) << "Initializing ddgraph_engine. This engine expects " << sizeof(EdgeDataType)
            << "-byte edge data. " << std::endl;
#else
            logstream(LOG_INFO) << "Initializing ddgraph_engine with dynamic edge-data. This engine expects " << sizeof(int)
            << "-byte edge data. " << std::endl;

#endif
            /* If number of shards is unspecified - discover */
            if (nshards < 1) {
                nshards = get_option_int("nshards", 0);
                if (nshards < 1) {
                    logstream(LOG_WARNING) << "Number of shards was not specified (command-line argument 'nshards'). Trying to detect. " << std::endl;
                    nshards = discover_shard_num();
                }
            }
            
            /* Initialize a plenty of fields */
            memoryshard = NULL;
            modifies_outedges = true;
            modifies_inedges = true;
            save_edgesfiles_after_inmemmode = false;
            preload_commit = true;
            only_adjacency = false;
            disable_outedges = false;
            reset_vertexdata = false;
            initialize_edges_before_run = false;
            blocksize = get_option_long("blocksize", 4096 * 1024);
#ifndef DYNAMICEDATA
            while (blocksize % sizeof(EdgeDataType) != 0) blocksize++;
#endif
            
            disable_vertexdata_storage = false;

            membudget_mb = get_option_int("membudget_mb", 1024);
            nupdates = 0;
            iter = 0;
            work = 0;
            nedges = 0;
            scheduler = NULL;
            store_inedges = true;
            degree_handler = NULL;
            vertex_data_handler = NULL;
            enable_deterministic_parallelism = true;
            load_threads = get_option_int("loadthreads", 2);
            exec_threads = get_option_int("execthreads", omp_get_max_threads());
            maxwindow = 40000000;

            /* Load graph shard interval information */
            _load_vertex_intervals();
            
            _m.set("file", _base_filename);
            _m.set("engine", "default");
            _m.set("nshards", (size_t)nshards);
        }
        
        virtual ~ddgraph_engine() {
            if (degree_handler != NULL) delete degree_handler;
            if (vertex_data_handler != NULL) delete vertex_data_handler;
            if (memoryshard != NULL) {
                delete memoryshard;
                memoryshard = NULL;
            }
            for(int i=0; i < (int)sliding_shards.size(); i++) {
                if (sliding_shards[i] != NULL) {
                    delete sliding_shards[i];
                }
                sliding_shards[i] = NULL;
            }
            degree_handler = NULL;
            vertex_data_handler = NULL;
            delete iomgr;
        }
        
        
        
    protected:
        
        virtual degree_data * create_degree_handler() {
            return new degree_data(base_filename, iomgr);
        }
        
        virtual bool disable_preloading() {
            return false;
        }
        
        
            
        /**
         * Try to find suitable shards by trying with different
         * shard numbers. Looks up to shard number 2000.
         */
        int discover_shard_num() {
#ifndef DYNAMICEDATA
            int _nshards = find_shards<EdgeDataType>(base_filename);
#else
            int _nshards = find_shards<int>(base_filename);
#endif
            if (_nshards == 0) {
                logstream(LOG_ERROR) << "Could not find suitable shards - maybe you need to run sharder to create them?" << std::endl;
                logstream(LOG_ERROR) << "Was looking with filename [" << base_filename << "]" << std::endl;
                logstream(LOG_ERROR) << "You need to create the shards with edge data-type of size " << sizeof(EdgeDataType) << " bytes." << std::endl;
                logstream(LOG_ERROR) << "To specify the number of shards, use command-line parameter 'nshards'" << std::endl;
                assert(0);
            }
            return _nshards;
        }
        
       //to be rewrited. 
        virtual void initialize_sliding_shards() {
            assert(sliding_shards.size() == 0);
            for(int p=0; p < nshards; p++) {
#ifndef DYNAMICEDATA
                std::string edata_filename = filename_block_edata<EdgeDataType>(base_filename, exec_interval, p, P, 0);
                std::string adj_filename = filename_block_adj(base_filename, exec_interval, p, P);
                /* Let the IO manager know that we will be reading these files, and
                 it should decide whether to preload them or not.
                 */
                iomgr->allow_preloading(edata_filename);
                iomgr->allow_preloading(adj_filename);
#else
                std::string edata_filename = filename_shard_edata<int>(base_filename, p, nshards); //todo:
                std::string adj_filename = filename_block_adj(base_filename, exec_interval, p, P);
#endif
                
                
                sliding_shards.push_back(
                                         new slidingshard_t(iomgr, edata_filename, 
                                                            adj_filename,
                                                            intervals[p].first, 
                                                            intervals[p].second, 
                                                            blocksize, 
                                                            m, 
                                                            !modifies_outedges, 
                                                            only_adjacency));
                //if (!only_adjacency) 
                //    nedges += sliding_shards[sliding_shards.size() - 1]->num_edges();
            }
            
        }
        
        /* rewrited by YongLi Cheng. */
        virtual void initialize_sliding_shards_XG() {
            for(int p=0; p < nshards; p++) {
#ifndef DYNAMICEDATA
                std::string edata_filename = filename_block_edata<EdgeDataType>(bname, exec_interval, p, P, 1);
                std::string adj_filename = filename_block_adj(bname, exec_interval, p, P);
                /* Let the IO manager know that we will be reading these files, and
                 it should decide whether to preload them or not.
                 */
                iomgr->allow_preloading(edata_filename);
                iomgr->allow_preloading(adj_filename);
#else
                std::string edata_filename = filename_shard_edata<int>(base_filename, p, nshards); //todo:
                std::string adj_filename = filename_block_adj(base_filename, exec_interval, p, P);
#endif
                
                
                sliding_shards.push_back(
                      new slidingshard_t(iomgr, edata_filename, 
                            adj_filename,
                            intervals[p].first, 
                            intervals[p].second, 
                            blocksize, 
                            m, 
                            !modifies_outedges, 
                            only_adjacency));
            }
            
        }

        virtual void initialize_scheduler() {
            if (use_selective_scheduling) {
                if (scheduler != NULL) delete scheduler;
                scheduler = new bitset_scheduler((vid_t) num_vertices());
                scheduler->add_task_to_all();
            } else {
                scheduler = NULL;
            }
        }
        
        /**
         * If the data is only in one shard, we can just
         * keep running from memory.
         */
        bool is_inmemory_mode() {
            return nshards == 1;
        }
        
        
        /**
         * Extends the window to fill the memory budget, but not over maxvid
         */
        virtual vid_t determine_next_window(vid_t iinterval, vid_t fromvid, vid_t maxvid, size_t membudget) {
            /* Load degrees */
            degree_handler->load(fromvid, maxvid);
            
            /* If is in-memory-mode, memory budget is not considered. */
            if (is_inmemory_mode() || svertex_t().computational_edges()) {
                return maxvid;
            } else {
                size_t memreq = 0;
                int max_interval = maxvid - fromvid;
                for(int i=0; i < max_interval; i++) {
                    degree deg = degree_handler->get_degree(fromvid + i);
                    int inc = deg.indegree;
                    int outc = deg.outdegree * (!disable_outedges);
                    
                    // Raw data and object cost included
                    memreq += sizeof(svertex_t) + (sizeof(EdgeDataType) + sizeof(vid_t) + sizeof(ddgraph_edge<EdgeDataType>))*(outc + inc);
                    if (memreq > membudget) {
                        logstream(LOG_DEBUG) << "Memory budget exceeded with " << memreq << " bytes." << std::endl;
                        return fromvid + i - 1;  // Previous was enough
                    }
                }
                return maxvid;
            }
        }
        
        /** 
         * Calculates the exact number of edges
         * required to load in the subinterval.
         */
        size_t num_edges_subinterval(vid_t st, vid_t en) {
            size_t num_edges = 0;
            int nvertices = en - st + 1;
            if (scheduler != NULL) {
                for(int i=0; i < nvertices; i++) {
                    bool is_sched = scheduler->is_scheduled(st + i);
                    if (is_sched) {
                        degree d = degree_handler->get_degree(st + i);
                        num_edges += d.indegree * store_inedges + d.outdegree;
                    }
                }
            } else {
                for(int i=0; i < nvertices; i++) {
                    degree d = degree_handler->get_degree(st + i);
                    num_edges += d.indegree * store_inedges + d.outdegree;
                }
            }
            return num_edges;
        }
        
        virtual void load_before_updates(std::vector<svertex_t> &vertices) {
            volatile int done = 0;
            omp_set_num_threads(load_threads);
#pragma omp parallel for schedule(dynamic, 1)
            for(int p=-1; p < nshards; p++)  {
                if (p==(-1)) {
                    /* Load memory shard */
                    if (!memoryshard->loaded()) {
                        memoryshard->load_XG();
                    }
                    
                    /* Load vertex edges from memory shard */
                    memoryshard->load_vertices_XG(sub_interval_st, sub_interval_en, vertices, true, !disable_outedges);
                    
                    /* Load vertices */ 
                    vertex_data_handler->load(sub_interval_st, sub_interval_en);

                    /* Load vertices */
                    if (!disable_vertexdata_storage) {
                        vertex_data_handler->load(sub_interval_st, sub_interval_en);
                    }
                } else {
                    /* Load edges from a sliding shard */
                    if (!disable_outedges) {
                        if (p != exec_interval) {
                            if (randomization) {
                              sliding_shards[p]->set_disable_async_writes(true);   
                            }

                            sliding_shards[p]->read_next_vertices_XG(p,(int) vertices.size(), sub_interval_st, vertices,
                                                                  scheduler != NULL && chicontext.iteration == 0);
                            
                        }
                    }
                    __sync_add_and_fetch(&done, 1);
                }
            }
            
            /* Wait for all reads to complete */
            while(done < nshards) {}
            obuf_loaded = 1;
            iomgr->wait_for_reads();
        }
        
        void exec_updates(GraphChiProgram<VertexDataType, EdgeDataType, svertex_t> &userprogram,
                          std::vector<svertex_t> &vertices) {
            metrics_entry me = m.start_time();
            size_t nvertices = vertices.size();
            if (!enable_deterministic_parallelism) {
                for(int i=0; i < (int)nvertices; i++) vertices[i].parallel_safe = true;
            }
            int sub_interval_len = sub_interval_en - sub_interval_st;

            std::vector<vid_t> random_order(randomization ? sub_interval_len + 1 : 0);
            if (randomization) {
                // Randomize vertex-vector
                for(int idx=0; idx <= (int)sub_interval_len; idx++) random_order[idx] = idx;
                std::random_shuffle(random_order.begin(), random_order.end());
            }
             
            do {
                omp_set_num_threads(exec_threads);
                
        #pragma omp parallel sections 
                    {
        #pragma omp section
                        {
        #pragma omp parallel for schedule(dynamic)
                            for(int idx=0; idx <= (int)sub_interval_len; idx++) {
                                vid_t vid = sub_interval_st + (randomization ? random_order[idx] : idx);
                                svertex_t & v = vertices[vid - sub_interval_st];
                                
                                if (exec_threads == 1 || v.parallel_safe) {
                                    if (!disable_vertexdata_storage)
                                        v.dataptr = vertex_data_handler->vertex_data_ptr(vid);
                                    if (v.scheduled){ 
                                        userprogram.update(v, chicontext);
				    }
                                }
                            }
                        }
        #pragma omp section
                        {
                            if (exec_threads > 1 && enable_deterministic_parallelism) {
                                int nonsafe_count = 0;
                                for(int idx=0; idx <= (int)sub_interval_len; idx++) {
                                    vid_t vid = sub_interval_st + (randomization ? random_order[idx] : idx);
                                    svertex_t & v = vertices[vid - sub_interval_st];
                                    if (!v.parallel_safe && v.scheduled) {
                                        if (!disable_vertexdata_storage)
                                            v.dataptr = vertex_data_handler->vertex_data_ptr(vid);
                                        userprogram.update(v, chicontext);
                                        nonsafe_count++;
                                    }
                                }
                                
                                m.add("serialized-updates", nonsafe_count);
                            }
                        }
                }
            } while (userprogram.repeat_updates(chicontext));
            m.stop_time(me, "execute-updates");
        }
        

        void prev_updates(GraphChiProgram<VertexDataType, EdgeDataType, svertex_t> &userprogram,
                          std::vector<svertex_t> &vertices) {
            size_t nvertices = vertices.size();
            if (!enable_deterministic_parallelism) {
                for(int i=0; i < (int)nvertices; i++) vertices[i].parallel_safe = true;
            }
            int sub_interval_len = sub_interval_en - sub_interval_st;

            std::vector<vid_t> random_order(randomization ? sub_interval_len + 1 : 0);
            if (randomization) {
                // Randomize vertex-vector
                for(int idx=0; idx <= (int)sub_interval_len; idx++) random_order[idx] = idx;
                std::random_shuffle(random_order.begin(), random_order.end());
            }
             
            do {
                omp_set_num_threads(exec_threads);
                
        #pragma omp parallel sections 
                    {
        #pragma omp section
                        {
        #pragma omp parallel for schedule(dynamic)
                            for(int idx=0; idx <= (int)sub_interval_len; idx++) {
                                vid_t vid = sub_interval_st + (randomization ? random_order[idx] : idx);
                                svertex_t & v = vertices[vid - sub_interval_st];


if(v.crucial_flag == false){
                                
                                if (exec_threads == 1 || v.parallel_safe) {
                                    if (!disable_vertexdata_storage)
                                        v.dataptr = vertex_data_handler->vertex_data_ptr(vid);
                                    if (v.scheduled){ 
                                        userprogram.update(v, chicontext);
				                    }
                                }

}                                


                            }
                        }
        #pragma omp section
                        {
                            if (exec_threads > 1 && enable_deterministic_parallelism) {
                                int nonsafe_count = 0;
                                for(int idx=0; idx <= (int)sub_interval_len; idx++) {
                                    vid_t vid = sub_interval_st + (randomization ? random_order[idx] : idx);
                                    svertex_t & v = vertices[vid - sub_interval_st];
if(v.crucial_flag==false){

                                    if (!v.parallel_safe && v.scheduled) {
                                        if (!disable_vertexdata_storage)
                                            v.dataptr = vertex_data_handler->vertex_data_ptr(vid);
                                        userprogram.update(v, chicontext);
                                        nonsafe_count++;
                                    }

}


                                }
                                
                                m.add("serialized-updates", nonsafe_count);
                            }
                        }
                }
            } while (userprogram.repeat_updates(chicontext));
        }
        
        void crucial_updates(GraphChiProgram<VertexDataType, EdgeDataType, svertex_t> &userprogram,
                          std::vector<svertex_t> &vertices) {
            metrics_entry me = m.start_time();
            size_t nvertices = vertices.size();
            if (!enable_deterministic_parallelism) {
                for(int i=0; i < (int)nvertices; i++) vertices[i].parallel_safe = true;
            }
            int sub_interval_len = sub_interval_en - sub_interval_st;

            std::vector<vid_t> random_order(randomization ? sub_interval_len + 1 : 0);
            if (randomization) {
                // Randomize vertex-vector
                for(int idx=0; idx <= (int)sub_interval_len; idx++) random_order[idx] = idx;
                std::random_shuffle(random_order.begin(), random_order.end());
            }
             
            do {
                omp_set_num_threads(exec_threads);
                
        #pragma omp parallel sections 
                    {
        #pragma omp section
                        {
        #pragma omp parallel for schedule(dynamic)
                            for(int idx=0; idx <= (int)sub_interval_len; idx++) {
                                vid_t vid = sub_interval_st + (randomization ? random_order[idx] : idx);
                                svertex_t & v = vertices[vid - sub_interval_st];


if(v.crucial_flag == true){
                                
                                if (exec_threads == 1 || v.parallel_safe) {
                                    if (!disable_vertexdata_storage)
                                        v.dataptr = vertex_data_handler->vertex_data_ptr(vid);
                                    if (v.scheduled){ 
                                        userprogram.update(v, chicontext);
				                    }
                                }

}                                


                            }
                        }
        #pragma omp section
                        {
                            if (exec_threads > 1 && enable_deterministic_parallelism) {
                                int nonsafe_count = 0;
                                for(int idx=0; idx <= (int)sub_interval_len; idx++) {
                                    vid_t vid = sub_interval_st + (randomization ? random_order[idx] : idx);
                                    svertex_t & v = vertices[vid - sub_interval_st];
if(v.crucial_flag==true){

                                    if (!v.parallel_safe && v.scheduled) {
                                        if (!disable_vertexdata_storage)
                                            v.dataptr = vertex_data_handler->vertex_data_ptr(vid);
                                        userprogram.update(v, chicontext);
                                        nonsafe_count++;
                                    }

}


                                }
                                
                                m.add("serialized-updates", nonsafe_count);
                            }
                        }
                }
            } while (userprogram.repeat_updates(chicontext));
            m.stop_time(me, "execute-updates");
        }
        


        /**
         Special method for running all iterations with the same vertex-vector.
         This is a hacky solution.

         FIXME:  this does not work well with deterministic parallelism. Needs a
         a separate analysis phase to check which vertices can be run in parallel, and
         then run it in chunks. Not difficult.
         **/
        void exec_updates_inmemory_mode(GraphChiProgram<VertexDataType, EdgeDataType, svertex_t> &userprogram,
                                        std::vector<svertex_t> &vertices) {
            work = nupdates = 0;
            for(iter=0; iter<niters; iter++) {
                logstream(LOG_INFO) << "In-memory mode: Iteration " << iter << " starts. (" << chicontext.runtime() << " secs)" << std::endl;
                chicontext.iteration = iter;
                if (iter > 0) // First one run before -- ugly
                    userprogram.before_iteration(iter, chicontext);
                userprogram.before_exec_interval(0, (int)num_vertices(), chicontext);
                
                if (use_selective_scheduling) {
                    scheduler->new_iteration(iter);
                    if (iter > 0 && !scheduler->has_new_tasks) {
                        logstream(LOG_INFO) << "No new tasks to run!" << std::endl;
                        break;
                    }
                    for(int i=0; i < (int)vertices.size(); i++) { // Could, should parallelize
                        if (iter == 0 || scheduler->is_scheduled(i)) {
                            vertices[i].scheduled =  true;
                            nupdates++;
                            work += vertices[i].inc + vertices[i].outc;
                        } else {
                            vertices[i].scheduled = false;
                        }
                    }
                    
                    scheduler->has_new_tasks = false; // Kind of misleading since scheduler may still have tasks - but no new tasks.
                } else {
                    nupdates += num_vertices();
                    //work += num_edges();
                }
                
                exec_updates(userprogram, vertices);
                load_after_updates(vertices);
                
                userprogram.after_exec_interval(0, (int)num_vertices(), chicontext);
                userprogram.after_iteration(iter, chicontext);
                if (chicontext.last_iteration > 0 && chicontext.last_iteration <= iter){
                   logstream(LOG_INFO)<<"Stopping engine since last iteration was set to: " << chicontext.last_iteration << std::endl;
                   break;
                }

            }
            
            if (save_edgesfiles_after_inmemmode) {
                logstream(LOG_INFO) << "Saving memory shard..." << std::endl;
                
            }
        }
        

        virtual void init_vertices(std::vector<svertex_t> &vertices, ddgraph_edge<EdgeDataType> * &edata) {
            size_t nvertices = vertices.size();
            
            /* Compute number of edges */
            size_t num_edges = num_edges_subinterval(sub_interval_st, sub_interval_en);
            
            /* Allocate edge buffer */
            edata = (ddgraph_edge<EdgeDataType>*) malloc(num_edges * sizeof(ddgraph_edge<EdgeDataType>));
            
            /* Assign vertex edge array pointers */
            size_t ecounter = 0;
            for(int i=0; i < (int)nvertices; i++) {
                degree d = degree_handler->get_degree(sub_interval_st + i);
                int inc = d.indegree;
                int outc = d.outdegree * (!disable_outedges);
                vertices[i] = svertex_t(sub_interval_st + i, &edata[ecounter], 
                                        &edata[ecounter + inc * store_inedges], inc, outc);
                if (scheduler != NULL) {
                    bool is_sched = ( scheduler->is_scheduled(sub_interval_st + i));
                    if (is_sched) {
                        vertices[i].scheduled =  true;
                        nupdates++;
                        ecounter += inc * store_inedges + outc;
                    }
                } else {
                    nupdates++; 
                    vertices[i].scheduled =  true;
                    ecounter += inc * store_inedges + outc;               
                }
            }                   
            work += ecounter;
            assert(ecounter <= num_edges);
        }
        
        
        void save_vertices(std::vector<svertex_t> &vertices) {
            if (disable_vertexdata_storage) return;
            size_t nvertices = vertices.size();
            bool modified_any_vertex = false;
            for(int i=0; i < (int)nvertices; i++) {
                if (vertices[i].modified) {
                    modified_any_vertex = true;
                    break;
                }
            }
            if (modified_any_vertex) {
                vertex_data_handler->save();
            }
        }
        
        virtual void load_after_updates(std::vector<svertex_t> &vertices) {
            // Do nothing.
        }   
        
        virtual void write_delta_log() {
            // Write delta log
            std::string deltafname = iomgr->multiplexprefix(0) + base_filename + ".deltalog";
            FILE * df = fopen(deltafname.c_str(), (chicontext.iteration == 0  ? "w" : "a"));
            fprintf(df, "%d,%lu,%lu,%lf\n", chicontext.iteration, nupdates, work, chicontext.get_delta()); 
            fclose(df);
        }
        
    public:
        
        virtual std::vector< std::pair<vid_t, vid_t> > get_intervals() {
            return intervals;
        }
        
        virtual std::pair<vid_t, vid_t> get_interval(int i) {
            return intervals[i];
        }
        
        /**
         * Returns first vertex of i'th interval.
         */
        vid_t get_interval_start(int i) {
            return get_interval(i).first;
        }
        
        /** 
         * Returns last vertex (inclusive) of i'th interval.
         */
        vid_t get_interval_end(int i) {
            return get_interval(i).second;
        }
        
        virtual size_t num_vertices() {
            return 1 + intervals[nshards - 1].second;
        }
        
        ddgraph_context &get_context() {
            return chicontext;
        }
        
        virtual int get_nshards() {
            return nshards;
        }
        
        size_t num_updates() {
            return nupdates;
        }
        
        /**
         * Thread-safe version of num_edges
         */
        virtual size_t num_edges_safe() {
            return num_edges();
        }
        
        virtual size_t num_buffered_edges() {
            return 0;
        }
        
        /** 
         * Counts the number of edges from shard sizes.
         */
        virtual size_t num_edges() {
            if (sliding_shards.size() == 0) {
                logstream(LOG_ERROR) << "engine.num_edges() can be called only after engine has been started. To be fixed later. As a workaround, put the engine into a global variable, and query the number afterwards in begin_iteration(), for example." << std::endl;
                assert(false);
            }
            if (only_adjacency) {
                // TODO: fix.
                logstream(LOG_ERROR) << "Asked number of edges, but engine was run without edge-data." << std::endl; 
                return 0;
            }
            return nedges;
        }
        
        /**
         * Checks whether any vertex is scheduled in the given interval.
         * If no scheduler is configured, returns always true.
         */
        // TODO: support for a minimum fraction of scheduled vertices
        bool is_any_vertex_scheduled(vid_t st, vid_t en) {
            if (scheduler == NULL) return true;
            for(vid_t v=st; v<=en; v++) {
                if (scheduler->is_scheduled(v)) {
                    return true;
                }
            }
            return false;
        }
        
        virtual void initialize_iter() {
            // Do nothing
        }
        
        virtual void initialize_before_run() {
            if (reset_vertexdata) {
                vertex_data_handler->clear(num_vertices());
            }
        }
        
        virtual memshard_t * create_memshard(vid_t interval_st, vid_t interval_en) {
#ifndef DYNAMICEDATA
            return new memshard_t(this->iomgr,
                                  filename_shard_edata<EdgeDataType>(base_filename, exec_interval, nshards),  
                                  filename_shard_adj(base_filename, exec_interval, nshards),  
                                  interval_st, 
                                  interval_en,
                                  blocksize,
                                  m);
#else
            return new memshard_t(this->iomgr,
                                  filename_shard_edata<int>(base_filename, exec_interval, nshards),
                                  filename_shard_adj(base_filename, exec_interval, nshards),
                                  interval_st,
                                  interval_en,
                                  blocksize,
                                  m);
#endif
        }
        
	/*
	 * Rewrited by Yongli Cheng 2014/04/10
	 */
	void run(GraphChiProgram<VertexDataType, EdgeDataType, svertex_t> &userprogram, int _niters) {
            m.start_time("runtime");
	        randomization = get_option_int("randomization", 0) == 1;
            
            std::cout<<"running...."<<std::endl;
            int fd_rw,c,m1;
	        int wret, len, ntop;
	        int c1=0;
	        char buf[1024 * 1024 + 1];
	        char workerIP[1024][15]; //save worker IP 
            int ret;
            while(fd_master < 0) {}//waiting for xgMaster connection. volatile int fd_master
            ret =  read(fd_master, buf, 1024 * 1024); //Receive 'S' Package
	        if(ret < 0) std::cout<<"xgMaster shutdown, exit...."<<std::endl;
            assert(ret > 0);

            memcpy(&c, buf+6, sizeof(int));
            memcpy(&M, buf+6+sizeof(int), sizeof(int));
            memcpy(&N, buf+6+2*sizeof(int), sizeof(int));
            memcpy(&P, buf+6+3*sizeof(int), sizeof(int));
            memcpy(&CO, buf+6+4*sizeof(int), sizeof(int));
            memcpy(&ntop, buf+6+5*sizeof(int), sizeof(int));

            T = (task*)malloc(sizeof(task) * c);
            iL = (intervallock*)malloc(sizeof(intervallock) * P);

	        //std::cout<<"C:"<<c<<"M: "<<M<<"N:  "<<N<<"P:  "<<P<<"  CO: "<<CO<<" ntop: "<<ntopstd::endl;

	        /* worker IP */
	        for(int k = 0; k < M; k++) memcpy(workerIP[k], buf+6+6*sizeof(int)+k*15, 15);
	        for(int k = 0; k < M; k++) std::cout<<"worker: "<<workerIP[k]<<std::endl;

	        /* Init iL */
            for(int i=0; i<P; i++){
		        iL[i].interval = -1;
		        for(int j=0;j<2*P-1;j++) iL[i].lock[j] = '0';
            }

	        /* Init interval_dirt_disk and interval_dirt_memory and interval_dirt_main */
	        for(int i = 0; i < P; i ++){
		        if(CO == 0){	//need not to conctrol
			        interval_dirt_disk[i] = 0;
			        interval_dirt_memory[i] = 0;
			        interval_dirt_main[i] = 0;
		        }else if(CO == 1 || CO == 2) {
			        interval_dirt_memory[i] = mini(M, i);
			        interval_dirt_disk[i] = maxi(0, i - M);
			        interval_dirt_main[i] = 0;
		        }else if(CO ==3){
			        interval_dirt_memory[i] = mini(M, i) * 2;
			        interval_dirt_disk[i] = maxi(0, i - M) * 2;
			        interval_dirt_main[i] = 0;

		        }		
	        }

            for(int j = 0; j < c; j++){
                T[j]= *(task*)(buf+6+6*sizeof(int)+15*M+j*sizeof(task));
                //std::cout<<T[j].machine<<"==== " <<T[j].iter<<" "<<T[j].interval<<std::endl;
		        if(T[j].iter == 0) {
			    iL[T[j].interval].interval = T[j].interval;
                   	//std::cout<<T[j].machine<<"==== " <<T[j].iter<<" "<<T[j].interval<<std::endl;
		        } 
            }

	        for(int i=0; i<P; i++){
		        //std::cout<<iL[i].interval<<std::endl;
		        //for(int j=0;j<2*P-1;j++) std::cout<<"  "<<iL[i].lock[j];
		       // std::cout<<std::endl;
	        }

            /* send a message to xgMaster. */
            char* msg = (char*)"XGACKS";
            if (send(fd_master, msg, strlen(msg), 0) == -1)
                   std::cout<<"send back XGACKS fail!"<<std::endl;

            if (degree_handler == NULL)
                degree_handler = create_degree_handler();

            niters = N;   //todo:need or not
            logstream(LOG_INFO) << "GX starting" << std::endl;
            logstream(LOG_INFO) << "Licensed under the Apache License 2.0" << std::endl;
            logstream(LOG_INFO) << "Copyright YongLi Cheng et al., HuaZhong Technology University (2014)" << std::endl;

            if (vertex_data_handler == NULL)
                vertex_data_handler = new vertex_data_store<VertexDataType>(base_filename, num_vertices(), iomgr);

            initialize_before_run();
            initialize_scheduler();
            /* Init slidingshards. */
            bname = base_filename;
            initialize_sliding_shards_XG();
            omp_set_nested(1);

            /* Install a 'mock'-scheduler to chicontext if scheduler
             is not used. */
            chicontext.scheduler = scheduler;
            if (scheduler == NULL) {
                chicontext.scheduler = new non_scheduler();
            }

            /* Print configuration */
            print_config();

            /* Main loop */
	        int iter_prev = -1;
            int ll_in_to_out = 0, ll_out_to_in = 0;//check data length
            std::string block_filename;
            for(int i=0; i < c; i++) {
		        int iter = T[i].iter;
		        interval = T[i].interval;
                exec_interval = interval % P;
                if(i == 0) exec_interval_buf = exec_interval;
                logstream(LOG_INFO) << "Start iteration: " << iter << "Interval: " << interval << std::endl;

                metrics_entry me2 = m.start_time();
                std::stringstream ss2;
                ss2<<"prepare time ";

		        if(iter > iter_prev) { //New iteration
			        iter_prev ++ ;
			        initialize_iter(); //do nothing.

			        /* Check vertex data file has the right size (number of vertices may change) */
			        if (!disable_vertexdata_storage)
				        vertex_data_handler->check_size(num_vertices());

			        /* Keep the context object updated */
			        chicontext.filename = base_filename;//to be a global on xgMaster
			        chicontext.iteration = iter;
			        chicontext.num_iterations = niters;
			        chicontext.nvertices = num_vertices();
			        //if (!only_adjacency) chicontext.nedges = num_edges();

			        chicontext.execthreads = exec_threads;
			        chicontext.reset_deltas(exec_threads);


			        /* Call iteration-begin event handler */
			        userprogram.before_iteration(iter, chicontext);

			        /* Check scheduler. If no scheduled tasks, terminate. */
			        if (use_selective_scheduling) {
				        if (scheduler != NULL) {
					        if (!scheduler->has_new_tasks) {
						        logstream(LOG_INFO) << "No new tasks to run!" << std::endl;
					        }
					        scheduler->has_new_tasks = false; // Kind of misleading since scheduler may still have tasks - but no new tasks.
				        }
			        }

			        /* Now clear scheduler bits for the interval */
			        if (scheduler != NULL)
				        scheduler->new_iteration(iter);
		        }

                /* Determine interval limits */
                vid_t interval_st = get_interval_start(exec_interval);
                vid_t interval_en = get_interval_end(exec_interval);
                                        

                /* for previous update */
                if(prev_bitset != NULL) delete prev_bitset;
                prev_bitset = new dense_bitset(interval_en - interval_st + 1);
                prev_bitset->setall();

                chicontext.range_st = interval_st;
                chicontext.range_en = interval_en;

                if (interval_st > interval_en) continue; // Can happen on very very small graphs.

                if (!is_inmemory_mode())
                    userprogram.before_exec_interval(interval_st, interval_en, chicontext);

                /* Initialize memory shard */
                if (memoryshard != NULL) delete memoryshard;
                memoryshard = create_memshard(interval_st, interval_en);
                memoryshard->only_adjacency = only_adjacency;
                memoryshard->set_disable_async_writes(randomization);

                sub_interval_st = interval_st;
		        sub_interval_en = interval_en;

/* Loading .....start     */
		        while((interval_dirt_disk[exec_interval] > 0 || interval_dirt_main[exec_interval] == 1) && ( CO != 0)) {
			        //std::cout<<"interval_dirt_disk "<<exec_interval<<"("<<interval<<")"<<" : "<<interval_dirt_disk[exec_interval]<<std::endl;
			        //std::cout<<"interval_dirt_main "<<exec_interval<<"("<<interval<<")"<<" : "<<interval_dirt_main[exec_interval]<<std::endl;
		        }
                if(CO == 0){
			        //do nothing;
                }else if(CO == 1 || CO == 2) {
                   interval_dirt_disk[exec_interval] = maxi(P - 1 - M, 0);
                   interval_dirt_main[exec_interval] = 1;
                }else if(CO ==3){
                   interval_dirt_disk[exec_interval] = maxi((P - 1 - M) * 2, 0);
                   interval_dirt_main[exec_interval] = 1;

                }
/* Loading .....end     */

		        degree_handler->load(sub_interval_st, sub_interval_en);

                logstream(LOG_INFO) << chicontext.runtime() << "s: Starting: "
                << sub_interval_st << " -- " << interval_en << std::endl;

                 /* Initialize vertices */
                 modification_lock.lock();
                 int nvertices = sub_interval_en - sub_interval_st + 1;
                 ddgraph_edge<EdgeDataType> * edata = NULL;
                 std::vector<svertex_t> vertices(nvertices, svertex_t());
                 init_vertices(vertices, edata);

                /* block from interval - M block */

                int pos = exec_interval;
                for(int j = 1; j <= M; j++){
                    pos --;
                    if(pos < 0) pos = P - 1;
                }
                if((CO == 1) && interval >= M){
                    if(inbuf[pos] != NULL) free(inbuf[pos]);
                    inbuf[pos] = (char*)malloc(ll_out_to_in+1);
                    memcpy(inbuf[pos], out_to_in_buf, ll_out_to_in);
                    inlength[pos] = ll_out_to_in;
                    if(out_to_in_buf != NULL) free(out_to_in_buf);
                    out_to_in_buf = NULL;
                    //std::cout<<"interval_dirt_memory"<<exec_interval<<"  "<<interval_dirt_memory[exec_interval]<<std::endl;
                    __sync_add_and_fetch(&interval_dirt_memory[exec_interval],-1);
                }else if((CO ==2)&& interval>= M){
                    if(obuf[pos]!=NULL) free(obuf[pos]);
                    obuf[pos] = (char*)malloc(ll_in_to_out+1);
                    memcpy(obuf[pos], in_to_out_buf, ll_in_to_out);
                    olength[pos] = ll_in_to_out;
                    if(in_to_out_buf != NULL) free(in_to_out_buf);
                    in_to_out_buf = NULL;
                    //std::cout<<"interval_dirt_memory"<<exec_interval<<"  "<<interval_dirt_memory[exec_interval]<<std::endl;
                    __sync_add_and_fetch(&interval_dirt_memory[exec_interval],-1);
                }else if(CO==3 && interval>=M){
                    if(inbuf[pos] != NULL) free(inbuf[pos]);
                    inbuf[pos] = (char*)malloc(ll_out_to_in+1);
                    memcpy(inbuf[pos], out_to_in_buf, ll_out_to_in);
                    inlength[pos] = ll_out_to_in;
                    if(out_to_in_buf != NULL) free(out_to_in_buf);
                    out_to_in_buf = NULL;

                    if(obuf[pos]!=NULL) free(obuf[pos]);
                    obuf[pos] = (char*)malloc(ll_in_to_out+1);
                    memcpy(obuf[pos], in_to_out_buf, ll_in_to_out);
                    olength[pos] = ll_in_to_out;
                    if(in_to_out_buf != NULL) free(in_to_out_buf);
                    in_to_out_buf = NULL;

                    __sync_add_and_fetch(&interval_dirt_memory[exec_interval],-2);
                }

                  /* Load data */
                 load_before_updates(vertices);

                 modification_lock.unlock();


                /* previous update */
                std::cout<<"exec:     prev_update"<<std::endl;
                metrics_entry me3 = m.start_time();
                std::stringstream ss3;
                ss3<<"prev_update time ";
                if(CO>0) prev_updates(userprogram, vertices);
                m.stop_time(me3,ss3.str(),false);
                 
                if(CO>0) logstream(LOG_INFO) << "Waiting for crucial blocks......" << std::endl;


                metrics_entry me5 = m.start_time();
                std::stringstream ss5;
                ss5<<"wait for crucial_block time ";
		        while((interval_dirt_memory[exec_interval] > 0) && (CO != 0)){
			        //std::cout<<"interval_dirt_memory "<<interval<<" : "<<interval_dirt_memory[exec_interval]<<std::endl;
		        } // Waiting for memory sliding block.
                m.stop_time(me5,ss5.str(),false);
                inbuf_loaded = 0;
                obuf_loaded = 0;
                if(CO == 0){
			        //do nothing
                }else if(CO == 1 || CO == 2) {
                        interval_dirt_memory[exec_interval] = mini(M, P-1);
                }else if(CO ==3){
                        interval_dirt_memory[exec_interval] = mini(M * 2, (P - 1)*2);
                }
		        /* Waiting for signal to update end. */	

                /* Waiting  signal from xgMaster. */
		        int ret;
                if(((T[i].interval != 0)&&(CO != 0))||((CO == 0)&&(M<P))) { //Waiting for xgMaster.
                    while((ret = read(fd_master, buf, 6)) < 6){
			            std::cout<<"Waiting for xgMaster signal to update..." <<std::endl;
                        if(ret == -1){
                                std::cout<<"Recv fail, exit..."<<std::endl;
                                exit(0);
                        }
		            }
                }


                /* Execute updates */
                metrics_entry me1 = m.start_time();
                std::stringstream ss1;
                ss1<<"crucial_updatetime ";
                logstream(LOG_INFO) << "Start updates" << std::endl;
                 if(CO == 0 && M == P){
                     for(int y = 0; y < c; y++){
                        ret = read(fd_master, buf, 6);
                        exec_updates(userprogram, vertices);
                        if (!disable_vertexdata_storage) {
                            save_vertices(vertices);
                        }
                       ret = send(fd_master, "XGREQU", 6, 0); 
                     }
                     break;
                 } else{
                      //exec_updates(userprogram, vertices);
                      if(CO>0) {
                          std::cout<<"exec:     crucial_update"<<std::endl;
                          crucial_updates(userprogram, vertices);
                      }
                      else if(CO==0) exec_updates(userprogram, vertices);
                 }
                 m.stop_time(me1,ss1.str(),false);
                /* Load phase after updates (used by the functional engine) */
                load_after_updates(vertices);
                logstream(LOG_INFO) << "Finished updates" << std::endl;
                if(i+1 < c) exec_interval_buf = T[i+1].interval % P;

		        if(CO == 0){
                      while((ret = send(fd_master, "XGREQU", 6, 0)) < 6){
                             std::cout<<"Send Successful signal to xgMaster..."<<std::endl;
                             if(ret == -1){
                                    std::cout<<"Recv fail, exit..."<<std::endl;
                                    exit(0);
                             }
                      }//this interval have been updated
		        }

                metrics_entry me6 = m.start_time();
                std::stringstream ss6;
                ss6<<"send blocks time ";
#ifndef DYNAMICEDATA
                /* release obuf */
                if(CO == 2 || CO == 0){
                    for(int n = 0; n < P; n++){
                        if(obuf[n] != NULL) free(obuf[n]);
                        olength[n] = -1;
                        obuf[n] = NULL;
                    }
                }

                /* release inbuf */
                if(CO == 1 || CO == 0){
                    //only commit main shard.
                    if(CO == 1) {
                        if(main_to_main_buf!=NULL) free(main_to_main_buf);
                        main_to_main_buf = (char*)malloc(inlength[exec_interval]);
                        memcpy(main_to_main_buf, inbuf[exec_interval], inlength[exec_interval]);
                        ll_main_to_main = inlength[exec_interval];
                    }
                   // writeedatablock(base_filename,(char*)inbuf[exec_interval], exec_interval, exec_interval , P, 0,inlength[exec_interval]);
                    for(int m = 0; m < P; m++){
                        if(inbuf[m] != NULL) free(inbuf[m]);
                        inbuf[m] = NULL;
                        inlength[m] = -1;
                    }
                }

                if(CO == 3) {
                    if(main_to_main_buf!=NULL) free(main_to_main_buf);
                    main_to_main_buf = (char*)malloc(inlength[exec_interval]);
                    memcpy(main_to_main_buf, inbuf[exec_interval], inlength[exec_interval]);
                    ll_main_to_main = inlength[exec_interval];
                    //writeedatablock(base_filename,(char*)inbuf[exec_interval], exec_interval, exec_interval , P, 0,inlength[exec_interval]);
                    if(inbuf[exec_interval] != NULL) free(inbuf[exec_interval]);
                    inbuf[exec_interval] = NULL;
                    inlength[exec_interval] = -1;
                }

		        /* send Related datat to other intervls */
                int ssize;
		        if(CO != 0){
		            for(int j = 1; j < P; j++){
                       int from, to;
                       from = T[i].interval;
                       to   = T[i].interval + j;
			            if(j == 1){//inform next interval to begin.
                	    	/* send to xgMaster to next */
                		    while((ret = send(fd_master, "XGREQU", 6, 0)) < 6){
                        		std::cout<<"Send Successful signal to xgMaster..."<<std::endl;
                        		if(ret == -1){
                              			std::cout<<"Recv fail, exit..."<<std::endl;
                              			exit(0);
                        		}
                		    }//this interval have been updated
			            }
                        if((to - from) % M == 0){//local
                           if( j == M ){
                                //in-memory, need not saving to disk
                                if(CO == 0){
                                    //do nothing
                                }else if(CO == 1) {
						            if(to > P * N -1){
                                        //check data length
                                        /*
                                        block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 0);
                                        if (file_exists(block_filename)){
                                            size_t fsize = get_filesize(block_filename);
                                            if(fsize != olength[to%P]){
                                                std::cout<<"edata is bad!"<<std::endl;
                                            }
                                            assert(fsize==olength[to%P]);
                                        }else{
                                            std::cout<<"can not find file: " << block_filename <<std::endl;
                                            assert(false);
                                        }
                                        */
                                        __sync_add_and_fetch(&cout_g, 1);
                                        writeedatablock(base_filename,(char*)obuf[to % P], to%P, from%P , P, 0,olength[to%P]);
                                    }else{
                                        out_to_in_buf = (char*)malloc(olength[to%P] + 1);
                                        ll_out_to_in = olength[to%P];
                                        __sync_add_and_fetch(&cout_g, 1);
                                        memcpy(out_to_in_buf, obuf[to % P], olength[to%P]);
                                    }
                               }else if(CO == 2) {
						            if(to > P * N -1){
                                        //check data length
                                        /*
                                        block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 1);
                                        if (file_exists(block_filename)){
                                            size_t fsize = get_filesize(block_filename);
                                            if(fsize != inlength[to%P]){
                                                std::cout<<"edata is bad!"<<std::endl;
                                            }
                                            assert(fsize==inlength[to%P]);
                                        }else{
                                            std::cout<<"can not find file: " << block_filename <<std::endl;
                                            assert(false);
                                        }
                                        */
                                        writeedatablock(base_filename, (char*)inbuf[to%P], to%P, from%P , P, 1,inlength[to%P]);
                                    }else{
                                        if(in_to_out_buf != NULL) free(in_to_out_buf);
                                        in_to_out_buf = (char*)malloc(inlength[to%P]+1);
                                        ll_in_to_out = inlength[to%P];
                                        memcpy(in_to_out_buf, inbuf[to % P],inlength[to%P]);
                                    }
                               }else if(CO ==3){
						            if(to > P * N -1){
                                        //check data length
                                        /*
                                        block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 0);
                                        if (file_exists(block_filename)){
                                            size_t fsize = get_filesize(block_filename);
                                            if(fsize != olength[to%P]){
                                                std::cout<<"edata is bad!"<<std::endl;
                                            }
                                            assert(fsize==olength[to%P]);
                                        }else{
                                            std::cout<<"can not find file: " << block_filename <<std::endl;
                                            assert(false);
                                        }
                                        */
                                        writeedatablock(base_filename,(char*)obuf[to % P], to%P, from%P , P, 0,olength[to%P]);
                                        //check data length
                                        /*
                                        block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 1);
                                        if (file_exists(block_filename)){
                                            size_t fsize = get_filesize(block_filename);
                                            if(fsize != inlength[to%P]){
                                                std::cout<<"edata is bad!"<<std::endl;
                                            }
                                            assert(fsize==inlength[to%P]);
                                        }else{
                                            std::cout<<"can not find file: " << block_filename <<std::endl;
                                            assert(false);
                                        }
                                        */
                                        writeedatablock(base_filename,(char*)inbuf[to%P], to%P, from%P , P, 1,inlength[to%P]);
                                    }else{
                                        if(out_to_in_buf != NULL) free(out_to_in_buf);
                                        out_to_in_buf = (char*)malloc(olength[to%P] + 1);
                                        memcpy(out_to_in_buf, obuf[to % P], olength[to%P]);
                                        ll_out_to_in = olength[to%P];
                                        if(obuf[to%P] != NULL) free(obuf[to%P]);
                                        obuf[to%P] = NULL;
                                        olength[to%P] = -1;

                                        if(in_to_out_buf != NULL) free(in_to_out_buf);
                                        in_to_out_buf = (char*)malloc(inlength[to%P]+1);
                                        memcpy(in_to_out_buf, inbuf[to % P],inlength[to%P]);
                                        ll_in_to_out = inlength[to%P];
                                        if(inbuf[to%P] != NULL) free(inbuf[to%P]);
                                        inbuf[to%P] = NULL;
                                        inlength[to%P] = -1;

                                        //interval_dirt_memory[(T[i].interval + j) % P] -= 2;
                                        //__sync_add_and_fetch(&interval_dirt_memory[(T[i].interval + j) % P],-1);
                                    }
                                }
                           }else {
                                //save to disk
                                if(CO == 0){
                                     //do nothing.
                                }else if(CO == 1) {
                                   //check data length
                                   /*
                                   block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 0);
                                   if (file_exists(block_filename)){
                                        size_t fsize = get_filesize(block_filename);
                                        if(fsize != olength[to%P]){
                                           std::cout<<"edata is bad!"<<std::endl;
                                        }
                                        assert(fsize==olength[to%P]);
                                    }else{
                                        std::cout<<"can not find file: " << block_filename <<std::endl;
                                        assert(false);
                                    }
                                    */
                                     writeedatablock(base_filename, (char*)obuf[to % P], to%P, from%P , P, 0,olength[to%P]);
                                      __sync_add_and_fetch(&cout_g, 1);
                                     __sync_add_and_fetch(&interval_dirt_disk[(T[i].interval + j) % P],-1);
                                }else if(CO == 2){
                                     //check data length
                                     /*
                                     block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 1);
                                     if (file_exists(block_filename)){
                                         size_t fsize = get_filesize(block_filename);
                                         if(fsize != inlength[to%P]){
                                             std::cout<<"edata is bad!"<<std::endl;
                                         }
                                         assert(fsize==inlength[to%P]);
                                     }else{
                                         std::cout<<"can not find file: " << block_filename <<std::endl;
                                         assert(false);
                                     }
                                     */
                                    writeedatablock(base_filename, inbuf[to%P], to%P, from%P , P, 1,inlength[to%P]);
                                    __sync_add_and_fetch(&interval_dirt_disk[(T[i].interval + j) % P],-1);
                                }else if(CO ==3){
                                   //check data length
                                   /*
                                   block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 0);
                                   if (file_exists(block_filename)){
                                        size_t fsize = get_filesize(block_filename);
                                        if(fsize != olength[to%P]){
                                           std::cout<<"edata is bad!"<<std::endl;
                                        }
                                        assert(fsize==olength[to%P]);
                                    }else{
                                        std::cout<<"can not find file: " << block_filename <<std::endl;
                                        assert(false);
                                    }
                                    */

                                    writeedatablock(base_filename, (char*)obuf[to % P], to%P, from%P , P, 0,olength[to%P]);
                                     //check data length
                                     /*
                                     block_filename = filename_block_edata<EdgeDataType>(bname, to%P, from%P, P, 1);
                                     if (file_exists(block_filename)){
                                         size_t fsize = get_filesize(block_filename);
                                         if(fsize != inlength[to%P]){
                                             std::cout<<"edata is bad!"<<std::endl;
                                         }
                                         assert(fsize==inlength[to%P]);
                                     }else{
                                         std::cout<<"can not find file: " << block_filename <<std::endl;
                                         assert(false);
                                     }
                                     */
                                    if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                    obuf[to%P] = NULL;
                                    olength[to%P] = -1;

                                    writeedatablock(base_filename, (char*)inbuf[to%P], to%P, from%P , P, 1,inlength[to%P]);
                                    if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                    inbuf[to%P] = NULL;
                                    inlength[to%P] = -1;
                                    __sync_add_and_fetch(&interval_dirt_disk[(T[i].interval + j) % P],-2);
                                    //interval_dirt_disk[(T[i].interval + j) % P] -= 2;
                                }
                            }
                       }else { //send to remote worker...
                                //metrics_entry me7;
                                //std::stringstream ss7;
                                if(j==1){
                                    //ss7<<"sent next time ";
                                    m.start_time("sent next time");
                                }
                                m1 = (T[i].interval +j) % M; //send to other worker.
                                if(CO == 0){
                                        //do nothing.
                                }else if(CO == 1) {
                                    write_to_worker(workerIP[m1], PORT1, (char*)obuf[to % P], from, to, 0, olength[to%P]);
                                }else if(CO == 2) {
                                    write_to_worker(workerIP[m1], PORT1, inbuf[to % P], from, to, 1, inlength[to%P]);
                                }else if(CO == 3){
                                #pragma omp parallel sections
                                    {
                                     #pragma omp section
                                     {
                                        write_to_worker(workerIP[m1], PORT1, (char*)obuf[to % P], from, to, 0, olength[to%P]);
                                        if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                        obuf[to%P] = NULL;
                                        olength[to%P] = -1;
                                      }

                                     #pragma omp section
                                     {
                                        write_to_worker(workerIP[m1], PORT2, inbuf[to % P], from, to, 1, inlength[to%P]);
                                        if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                        inbuf[to%P] = NULL;
                                        inlength[to%P] = -1;
                                     }
                                    }
                                }      

                                if(j==1) m.stop_time("sent next time",false);
                           }
		              }
		         }
		
#else
                /* release obuf */
                if(CO == 2 || CO == 0){
                    for(int n = 0; n < P; n++){
                        if(obuf[n] != NULL) free(obuf[n]);
                        olength[n] = -1;
                        obuf[n] = NULL;
                    }
                }

                /* release inbuf */
                uint8_t *str;
                int len;
                if(CO == 1 || CO == 0){
                    //only commit main shard.
                    memoryshard->dynamicblocks[exec_interval]->write(&str,len);
                    if(main_to_main_buf!=NULL) free(main_to_main_buf);
                    main_to_main_buf = (char*)malloc(len+1);
                    memcpy(main_to_main_buf,str,len);
                    ll_main_to_main = len;
                    free(str);
                    //if(CO == 1) writeedatablock(base_filename,(char*)inbuf[exec_interval], exec_interval, exec_interval , P, 0,len);
                    for(int m = 0; m < P; m++){
                        if(inbuf[m] != NULL) free(inbuf[m]);
                        delete memoryshard->dynamicblocks[m];
                        memoryshard->dynamicblocks[m] = NULL;
                        inbuf[m] = NULL;
                        inlength[m] = -1;
                    }
                }

                if(CO == 3) {
                    memoryshard->dynamicblocks[exec_interval]->write(&str,len);
                    if(main_to_main_buf!=NULL) free(main_to_main_buf);
                    main_to_main_buf = (char*)malloc(len+1);
                    memcpy(main_to_main_buf,str,len);
                    free(str);
                    ll_main_to_main = len;

                    //writeedatablock(base_filename,(char*)str, exec_interval, exec_interval , P, 0,len);
                    if(inbuf[exec_interval] != NULL) free(inbuf[exec_interval]);
                    delete memoryshard->dynamicblocks[exec_interval];
                    memoryshard->dynamicblocks[exec_interval] = NULL;
                    inbuf[exec_interval] = NULL;
                    inlength[exec_interval] = -1;
                }
                int ssize;
		        if(CO != 0){
		            for(int j = 1; j < P; j++){
                       int from, to;
                       from = T[i].interval;
                       to   = T[i].interval + j;
			            if(j == 1){//inform next interval to begin.
                	    	/* send to xgMaster to next */
                		    while((ret = send(fd_master, "XGREQU", 6, 0)) < 6){
                        		std::cout<<"Send Successful signal to xgMaster..."<<std::endl;
                        		if(ret == -1){
                              			std::cout<<"Recv fail, exit..."<<std::endl;
                              			exit(0);
                        		}
                		    }//this interval have been updated
			            }
                        if((to - from) % M == 0){//local
                           if( j == M ){
                                //in-memory, need not saving to disk
                                if(CO == 0){
                                    //do nothing
                                }else if(CO == 1) {
						            if(to > P * N -1){
                                        __sync_add_and_fetch(&cout_g, 1);
                                        uint8_t *str;
                                        int len;
                                        if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                        sliding_shards[to%P]->dyb->write(&str, len);
                                        obuf[to%P] = (char*)malloc(len+1);
                                        memcpy(obuf[to%P],str,len);
                                        olength[to%P] = len;
                                        writeedatablock(base_filename,(char*)obuf[to%P], to%P, from%P , P, 0,olength[to%P]);
                                        free(str);
                                    }else{
                                        __sync_add_and_fetch(&cout_g, 1);
                                        int len;
                                        uint8_t *str;
                                        sliding_shards[to%P]->dyb->write(&str, len);
                                        ll_out_to_in = len;
                                        olength[to%P] = len;
                                        if(out_to_in_buf != NULL) free(out_to_in_buf);
                                        out_to_in_buf = (char*)malloc(len+1);
                                        memcpy(out_to_in_buf, str, len);
                                        if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                        obuf[to%P] = (char*)malloc(len+1);
                                        memcpy(obuf[to%P], str, len);
                                        free(str);
                                    }
                               }else if(CO == 2) {
						            if(to > P * N -1){
                                        uint8_t *str;
                                        int len;
                                        memoryshard->dynamicblocks[to%P]->write(&str,len);
                                        if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                        inbuf[to%P] = (char*)malloc(len+1);
                                        inlength[to%P] = len;
                                        memcpy(inbuf[to%P],str,len);
                                        writeedatablock(base_filename, (char*)str, to%P, from%P , P, 1,len);
                                        free(str);
                                    }else{
                                        uint8_t *str;
                                        int len;
                                        memoryshard->dynamicblocks[to%P]->write(&str,len);
                                        ll_in_to_out = len;
                                        inlength[to%P] = len;
                                        if(in_to_out_buf!=NULL) free(in_to_out_buf);
                                        in_to_out_buf = (char*)malloc(len+1);
                                        memcpy(in_to_out_buf,str,len);
                                        if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                        inbuf[to%P] = (char*)malloc(len+1);
                                        memcpy(inbuf[to%P],str,len);
                                        free(str);
                                    }
                               }else if(CO ==3){
						            if(to > P * N -1){
                                        uint8_t *str;
                                        int len;
                                        if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                        obuf[to%P] = NULL;
                                        sliding_shards[to%P]->dyb->write(&str, len);
                                        olength[to%P] = -1;
                                        writeedatablock(base_filename,(char*)str, to%P, from%P , P, 0,len);
                                        free(str);

                                        memoryshard->dynamicblocks[to%P]->write(&str,len);
                                        if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                        inbuf[to%P] = NULL;
                                        writeedatablock(base_filename, (char*)str, to%P, from%P , P, 1,len);
                                        inlength[to%P] = -1;
                                        free(str);

                                    }else{
                                        int len;
                                        uint8_t *str;
                                        sliding_shards[to%P]->dyb->write(&str, len);
                                        ll_out_to_in = len;
                                        olength[to%P] = -1;
                                        if(out_to_in_buf != NULL) free(out_to_in_buf);
                                        out_to_in_buf = (char*)malloc(len+1);
                                        memcpy(out_to_in_buf, str, len);
                                        if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                        obuf[to%P] = NULL;
                                        free(str);

                                        memoryshard->dynamicblocks[to%P]->write(&str,len);
                                        ll_in_to_out = len;
                                        inlength[to%P] = -1;
                                        if(in_to_out_buf!=NULL) free(in_to_out_buf);
                                        in_to_out_buf = (char*)malloc(len+1);
                                        memcpy(in_to_out_buf,str,len);
                                        if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                        inbuf[to%P] = NULL;
                                        free(str);

                                    }
                                }
                           }else {
                                //save to disk
                                if(CO == 0){
                                     //do nothing.
                                }else if(CO == 1) {
                                    uint8_t *str;
                                    int len;
                                    sliding_shards[to%P]->dyb->write(&str, len);
                                    if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                    obuf[to%P] = (char*)malloc(len+1);
                                    memcpy(obuf[to%P],str, len);
                                    writeedatablock(base_filename,(char*)obuf[to%P], to%P, from%P , P, 0,len);
                                    olength[to%P] = len;
                                    free(str);
                                      __sync_add_and_fetch(&cout_g, 1);
                                     __sync_add_and_fetch(&interval_dirt_disk[(T[i].interval + j) % P],-1);
                                }else if(CO == 2){
                                    uint8_t *str;
                                    int len;
                                    memoryshard->dynamicblocks[to%P]->write(&str,len);
                                    if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                    inbuf[to%P] = (char*)malloc(len+1);
                                    memcpy(inbuf[to%P],str,len);
                                    inlength[to%P] = len;
                                    writeedatablock(base_filename, (char*)str, to%P, from%P , P, 1,len);
                                    free(str);
                                    __sync_add_and_fetch(&interval_dirt_disk[(T[i].interval + j) % P],-1);
                                }else if(CO ==3){
                                    uint8_t *str;
                                    int len;
                                    sliding_shards[to%P]->dyb->write(&str, len);
                                    if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                    obuf[to%P] = NULL;
                                    writeedatablock(base_filename,(char*)str, to%P, from%P , P, 0,len);
                                    olength[to%P] = -1;
                                    free(str);

                                    memoryshard->dynamicblocks[to%P]->write(&str,len);
                                    if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                    inbuf[to%P] = NULL;
                                    inlength[to%P] = -1;
                                    writeedatablock(base_filename, (char*)str, to%P, from%P , P, 1,len);
                                    free(str);

                                    __sync_add_and_fetch(&interval_dirt_disk[(T[i].interval + j) % P],-2);
                                }
                            }
                       }else { //send to remote worker...
                                m1 = (T[i].interval +j) % M; //send to other worker.
                                if(CO == 0){
                                        //do nothing.
                                }else if(CO == 1) {
                                    uint8_t *str;
                                    int len;
                                    sliding_shards[to%P]->dyb->write(&str, len);
                                    if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                    obuf[to%P] = (char*)malloc(len+1);
                                    memcpy(obuf[to%P],str,len);
                                    free(str);
                                    olength[to%P] = len;
                                    std::cout<<"len:"<<len<<std::endl;
                                    write_to_worker(workerIP[m1], PORT1, (char*)obuf[to%P], from, to, 0, len);
                                }else if(CO == 2) {
                                    uint8_t *str;
                                    int len;
                                    memoryshard->dynamicblocks[to%P]->write(&str,len);
                                    if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                    inbuf[to%P] = (char*)malloc(len+1);
                                    memcpy(inbuf[to%P],str,len);
                                    free(str);
                                    inlength[to%P] = len;
                                    write_to_worker(workerIP[m1], PORT1, (char*)inbuf[to%P], from, to, 1, len);
                                }else if(CO == 3){
                                #pragma omp parallel sections
                                    {
                                     #pragma omp section
                                        {
                                            uint8_t *str;
                                            int len;
                                            sliding_shards[to%P]->dyb->write(&str, len);
                                            if(obuf[to%P]!=NULL) free(obuf[to%P]);
                                            write_to_worker(workerIP[m1], PORT1, (char*)str, from, to, 0, len);
                                            obuf[to%P] = NULL;
                                            olength[to%P] = -1;
                                            free(str);
                                        }

                                     #pragma omp section
                                        {
                                            uint8_t *str;
                                            int len;
                                            memoryshard->dynamicblocks[to%P]->write(&str,len);
                                            if(inbuf[to%P]!=NULL) free(inbuf[to%P]);
                                            write_to_worker(workerIP[m1], PORT2, (char*)str, from, to, 1, len);
                                            inbuf[to%P] = NULL;
                                            inlength[to%P] = -1;
                                            free(str);
                                        }
                                    }
                                }      

                           }
		              }
		         }
#endif
                 m.stop_time(me6,ss6.str(),false);

                /* Save vertices */

                if (!disable_vertexdata_storage) {
                    save_vertices(vertices);
                }

				/* Delete edge buffer. TODO: reuse. */
                if((CO == 0 && i == c && P==M)||CO > 0||(CO==0 && M<P)){
                    if (edata != NULL) {
                        free(edata);
                        edata = NULL;
                    }
                }

                /* write related shard blocks to disk. */
                if(CO==1 || CO==3){
                    writeedatablock(base_filename,main_to_main_buf, exec_interval, exec_interval , P, 0, ll_main_to_main);
                    free(main_to_main_buf);
                    main_to_main_buf = NULL;
                }

                if ((memoryshard->loaded() && (save_edgesfiles_after_inmemmode || !is_inmemory_mode()))&&(CO > 0)) {
                    logstream(LOG_INFO) << "Commit memshard" << std::endl;

                    memoryshard->commit_XG(modifies_inedges, modifies_outedges & !disable_outedges); //todo:need it or not??
                    delete memoryshard;
                    memoryshard = NULL;
                }
                if (!is_inmemory_mode())
                    userprogram.after_exec_interval(interval_st, interval_en, chicontext); //do nothing
         
		        if(((T[i].iter < T[i+1].iter)&&(i < c -1)) || (i == c -1)) { //Inter finished.
			        if (!is_inmemory_mode())  // Run sepately
				        userprogram.after_iteration(iter, chicontext);
			        /* Write progress log */
			        write_delta_log();
			        iteration_finished();
		        }

	    	    /* save outedge blocks to disk */
                if(CO == 1){
                int volatile c1 = 0;
                #pragma omp parallel for
                    for(int i = 0; i < P; i++){
                        uint8_t *str;
                        int len;
                        if(i != exec_interval){
                            writeedatablock(base_filename, (char*)obuf[i], exec_interval, i , P, 1,olength[i]); 
#ifdef DYNAMICEDATA
                            delete sliding_shards[i]->dyb;
                            sliding_shards[i]->dyb = NULL;
#endif
                        } 
                        if(obuf[i] != NULL) free(obuf[i]);
                        obuf[i] = NULL;
                        olength[i] = -1;
                        __sync_add_and_fetch(&c1, 1);
                    }
                    while(c1 < P);
                }
                /* End */
		        interval_dirt_main[exec_interval] = 0;	

                m.stop_time(me2,ss2.str(),false);

            } // Intervals

            /* release slidingshards */
           for(int i = 0; i < P; i++){
               //if(sliding_shards[i] != NULL) delete sliding_shards[i];
               //sliding_shards[i] = NULL;
           }
            m.stop_time("runtime");

            m.set("updates", nupdates);
            m.set("work", work);
            m.set("nvertices", num_vertices());
            m.set("execthreads", (size_t)exec_threads);
            m.set("loadthreads", (size_t)load_threads);
#ifndef GRAPHCHI_DISABLE_COMPRESSION
            m.set("compression", 1);
#else
            m.set("compression", 0);
#endif

            m.set("scheduler", (size_t)use_selective_scheduling);
            m.set("niters", niters);

            // Close outputs
            for(int i=0; i< (int)outputs.size(); i++) {
                outputs[i]->close();
            }
            outputs.clear();
            /* Send to ntoplist to xgMaster */
#ifndef VERTEX_TYPE_INT
            char *topbuf;
            int i=0;
            topbuf = (char*)malloc(sizeof(int)*ntop + sizeof(float)*ntop +1+sizeof(int));
            memcpy(topbuf,&i,sizeof(int));
            std::vector< vertex_value<float> > top = get_top_vertices<float>(base_filename, ntop);
            for(int i=0; i < (int)top.size(); i++) {
                 memcpy(topbuf + sizeof(int) + i*(sizeof(int)+sizeof(float)),&top[i].vertex,sizeof(int));
                 memcpy(topbuf + sizeof(int) + i*(sizeof(int)+sizeof(float))+sizeof(int),&top[i].value,sizeof(float));
            }
            send(fd_master, topbuf, ntop*(sizeof(int)+sizeof(float))+sizeof(int), 0); 
#else
            /*
        std::vector< vertex_value<VertexDataType> > top1 = get_top_vertices<VertexDataType>(base_filename, 20);
            std::cout << "Print top 20 vertices: " << std::endl;
                for(int i=0; i < (int) top1.size(); i++) {
                        std::cout << (i+1) << ". " << top1[i].vertex << "\t" << top1[i].value << std::endl;
                    
                }
                */
            char *topbuf;
            int i=1;
            topbuf = (char*)malloc(sizeof(int)*ntop + sizeof(int)*ntop +1+sizeof(int));
            memcpy(topbuf,&i,sizeof(int));
            std::vector< vertex_value<int> > top = get_top_vertices<int>(base_filename, ntop);
            for(int i=0; i < (int)top.size(); i++) {
                 memcpy(topbuf + sizeof(int) + i*(sizeof(int)+sizeof(int)),&top[i].vertex,sizeof(int));
                 memcpy(topbuf + sizeof(int) + i*(sizeof(int)+sizeof(int))+sizeof(int),&top[i].value,sizeof(int));
                 //std::cout<<" v: "<<top[i].vertex<<"value: "<<top[i].value<<std::endl;
            }
            send(fd_master, topbuf, ntop*(sizeof(int)+sizeof(int))+sizeof(int), 0); 
#endif
	        recv(fd_master, buf, 1024,0);
	        close(fd_master);
	        close(sockfd1);
	        close(sockfd2);

            std::cout<<"Finished..."<<cout_g<<std::endl;
            std::cout<<"prev vertexs : "<<cc_prev<<std::endl;
        }
        
        virtual void iteration_finished() {
            // Do nothing
        }
        
        stripedio * get_iomanager() {
            return iomgr;
        }
        
        virtual void set_modifies_inedges(bool b) {
            modifies_inedges = b;
        }
        
        virtual void set_modifies_outedges(bool b) {
            modifies_outedges = b;
        }
        
        virtual void set_only_adjacency(bool b) {
            only_adjacency = b;
        }

        virtual void set_preload_commit(bool b){
            preload_commit = b;
        }
        
        virtual void set_disable_outedges(bool b) {
            disable_outedges = b;
        }
        
        /**
         * Configure the blocksize used when loading shards.
         * Default is one megabyte.
         * @param blocksize_in_bytes the blocksize in bytes
         */
        void set_blocksize(size_t blocksize_in_bytes) {
            blocksize = blocksize_in_bytes;
        }
        
        /**
         * Set the amount of memory available for loading graph
         * data. Default is 1000 megabytes.
         * @param mbs amount of memory to be used.
         */
        void set_membudget_mb(int mbs) {
            membudget_mb = mbs;
        }
        
        
        void set_load_threads(int lt) {
            load_threads = lt;
        }
        
        void set_exec_threads(int et) {
            exec_threads = et;
        }
        
        /**
         * Sets whether the engine is run in the deterministic
         * mode. Default true.
         */
        void set_enable_deterministic_parallelism(bool b) {
#ifdef DYNAMICEDATA
            if (!b) {
                logstream(LOG_ERROR) << "With dynamic edge data, you cannot disable determinic parallelism." << std::endl;
                logstream(LOG_ERROR) << "Otherwise race conditions would corrupt the structure of the data." << std::endl;
                assert(b);
                return;
            }
#endif
            enable_deterministic_parallelism = b;
        }
      
    public:
        void set_disable_vertexdata_storage() {
            this->disable_vertexdata_storage = true;
        }
        
        void set_enable_vertexdata_storage() {
            this->disable_vertexdata_storage = false;
        }
       
        void set_maxwindow(unsigned int _maxwindow){ 
            maxwindow = _maxwindow;
        }; 
        
        
        /* Outputs */
        size_t add_output(ioutput<VertexDataType, EdgeDataType> * output) {
            outputs.push_back(output);
            return (outputs.size() - 1);
        }
         
        ioutput<VertexDataType, EdgeDataType> * output(size_t idx) {
            if (idx >= outputs.size()) {
                logstream(LOG_FATAL) << "Tried to get output with index " << idx << ", but only " << outputs.size() << " outputs were initialized!" << std::endl;
            }
            assert(idx < outputs.size());
            return outputs[idx];
        }
        
    protected:
              
        virtual void _load_vertex_intervals() {
            load_vertex_intervals(base_filename, nshards, intervals);
        }
        
    protected:
        mutex httplock;
        std::map<std::string, std::string> json_params;
        
    public:
        
        /**
         * Replace all shards with zero values in edges.
         */
        template<typename ET>
        void reinitialize_edge_data(ET zerovalue) {
            
            for(int p=0; p < nshards; p++) {
                std::string edatashardname =  filename_shard_edata<ET>(base_filename, p, nshards);
                std::string dirname = dirname_shard_edata_block(edatashardname, blocksize);
                size_t edatasize = get_shard_edata_filesize<ET>(edatashardname);
                logstream(LOG_INFO) << "Clearing data: " << edatashardname << " bytes: " << edatasize << std::endl;
                int nblocks = (int) ((edatasize / blocksize) + (edatasize % blocksize == 0 ? 0 : 1));
                for(int i=0; i < nblocks; i++) {
                    std::string block_filename = filename_shard_edata_block(edatashardname, i, blocksize);
                    int len = (int) std::min(edatasize - i * blocksize, blocksize);
                    int f = open(block_filename.c_str(), O_RDWR | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
                    ET * buf =  (ET *) malloc(len);
                    for(int i=0; i < (int) (len / sizeof(ET)); i++) {
                        buf[i] = zerovalue;
                    }
                    write_compressed(f, buf, len);
                    close(f);
                    
#ifdef DYNAMICEDATA
                    write_block_uncompressed_size(block_filename, len);
#endif
                    
                }
            }
        }
        
        
        /**
          * If true, the vertex data is initialized before
          * the engineis started. Default false.
          */
        void set_reset_vertexdata(bool reset) {
            reset_vertexdata = reset;
        }
        
        
        /**
         * Whether edges should be saved after in-memory mode
         */
        virtual void set_save_edgesfiles_after_inmemmode(bool b) {
            this->save_edgesfiles_after_inmemmode = b;
        }

        virtual void set_initialize_edges_before_run(bool b) {
            this->initialize_edges_before_run = b;
        }
        
        
        /**
         * HTTP admin management
         */
        
        void set_json(std::string key, std::string value) {
            httplock.lock();
            json_params[key] = value;
            httplock.unlock();
        }
        
        template <typename T>
        void set_json(std::string key, T val) {
            std::stringstream ss;
            ss << val;
            set_json(key, ss.str());
        }
        
        std::string get_info_json() {
            std::stringstream json;
            json << "{";
            json << "\"file\" : \"" << base_filename << "\",\n";
            json << "\"numOfShards\": " << nshards << ",\n";
            json << "\"iteration\": " << chicontext.iteration << ",\n";
            json << "\"numIterations\": " << chicontext.num_iterations << ",\n";
            json << "\"runTime\": " << chicontext.runtime() << ",\n";
            
            json << "\"updates\": " << nupdates << ",\n";
            json << "\"nvertices\": " << chicontext.nvertices << ",\n";
            json << "\"interval\":" << exec_interval << ",\n";
            json << "\"windowStart\":" << sub_interval_st << ",";
            json << "\"windowEnd\": " << sub_interval_en << ",";
            json << "\"shards\": [";
            
            for(int p=0; p < (int)nshards; p++) {
                if (p>0) json << ",";
                
                json << "{";
                json << "\"p\": " << p << ", ";
                json << sliding_shards[p]->get_info_json();
                json << "}";
            }
            
            json << "]";
            json << "}";
            return json.str();
        }
        
    };
    
    
};



#endif


