
/**
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
 * The memory shard. This class should only be accessed internally by the GraphChi engine.
 */

#ifdef DYNAMICEDATA
#include "shards/dynamicdata/memoryshard.hpp"
#else

#ifndef DEF_GRAPHCHI_MEMSHARD
#define DEF_GRAPHCHI_MEMSHARD


#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <string>
#include <pthread.h>

#include "api/graph_objects.hpp"
#include "metrics/metrics.hpp"
#include "io/stripedio.hpp"
#include "ddgraph_types.hpp"

//static int exec_interval;
namespace ddgraph {
    
    
    template <typename VT, typename ET, typename svertex_t = ddgraph_vertex<VT, ET> >
    class memory_shard {
        
        stripedio * iomgr;
        
        std::string filename_edata;
        std::string filename_adj;
        
        vid_t range_st;
        vid_t range_end;
        size_t adjfilesize;
        size_t edatafilesize;
        
        size_t edgeptr;
        vid_t streaming_offset_vid;
        size_t streaming_offset; // The offset where streaming should continue
        size_t range_start_offset; // First byte for this range's vertices (used for writing only outedges)
        size_t range_start_edge_ptr;
        size_t streaming_offset_edge_ptr;
        uint8_t * adjdata;
        int * doneptr;
        std::vector<size_t> blocksizes;
        uint64_t chunkid;
        
        std::vector<int> block_edatasessions;
        int adj_session;
        
        bool async_edata_loading;
        bool is_loaded;
        bool disable_async_writes;
        size_t blocksize;
        metrics &m;
        
    public:
        bool only_adjacency;
        char ** edgedata;
        
        memory_shard(stripedio * iomgr,
                     std::string _filename_edata,
                     std::string _filename_adj,
                     vid_t _range_start,
                     vid_t _range_end,
                     size_t _blocksize,
                     metrics &_m) : iomgr(iomgr), filename_edata(_filename_edata),
        filename_adj(_filename_adj),
        range_st(_range_start), range_end(_range_end), blocksize(_blocksize),  m(_m) {
            adjdata = NULL;
            only_adjacency = false;
            is_loaded = false;
            adj_session = -1;
            edgedata = NULL;
            doneptr = NULL;
            disable_async_writes = false;
            async_edata_loading = !svertex_t().computational_edges();
#ifdef SUPPORT_DELETIONS
            async_edata_loading = false; // See comment above for memshard, async_edata_loading = false;
#endif
        }
        
        ~memory_shard() {
            int nblocks = (int) block_edatasessions.size();
            
            for(int i=0; i < nblocks; i++) {
                if (inbuf[i] != NULL) {
                    iomgr->managed_release(block_edatasessions[i], &inbuf[i]);
                    iomgr->close_session(block_edatasessions[i]);
                }
            }
            if (adj_session >= 0) {
                if (adjdata != NULL) iomgr->managed_release(adj_session, &adjdata);
                iomgr->close_session(adj_session);
            }
            /*
            if (inbuf != NULL)
                free(inbuf);
            inbuf = NULL;
            */
            if (doneptr != NULL) {
                free(doneptr);
            }
        }
       
        int check_segment(size_t index){
            for(int i = 0; i < P; i++){
                    if(index < blocksizes[i]) return(i);
                     index -= blocksizes[i];
            }
            return(-1);
       }

        void set_disable_async_writes(bool b) {
            disable_async_writes = b;
        }
        
        void commit(bool commit_inedges, bool commit_outedges) {
            if (block_edatasessions.size() == 0 || only_adjacency) return;
            assert(is_loaded);
            metrics_entry cm = m.start_time();
            
            /**
             * This is an optimization that is relevant only if memory shard
             * has been used in a case where only out-edges are considered.
             * Out-edges are in a continuous "window", while in-edges are
             * scattered all over the shard
             */
            int nblocks = (int) block_edatasessions.size();
            
            if (commit_inedges) {
                int start_stream_block = (int) (range_start_edge_ptr / blocksize);
                
        #pragma omp parallel for
                for(int i=0; i < nblocks; i++) {
                    /* Write asynchronously blocks that will not be needed by the sliding windows on
                     this iteration. */
                    if (i >= start_stream_block || disable_async_writes) {
                        iomgr->managed_pwritea_now(block_edatasessions[i], &edgedata[i], blocksizes[i], 0);
                        iomgr->managed_release(block_edatasessions[i], &edgedata[i]);
                        iomgr->close_session(block_edatasessions[i]);
                        
                        edgedata[i] = NULL;
                        
                    } else {
                        iomgr->managed_pwritea_async(block_edatasessions[i], &edgedata[i], blocksizes[i], 0, true, true);
                        edgedata[i] = NULL;
                    }
                }
            } else if (commit_outedges) {
                size_t last = streaming_offset_edge_ptr;
                if (last == 0){
                    // rollback
                    last = edatafilesize;
                }
                //char * bufp = ((char*)edgedata + range_start_edge_ptr);
                int startblock = (int) (range_start_edge_ptr / blocksize);
                int endblock = (int) (last / blocksize);
#pragma omp parallel for
                for(int i=0; i < nblocks; i++) {
                    if (i >= startblock && i <= endblock) {
                        iomgr->managed_pwritea_now(block_edatasessions[i], &edgedata[i], blocksizes[i], 0);
                    }
                    iomgr->managed_release(block_edatasessions[i], &edgedata[i]);
                    edgedata[i] = NULL;
                    iomgr->close_session(block_edatasessions[i]);
                }
            } else {
                for(int i=0; i < nblocks; i++) {
                    iomgr->close_session(block_edatasessions[i]);
                }
            }
            
            m.stop_time(cm, "memshard_commit");
            
            iomgr->managed_release(adj_session, &adjdata);
            // FIXME: this is duplicated code from destructor
            for(int i=0; i < nblocks; i++) {
                if (edgedata[i] != NULL) {
                    iomgr->managed_release(block_edatasessions[i], &edgedata[i]);
                }
            }
            block_edatasessions.clear();
            is_loaded = false;
        }
        
        void commit_XG(bool commit_inedges, bool commit_outedges) {
            metrics_entry cm = m.start_time();
            int volatile cc = 0;
            if(CO == 2){
            #pragma omp parallel for
                for(int i=0; i < P; i++) {
                    assert(inbuf[i]!=NULL);
                    iomgr->managed_pwritea_now(block_edatasessions[i], &inbuf[i], blocksizes[i], 0);
                    iomgr->managed_release(block_edatasessions[i], &inbuf[i]);
                    iomgr->close_session(block_edatasessions[i]);
                    inbuf[i] = NULL;
                    inlength[i] = -1;
                    __sync_add_and_fetch(&cc, 1);
                }
                while(cc < P);
            }
            
            for(int j = 0; j < P; j ++) iomgr->close_session(block_edatasessions[j]); //close sessions.
            assert(adjdata!=NULL);
            iomgr->managed_release(adj_session, &adjdata);
            iomgr->close_session(adj_session);
            // FIXME: this is duplicated code from destructor
            
            block_edatasessions.clear();
            is_loaded = false;
            m.stop_time(cm, "memshard_commit");
        }
        
        bool loaded() {
            return is_loaded;
        }
        
    private:
        void load_edata() {
            assert(blocksize % sizeof(ET) == 0);
            int nblocks = (int) (edatafilesize / blocksize + (edatafilesize % blocksize != 0));
            edgedata = (char **) calloc(nblocks, sizeof(char*));
            size_t compressedsize = 0;
            int blockid = 0;
            
            if (!async_edata_loading) {
                doneptr = (int *) malloc(nblocks * sizeof(int));
                for(int i=0; i < nblocks; i++) doneptr[i] = 1;
            }
            
            while(blockid < nblocks) {
                std::string block_filename = filename_shard_edata_block(filename_edata, blockid, blocksize);
                if (file_exists(block_filename)) {
                    size_t fsize = std::min(edatafilesize - blocksize * blockid, blocksize);
                    
                    compressedsize += get_filesize(block_filename);
                    int blocksession = iomgr->open_session(block_filename, false, true); // compressed
                    block_edatasessions.push_back(blocksession);
                    blocksizes.push_back(fsize);
                    
                    edgedata[blockid] = NULL;
                    iomgr->managed_malloc(blocksession, &edgedata[blockid], fsize, 0);
                    if (async_edata_loading) {
                        iomgr->managed_preada_async(blocksession, &edgedata[blockid], fsize, 0);
                    } else {
                        iomgr->managed_preada_async(blocksession, &edgedata[blockid], fsize, 0, (volatile int *)&doneptr[blockid]);
                    }
                    blockid++;
                    
                } else {
                    if (blockid == 0) {
                        logstream(LOG_ERROR) << "Shard block file did not exists:" << block_filename << std::endl;
                    }
                    if (blockid < nblocks) {
                        logstream(LOG_ERROR) << "Did not find block " << block_filename << std::endl;
                        logstream(LOG_ERROR) << "Going to exit..." << std::endl;
                    }
                    break;
                }
            }
            
            logstream(LOG_DEBUG) << "Compressed/full size: " << compressedsize * 1.0 / edatafilesize <<
            " number of blocks: " << nblocks << std::endl;
            assert(blockid == nblocks);
            
        }
        
        /* add by YouLi Cheng. 2014/05/08 */
        void load_edata_XG() {
             int nblocks;
             int dirt;
             if(CO == 0){
                 nblocks = P;
             }else if(CO == 1 || CO == 3) { 
                  if(interval<M) dirt = interval;
                  else dirt = mini(M, P-1);
                  nblocks = P - dirt;
             }else if(CO == 2){
                   nblocks = P;
             }    
                    
             size_t compressedsize = 0;
             int blockid = 0;
                   
             if (!async_edata_loading) {
                doneptr = (int *) malloc(nblocks * sizeof(int));//int * doneptr;
                for(int i=0; i < nblocks; i++) doneptr[i] = 1;
             }
                    
             int pos;
             for(int i = 0; i < P; i++){
                 std::string block_filename = filename_block_edata<ET>(bname, exec_interval, i, P, 0); //todo:filename_edata = basename
                 if (file_exists(block_filename)) {
                    size_t fsize = get_filesize(block_filename);                    
                    bool isload=false;
                    pos = exec_interval -1;
                    for(int y = 0; y < nblocks; y++){
                        pos++;
                        if(pos>=P) pos = 0;
                        if(pos == i) {
                            isload=true;
                            break;
                        }
                    }
                    assert(fsize%sizeof(ET)==0);
                    compressedsize += fsize;
                    int blocksession = iomgr->open_session(block_filename, false, true);                       
                    block_edatasessions.push_back(blocksession);
                    blocksizes.push_back(fsize);                    
                    if(isload){
                        iomgr->managed_malloc(blocksession, &(inbuf[i]), fsize+1, 0);
                        //std::cout<<"main:: "<<i<<"  nblocks : "<<nblocks<<"interval_dirt_memory[exec_interval];"<<interval_dirt_memory[exec_interval]<<" M "<<M<<std::endl;
                        inlength[i] = fsize;
                    }
                    else{
                        pthread_mutex_lock(&tlock);
                        if(inlength[i]==-1){
                             //std::cout<<"local i: "<<i<<" inlength[i]: "<< inlength[i]<<" fsize: " << fsize<<std::endl;
                             inlength[i] = fsize;
                             inbuf[i] = NULL;
                             iomgr->managed_malloc(blocksession, &(inbuf[i]), fsize+1, 0);
                             assert(inbuf != NULL);
                        }else{
                            //std::cout<<"remote i: "<<i<<" inlength[i]: "<< inlength[i]<<" fsize: " << fsize<<std::endl;
                            assert(inlength[i]==fsize);
                        }
                        pthread_mutex_unlock(&tlock);
                    }
                    assert(inbuf[i]!=NULL);
                    /* Init edgedata */
                    memset(inbuf[i],'\0',fsize+1);
                    for(int j=0; j < (int) (fsize/sizeof(ET)); j++) *((ET*)((inbuf[i])+sizeof(ET)*j)) = ET();
                  } else {
                        if (blockid == 0) {
                           logstream(LOG_ERROR) << "Shard block file did not exists:" << block_filename << std::endl;
                        }
                        if (blockid < nblocks) {
                           logstream(LOG_ERROR) << "Did not find block " << block_filename << std::endl;
                           logstream(LOG_ERROR) << strerror(errno)<<std::endl;
                           logstream(LOG_ERROR) << "Going to exit..." <<std::endl;
                        }
                         assert(false);
                  }
            } 
             inbuf_loaded = 1;
             for(int j = exec_interval; j < exec_interval + nblocks; j++){
                    if (async_edata_loading) {
                        iomgr->managed_preada_async(block_edatasessions[j%P], &inbuf[j%P], blocksizes[j%P], 0);
                    } else {
                            iomgr->managed_preada_async(block_edatasessions[j%P], &inbuf[j%P], blocksizes[j%P], 0, (volatile int *)&doneptr[j-exec_interval]);
                    }
                    
             }     
    } 
    public:
        
        // TODO: recycle ptr!
        void load() {
            is_loaded = true;
            adjfilesize = get_filesize(filename_adj);
            
#ifdef SUPPORT_DELETIONS
            async_edata_loading = false;  // Currently we encode the deleted status of an edge into the edge value (should be changed!),
            // so we need the edge data while loading
#endif
            
            //preada(adjf, adjdata, adjfilesize, 0);
            
            adj_session = iomgr->open_session(filename_adj, true);
            iomgr->managed_malloc(adj_session, &adjdata, adjfilesize, 0);
            
            /* Load in parallel: replaces older stream solution */
            size_t bufsize = 16 * 1024 * 1024;
            int n = (int) (adjfilesize / bufsize + 1);

#pragma omp parallel for
            for(int i=0; i < n; i++) {
                size_t toread = std::min(adjfilesize - i * bufsize, (size_t)bufsize);
                iomgr->preada_now(adj_session, adjdata + i * bufsize, toread, i * bufsize, true);
            }

            
            /* Initialize edge data asynchonous reading */
            if (!only_adjacency) {
                std::string fn;
                fn = filename_block_edata<ET>(bname, exec_interval, exec_interval, P, 3);
                std::cout<<"filename memory.........."<<std::endl; 
                edatafilesize = get_shard_edata_filesize<ET>(fn);
                load_edata();
            }
        }
        
        void load_XG() {
            is_loaded = true;
            filename_adj = filename_block_adj(bname, exec_interval, exec_interval, P);
            adjfilesize = get_filesize(filename_adj);
            
#ifdef SUPPORT_DELETIONS
            async_edata_loading = false;  // Currently we encode the deleted status of an edge into the edge value (should be changed!),
            // so we need the edge data while loading
#endif
            
            
            adj_session = iomgr->open_session(filename_adj, true);
            adjdata = NULL;
            iomgr->managed_malloc(adj_session, &adjdata, adjfilesize+1, 0);
            assert(adjdata != NULL);
            memset(adjdata, '\0', adjfilesize+1);
            
            /* Load in parallel: replaces older stream solution */
            size_t bufsize = 16 * 1024 * 1024;
            int n = (int) (adjfilesize / bufsize + 1);

#pragma omp parallel for
            for(int i=0; i < n; i++) {
                size_t toread = std::min(adjfilesize - i * bufsize, (size_t)bufsize);
                iomgr->preada_now(adj_session, adjdata + i * bufsize, toread, i * bufsize, true);
            }

            
            /* Initialize edge data asynchonous reading */
            if (!only_adjacency) {
                load_edata_XG();
            }
        }
        
    
        void load_vertices(vid_t window_st, vid_t window_en, std::vector<svertex_t> & prealloc, bool inedges=true, bool outedges=true) {
            //do nothing...
        }
        void load_vertices_XG(vid_t window_st, vid_t window_en, std::vector<svertex_t> & prealloc, bool inedges=true, bool outedges=true) {
            /* Find file size */
            m.start_time("memoryshard_create_edges");
            
            assert(adjdata != NULL);
            
            // Now start creating vertices
            uint8_t * ptr = adjdata;
            uint8_t * end = ptr + adjfilesize;
            vid_t vid = 0;
            edgeptr = 0;
            
            streaming_offset = 0;
            streaming_offset_vid = 0;
            streaming_offset_edge_ptr = 0;
            range_start_offset = adjfilesize;
            range_start_edge_ptr = edatafilesize;
            
            bool setoffset = false;
            bool setrangeoffset = false;
            size_t bid = 0;
            size_t brange = blocksizes[0];
            size_t edatasize = 0;
            for(int i = 0; i < P; i ++){
                edatasize += blocksizes[i];
            }
            while (ptr < end) {
                if (!setoffset && vid > range_end) {
                    // This is where streaming should continue. Notice that because of the
                    // non-zero counters, this might be a bit off.
                    streaming_offset = ptr-adjdata;
                    streaming_offset_vid = vid;
                    streaming_offset_edge_ptr = edgeptr;
                    setoffset = true;
                }
                if (!setrangeoffset && vid>=range_st) {
                    range_start_offset = ptr-adjdata;
                    range_start_edge_ptr = edgeptr;
                    setrangeoffset = true;
                }
                
                uint8_t ns = *ptr;
                int n;
                
                ptr += sizeof(uint8_t);
                
                if (ns == 0x00) {
                    // next value tells the number of vertices with zeros
                    uint8_t nz = *ptr;
                    ptr += sizeof(uint8_t);
                    vid++;
                    vid += nz;
                    continue;
                }
                
                if (ns == 0xff) {  // If 255 is not enough, then stores a 32-bit integer after.
                    n = *((uint32_t*)ptr);
                    ptr += sizeof(uint32_t);
                } else {
                    n = ns;
                }
                svertex_t* vertex = NULL;
                
                if (vid>=window_st && vid <=window_en) { // TODO: Make more efficient
                    vertex = &prealloc[vid-window_st];
                    if (!vertex->scheduled) vertex = NULL;
                }
                bool any_edges = false;
                while(--n>=0) {
                    int blockid = (int) (edgeptr / blocksize);

                    if (!async_edata_loading && !only_adjacency) {
                        /* Wait until blocks loaded (non-asynchronous version) */
                        while(doneptr[edgeptr / blocksize] != 0) { usleep(10); }
                    }
                    
                    vid_t target = *((vid_t*) ptr);
                    ptr += sizeof(vid_t);
                    if (vertex != NULL && outedges)
                    {
                        char * eptr = (only_adjacency ? NULL  : inbuf[bid]+(edgeptr - (brange-blocksizes[bid])));
                        vertex->add_outedge(target, (only_adjacency ? NULL : (ET*) eptr), false);
                    }
                    
                    if (target >= window_st)  {
                        if (target <= window_en) {                        /* In edge */
                            if (inedges) {
                                svertex_t & dstvertex = prealloc[target - window_st];
                                /* for previos updaete */
                                if((vid<window_st || vid>window_en)&&(CO==1 || CO==3)) dstvertex.crucial_flag = true;
                                //prev_bitset->clear_bit(target-window_st);

                                if (dstvertex.scheduled) {
                                    any_edges = true;
                                    //  assert(only_adjacency ||  edgeptr < edatafilesize);
                                    char * eptr = (only_adjacency ? NULL  : inbuf[bid]+(edgeptr - (brange-blocksizes[bid])));
                                    
                                    dstvertex.add_inedge(vid,  (only_adjacency ? NULL : (ET*) eptr), false);
                                    dstvertex.parallel_safe = dstvertex.parallel_safe && (vertex == NULL); // Avoid if
                                }
                            }
                        } else { // Note, we cannot skip if there can be "special edges". FIXME so dirty.
                            // This vertex has no edges any more for this window, bail out
                            if (vertex == NULL) {
                                ptr += sizeof(vid_t) * n;
                                edgeptr += (n + 1) * sizeof(ET);
                                /* Add by Cyl */
                                while(edgeptr>=brange) {
                                    bid++;
                                    brange += blocksizes[bid];
                                }
                                break;
                            }
                        }
                    }

                    edgeptr += sizeof(ET);
                    if(edgeptr>=brange) {
                        bid++;
                        brange += blocksizes[bid];
                    } 
                }
                
                if (any_edges && vertex != NULL) {
                    vertex->parallel_safe = false;
                }
                vid++;
            }
            m.stop_time("memoryshard_create_edges", false);
        }
        
        size_t offset_for_stream_cont() {
            return streaming_offset;
        }
        vid_t offset_vid_for_stream_cont() {
            return streaming_offset_vid;
        }
        size_t edata_ptr_for_stream_cont() {
            return streaming_offset_edge_ptr;
        }
        
        
        
        
    };
};

#endif
#endif
