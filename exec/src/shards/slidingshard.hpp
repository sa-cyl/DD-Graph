

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
 * The sliding shard.
 */

#ifdef DYNAMICEDATA
#include "shards/dynamicdata/slidingshard.hpp"
#else

#ifndef DEF_GRAPHCHI_SLIDINGSHARD
#define DEF_GRAPHCHI_SLIDINGSHARD


#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <string>

#include "api/graph_objects.hpp"
#include "metrics/metrics.hpp"
#include "logger/logger.hpp"
#include "io/stripedio.hpp"
#include "ddgraph_types.hpp"


namespace ddgraph {
    
    
    
    
    /**
     * A streaming block.
     */
    struct sblock {
        
        int writedesc;
        int readdesc;
        size_t offset;
        size_t end;
        uint8_t * data;
        uint8_t * ptr;
        bool active;
        bool is_edata_block;
        
        sblock() : writedesc(0), readdesc(0), active(false) { data = NULL; }
        sblock(int wdesc, int rdesc, bool is_edata_block=false) : writedesc(wdesc), readdesc(rdesc), active(false),
        is_edata_block(is_edata_block){ data = NULL; }
        
        void commit_async(stripedio * iomgr) {
            if (active && data != NULL && writedesc >= 0) {
                if (is_edata_block) {
                    iomgr->managed_pwritea_async(writedesc, &data, end-offset, 0, true, true);
                    data = NULL;
                } else {
                    iomgr->managed_pwritea_async(writedesc, &data, end-offset, offset, true);
                }
            }
        }
        
        void commit_now(stripedio * iomgr) {
            if (active && data != NULL && writedesc >= 0) {
                size_t len = ptr-data;
                if (len > end-offset) len = end-offset;
                if (is_edata_block) {
                    iomgr->managed_pwritea_now(writedesc, &data, end - offset, 0); /* Need to write whole block in the compressed regime */
                } else {
                    iomgr->managed_pwritea_now(writedesc, &data, len, offset);
                }
            }
        }
        void read_async(stripedio * iomgr) {
            if (is_edata_block) {
                iomgr->managed_preada_async(readdesc, &data, (end - offset), 0);
                
            } else {
                iomgr->managed_preada_async(readdesc, &data, end - offset, offset);
            }
        }
        void read_now(stripedio * iomgr) {
            if (is_edata_block) {
                iomgr->managed_preada_now(readdesc, &data, end-offset, 0);
            } else {
                iomgr->managed_preada_now(readdesc, &data, end-offset, offset);
            }
        }
        
        void release(stripedio * iomgr) {
            if (data != NULL) {
                iomgr->managed_release(readdesc, &data);
                if (is_edata_block) {
                    iomgr->close_session(readdesc);
                }
            }
            data = NULL;
            
        }
    };
    
    
    struct indexentry {
        size_t adjoffset, edataoffset;
        indexentry(size_t a, size_t e) : adjoffset(a), edataoffset(e) {}
    };
    
    /*
     * Graph shard that is streamed. I.e, it can only read in one direction, a chunk
     * a time.
     */
    template <typename VT, typename ET, typename svertex_t = ddgraph_vertex<VT, ET>, typename ETspecial = ET>
    class sliding_shard {
        
        stripedio * iomgr;
        
        std::string filename_edata;
        std::string filename_adj;
        vid_t range_st, range_end;
        size_t blocksize;
        
        vid_t curvid;
        size_t adjoffset, edataoffset, adjfilesize, edatafilesize;
        size_t window_start_edataoffset;
        
        std::vector<sblock> activeblocks;
        int adjfile_session;
        int writedesc;
        sblock * curblock;
        sblock * curadjblock;
        metrics &m;
        
        std::map<int, indexentry> sparse_index; // Sparse index that can be created in the fly
        bool disable_writes;
        bool async_edata_loading;
        bool disable_async_writes;
        // bool need_read_outedges; // Disabled - does not work with compressed data: whole block needs to be read.
        
        
    public:
        bool only_adjacency;
        size_t bsize;
        
        sliding_shard(stripedio * iomgr, std::string _filename_edata, std::string _filename_adj, vid_t _range_st, vid_t _range_en, size_t _blocksize, metrics &_m,
                      bool _disable_writes=false, bool onlyadj = false) :
        iomgr(iomgr),
        filename_edata(_filename_edata),
        filename_adj(_filename_adj),
        range_st(_range_st),
        range_end(_range_en),
        blocksize(_blocksize),
        m(_m),
        disable_writes(_disable_writes) {
            curvid = 0;
            adjoffset = 0;
            edataoffset = 0;
            disable_writes = false;
            only_adjacency = onlyadj;
            curblock = NULL;
            curadjblock = NULL;
            window_start_edataoffset = 0;
            disable_async_writes = false;
            
            while(blocksize % sizeof(ET) != 0) blocksize++;
            assert(blocksize % sizeof(ET)==0);
            
            //adjfilesize = get_filesize(filename_adj);
            if (!only_adjacency) {
                 std::string fn;
                 fn = filename_block_edata<ET>(bname, exec_interval, exec_interval, P, 3);  
                 //logstream(LOG_DEBUG) << "Total edge data size: " << edatafilesize  << ", " << filename_edata << "sizeof(ET): " << sizeof(ET) << std::endl;
            } else {
                // Nothing
            }
            
            //adjfile_session = iomgr->open_session(filename_adj, true);
            save_offset();
            
            async_edata_loading = !svertex_t().computational_edges();
#ifdef SUPPORT_DELETIONS
            async_edata_loading = false; // See comment above for memshard, async_edata_loading = false;
#endif
        }
        
        ~sliding_shard() {
            release_prior_to_offset(true);
            if (curblock != NULL) {
                curblock->release(iomgr);
                delete curblock;
                curblock = NULL;
            }
            if (curadjblock != NULL) {
                curadjblock->release(iomgr);
                delete curadjblock;
                curadjblock = NULL;
            }            
            iomgr->close_session(adjfile_session);
        }
        
        
        size_t num_edges() {
            return edatafilesize / sizeof(ET);
        }
        
        // Init edge data blocks
        void initdata() {
            logstream(LOG_DEBUG) << "Initialize edge data: " << filename_edata << std::endl;
            ET * initblock = (ET *) malloc(blocksize);
            for(int i=0; i < (int) (blocksize/sizeof(ET)); i++) initblock[i] = ET();
            for(size_t off=0; off < edatafilesize; off+=blocksize) {
                std::string blockfilename = filename_shard_edata_block(filename_edata, (int) (off / blocksize), blocksize);
                size_t len = std::min(blocksize, edatafilesize - off);
                int f = open(blockfilename.c_str(), O_WRONLY);
                pwritea(f, initblock, len, 0);
                close(f);
            }
            free(initblock);
        }
        
    protected:
        size_t get_adjoffset() { return adjoffset; }
        size_t get_edataoffset() { return edataoffset; }
        
        void save_offset() {
            // Note, so that we can use the lower bound operation in map, we need
            // to insert indices in reverse order
            sparse_index.insert(std::pair<int, indexentry>(-((int)curvid), indexentry(adjoffset, edataoffset)));
        }
        
        void move_close_to(vid_t v) {
            if (curvid >= v) return;
            
            std::map<int,indexentry>::iterator lowerbd_iter = sparse_index.lower_bound(-((int)v));
            int closest_vid = -((int)lowerbd_iter->first);
            assert(closest_vid>=0);
            indexentry closest_offset = lowerbd_iter->second;
            assert(closest_vid <= (int)v);
            if (closest_vid > (int)curvid) {   /* Note: this will fail if we have over 2B vertices! */
                if (curblock != NULL) // Move the pointer - this may invalidate the curblock, but it is being checked later
                    curblock->ptr += closest_offset.edataoffset - edataoffset;
                if (curadjblock != NULL)
                    curadjblock->ptr += closest_offset.adjoffset - adjoffset;
                curvid = (vid_t)closest_vid;
                adjoffset = closest_offset.adjoffset;
                edataoffset = closest_offset.edataoffset;
                return;
            } else {
                // Do nothing - just continue from current pos.
                return;
            }
            
        }
        
        inline void check_curblock(size_t toread) {
            if (curblock == NULL || curblock->end < edataoffset+toread) {
                if (curblock != NULL) {
                    if (!curblock->active) {
                        curblock->release(iomgr);
                    }
                }
                // Load next
                std::string blockfilename = filename_shard_edata_block(filename_edata, (int) (edataoffset / blocksize), blocksize);
                int edata_session = iomgr->open_session(blockfilename, false, true);
                sblock newblock(edata_session, edata_session, true);
                
                // We align blocks always to the blocksize, even if that requires
                // allocating and reading some unnecessary data.
                newblock.offset = (edataoffset / blocksize) * blocksize; // Align
                size_t correction = edataoffset - newblock.offset;
                newblock.end = std::min(edatafilesize, newblock.offset + blocksize);
                assert(newblock.end >= newblock.offset);
                iomgr->managed_malloc(edata_session, &newblock.data, newblock.end - newblock.offset, newblock.offset);
                newblock.ptr = newblock.data + correction;
                activeblocks.push_back(newblock);
                curblock = &activeblocks[activeblocks.size()-1];                
            }
        }
        
        inline void check_adjblock(size_t toread) {
            if (curadjblock == NULL || curadjblock->end <= adjoffset + toread) {
                if (curadjblock != NULL) {
                    curadjblock->release(iomgr);
                    delete curadjblock;
                    curadjblock = NULL;
                }
                sblock * newblock = new sblock(0, adjfile_session);
                newblock->offset = adjoffset;
                newblock->end = std::min(adjfilesize, adjoffset+blocksize);
                assert(newblock->end > 0);
                assert(newblock->end >= newblock->offset);
                iomgr->managed_malloc(adjfile_session, &newblock->data, newblock->end - newblock->offset, adjoffset);
                newblock->ptr = newblock->data;
                metrics_entry me = m.start_time();
                iomgr->managed_preada_now(adjfile_session, &newblock->data, newblock->end - newblock->offset, adjoffset);
                m.stop_time(me, "blockload");
                curadjblock = newblock;
            }
        }
        
        template <typename U>
        inline U read_val() {
            check_adjblock(sizeof(U));
            U res = *((U*)curadjblock->ptr);
            adjoffset += sizeof(U);
            curadjblock->ptr += sizeof(U);
            return res;
        }
        
        template <typename U>
        inline U * read_edgeptr() {
            if (only_adjacency) return NULL;
            check_curblock(sizeof(U));
            U * resptr = ((U*)curblock->ptr);
            edataoffset += sizeof(U);
            curblock->ptr += sizeof(U);
            return resptr;
        }
        
        inline void skip(int n, int sz) {
            size_t tot = n * sz;
            adjoffset += tot;
            if (curadjblock != NULL)
                curadjblock->ptr += tot;
            edataoffset += sizeof(ET)*n;
            if (curblock != NULL)
                curblock->ptr += sizeof(ET)*n;
        }
        
    public:
        /**
         * Read out-edges for vertices.
         */
        void read_next_vertices_XG(int p, int nvecs, vid_t start,  std::vector<svertex_t> & prealloc, bool record_index=false, bool disable_writes=false)  {
             metrics_entry me = m.start_time();
             /* load adj */
             char *adj;
             char *adjoffset;
             char *ed;
             filename_adj = filename_block_adj(bname, exec_interval, p, P);
             adjfile_session = iomgr->open_session(filename_adj, true);
             size_t jsize = get_filesize(filename_adj);
             adj = NULL;
             iomgr->managed_malloc(adjfile_session, &adj, jsize+1, 0);
             assert(adj != NULL);
             memset(adj,'\0', jsize + 1);
             iomgr->managed_preada_now(adjfile_session, &adj, jsize, 0);

             int nblocks,dirt;
             int esize, edatafile_session;

             if(CO == 0){ 
                   nblocks = P - 1;
             }else if(CO == 2 || CO == 3) { 
                if(interval < M) dirt = interval;
                else dirt = mini(M, P-1);
                nblocks = P - 1 - dirt;
             }else if(CO == 1){ 
                   nblocks = P -1;
             }

             bool toload = false;
             int pos = exec_interval;
             for(int i = 1; i <= nblocks; i++){
                 pos++; 
                 if(pos >= P) pos = 0;
                 if(p == pos){
                     toload = true;
                     break;
                 } 
             }
             filename_edata = filename_block_edata<ET>(bname, exec_interval, p, P, 1);
             esize = get_filesize(filename_edata);
             if(toload){ //from disk file.
                olength[p] = esize;
                edatafile_session = iomgr->open_session(filename_edata, true);
                obuf[p] = NULL;
                iomgr->managed_malloc(edatafile_session, &obuf[p], esize+1, 0);
                assert(obuf[p] != NULL);
                iomgr->managed_preada_now(edatafile_session, &obuf[p], esize, 0);
                olength[p] = esize;
            }else{
                 pthread_mutex_lock(&olock);
                 if(olength[p]==-1){
                      obuf[p] = (char*)malloc(esize+1);
                      olength[p] = esize;
                }else{
                     assert(olength[p]==esize);                                                                                                                                                   
                }   
                pthread_mutex_unlock(&olock);
            }
             

             adjoffset = adj;
             int curvid = start;
             ed = obuf[p];

             for(int i=0; i<nvecs && adjoffset!=NULL; i++) {
                int n;
                if(adjoffset-adj+sizeof(uint8_t)>jsize){
                    break;
                }
                uint8_t ns = *((uint8_t*)adjoffset);
                adjoffset += sizeof(uint8_t);

                if (ns == 0x00) {
                    if(adjoffset-adj+sizeof(uint8_t)>jsize){
                        std::cout<<"EROR:adjdata is bad, 0 : " <<adjoffset-adj+sizeof(uint8_t)<<" > "<<jsize<<std::endl;
                        assert(false);
                    }
                    curvid++;
                    uint8_t nz = *((uint8_t*)adjoffset);
                    adjoffset += sizeof(uint8_t);

                    curvid += nz;
                    i += nz;
                    continue;
                }
                
               if (ns == 0xff) {
                    if(adjoffset-adj+sizeof(uint32_t)>jsize){
                        std::cout<<"error:adjdata is bad 1 : " <<adjoffset-adj+sizeof(uint32_t)<<" > "<<jsize<<std::endl;
                        assert(false);
                    }
                    n = *((uint32_t*)adjoffset);
                    adjoffset += sizeof(uint32_t);
               } else {
                   n = ns;

               }

                svertex_t& vertex = prealloc[i];
                assert(vertex.id() == curvid);

                if (vertex.scheduled){
                    while(--n >= 0) {
                        bool special_edge = false;
                        if(adjoffset-adj+sizeof(vid_t)>jsize){
                            std::cout<<"error:adjdata is bad 2 : " <<adjoffset-adj+sizeof(vid_t)<<" > "<<jsize<<std::endl;
                            assert(false);
                        }
                        vid_t target = *((vid_t*)adjoffset);
                        adjoffset += sizeof(vid_t);

                        if(ed-obuf[p]+sizeof(ET)>esize){
                            std::cout<<"error: edata is bad!: " <<ed-obuf[p]+sizeof(ET)<<" > "<<esize<<std::endl;
                            assert(false);
                        }
                        ET * evalue = (ET*)ed;
                        ed += sizeof(ET);

                        if(target>=range_st && target<=range_end){
                            vertex.add_outedge(target, evalue, special_edge);
                            if(CO==2 || CO==3) vertex.crucial_flag=true;//prev_bitset->clear_bit(i); //for previous update
                        }

                        if (!((target >= range_st && target <= range_end))) {

                              logstream(LOG_ERROR) << "Error exec: "<<exec_interval<<"  " << target << " not in [" << range_st << " - " << range_end << "]" << std::endl;
                              iomgr->print_session(adjfile_session);
                        }
                       assert(target >= range_st && target <= range_end);
                   }    
                }else{
                    if(adjoffset-adj+n*sizeof(vid_t)>jsize){
                        std::cout<<"error:adjdata is bad 3  : " <<adjoffset-adj+n*sizeof(vid_t)<<"  > "<<jsize<<std::endl;
                        assert(false);
                    }
                    adjoffset += n * sizeof(vid_t);
                    ed += n * sizeof(ET);
                }
                    curvid++;
            }

             iomgr->managed_release(adjfile_session, &adj);
             iomgr->close_session(adjfile_session);

             if(toload) iomgr->close_session(edatafile_session);
             m.stop_time(me, "read_next_vertices");                                                                                                                 
        }

        void read_next_vertices(int nvecs, vid_t start,  std::vector<svertex_t> & prealloc, bool record_index=false, bool disable_writes=false)  {
            metrics_entry me = m.start_time();
            if (!record_index)
                move_close_to(start);
            
            /* Release the blocks we do not need anymore */
            curblock = NULL;
            release_prior_to_offset(false, disable_writes);
            assert(activeblocks.size() <= 1);
            
            /* Read next */
            if (!activeblocks.empty() && !only_adjacency) {
                curblock = &activeblocks[0];
            }
            vid_t lastrec = start;
            window_start_edataoffset = edataoffset;
            
            for(int i=((int)curvid) - ((int)start); i<nvecs; i++) {
                if (adjoffset >= adjfilesize) break;
                
                // TODO: skip unscheduled vertices.
                
                int n;
                if (record_index && (size_t)(curvid - lastrec) >= (size_t) std::max((int)100000, nvecs/16)) {
                    save_offset();
                    lastrec = curvid;
                }
                uint8_t ns = read_val<uint8_t>();
                if (ns == 0x00) {
                    curvid++;
                    uint8_t nz = read_val<uint8_t>();
                    curvid += nz;
                    i += nz;
                    continue;
                }
                
                if (ns == 0xff) {
                    n = read_val<uint32_t>();
                } else {
                    n = ns;
                }
                
                if (i<0) {
                    // Just skipping
                    skip(n, sizeof(vid_t));
                } else {
                    svertex_t& vertex = prealloc[i];
                    assert(vertex.id() == curvid);
                    
                    if (vertex.scheduled) {
                        
                        while(--n >= 0) {
                            bool special_edge = false;
                            vid_t target = (sizeof(ET) == sizeof(ETspecial) ? read_val<vid_t>() : translate_edge(read_val<vid_t>(), special_edge));
                            ET * evalue = (special_edge ? (ET*)read_edgeptr<ETspecial>(): read_edgeptr<ET>());
                            
                            if (!only_adjacency) {
                                if (!curblock->active) {
                                    if (async_edata_loading) {
                                        curblock->read_async(iomgr);
                                    } else {
                                        curblock->read_now(iomgr);
                                    }
                                }
                                // Note: this needs to be set always because curblock might change during this loop.
                                curblock->active = true; // This block has an scheduled vertex - need to commit
                            }
                            vertex.add_outedge(target, evalue, special_edge);
                            
                            if (!((target >= range_st && target <= range_end))) {
                                logstream(LOG_ERROR) << "Error : " << target << " not in [" << range_st << " - " << range_end << "]" << std::endl;
                                iomgr->print_session(adjfile_session);
                            }
                            assert(target >= range_st && target <= range_end);
                        }
                        
                    } else {
                        // This vertex was not scheduled, so we can just skip its edges.
                        skip(n, sizeof(vid_t));
                    }
                }
                curvid++;
            }
            m.stop_time(me, "read_next_vertices");
            curblock = NULL;
        }
        
        
        /**
         * Commit modifications.
         */
        void commit(sblock &b, bool synchronously, bool disable_writes=false) {
            if (disable_async_writes) synchronously = true;
            if (synchronously) {
                metrics_entry me = m.start_time();
                if (!disable_writes) b.commit_now(iomgr);
                m.stop_time(me, "commit");
                b.release(iomgr);
            } else {
                if (!disable_writes) b.commit_async(iomgr);
                else b.release(iomgr);
            }
        }
        
        /**
         * Release all buffers
         */
        void flush() {
            release_prior_to_offset(true);
            if (curadjblock != NULL) {
                curadjblock->release(iomgr);
                delete curadjblock;
                curadjblock = NULL;
            }
        }
        
        /**
         * Set the position of the sliding shard.
         */
        void set_offset(size_t newoff, vid_t _curvid, size_t edgeptr) {
            this->adjoffset = newoff;
            this->curvid = _curvid;
            this->edataoffset = edgeptr;
            if (curadjblock != NULL) {
                curadjblock->release(iomgr);
                delete curadjblock;
                curadjblock = NULL;
            }
        }
        
        /**
         * Release blocks that come prior to the current offset/
         */
        void release_prior_to_offset(bool all=false, bool disable_writes=false) { // disable writes is for the dynamic case
            for(int i=(int)activeblocks.size() - 1; i >= 0; i--) {
                sblock &b = activeblocks[i];
                if (b.end <= edataoffset || all) {
                    commit(b, all, disable_writes);
                    activeblocks.erase(activeblocks.begin() + (unsigned int)i);
                }
            }
        }
        
        void set_disable_async_writes(bool b) {
            disable_async_writes = b;
        }
        
        std::string get_info_json() {
            std::stringstream json;
            json << "\"size\": ";
            json << edatafilesize << std::endl;
            json << ", \"windowStart\": ";
            json << window_start_edataoffset;
            json << ", \"windowEnd\": ";
            json << edataoffset;
            json << ", \"intervalStart\": ";
            json << range_st;
            json << ", \"intervalEnd\": ";
            json << range_end;
            return json.str();
        }
        
    };
    
    
    
    
};


#endif
#endif

