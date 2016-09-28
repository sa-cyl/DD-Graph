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
 * Simple pagerank implementation. Uses the basic vertex-based API for
 * demonstration purposes. A faster implementation uses the functional API,
 * "pagerank_functional".
 */

#include <string>
#include <fstream>
#include <cmath>

#define GRAPHCHI_DISABLE_COMPRESSION


#include "ddgraph_basic_includes.hpp"

using namespace ddgraph;
 
#define THRESHOLD 1e-1    
#define RANDOMRESETPROB 0.15


typedef float VertexDataType;
typedef float EdgeDataType;
//for a test
int c;
struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, ddgraph_context &info) {
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, ddgraph_context &ginfo) {
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
    void before_exec_interval(vid_t window_st, vid_t window_en, ddgraph_context &ginfo) {        
    }
    
    
    /**
      * Pagerank update function.
      */
    void update(ddgraph_vertex<VertexDataType, EdgeDataType> &v, ddgraph_context &ginfo) {
        float sum=0;
	//for a test
	//int have;
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
            for(int i=0; i < v.num_outedges(); i++) {
                ddgraph_edge<float> * edge = v.outedge(i);
                int v1=v.vertexid;
                edge->set_data(1.0 / v.num_outedges());
            }
            v.set_data(RANDOMRESETPROB); 
        } else {
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
	    //have=0;
            for(int i=0; i < v.num_inedges(); i++) {
                float val = v.inedge(i)->get_data();
                sum += val;
				//for a test
				//if(((v.inedge(i)->vertexid<ginfo.range_st) || (v.inedge(i)->vertexid>ginfo.range_en))&&(ginfo.iteration==1)) have=1;
				//for a test
            }
           //if(have==1) c++; 
	   //if((c%1000)==0)std::cout<<"Crucial node number:======================"<<c<<std::endl;
            /* Compute my pagerank */
            float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            /* Write my pagerank divided by the number of out-edges to
               each of my out-edges. */
            if (v.num_outedges() > 0) {
                float pagerankcont = pagerank / v.num_outedges();
                for(int i=0; i < v.num_outedges(); i++) {
                    ddgraph_edge<float> * edge = v.outedge(i);
                    edge->set_data(pagerankcont);
                }
            }
                
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
            ginfo.log_change(std::abs(pagerank - v.get_data()));
            
            /* Set my new pagerank as the vertex value */
            v.set_data(pagerank); 
        }
    }
    
};
int main(int argc, const char ** argv) {
    ddgraph_init(argc, argv);
    metrics m("pagerank");
    
    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 4);
    bool scheduler          = false;                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 20);
    
    /* Process input file - if not already preprocessed */
    int nshards             = get_option_int("nshards", -1);
    if(nshards == -1){
        std::cout<<"Please input nshards...." <<std::endl;
        assert(nshards != -1);
    }

    /* Run */
    ddgraph_engine<float, float> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    PagerankProgram program;
    engine.run(program, niters);
        
    /* Output top ranked vertices */
    /*
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
   */ 
    metrics_report(m);    
    return 0;
}

