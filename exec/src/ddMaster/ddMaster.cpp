/*
 * @file
 * @author  YongLi Cheng <sa_cyl@163.com>
 * @version 1.0
 *
 * Copyright [2014] [YongLI Cheng / HuaZhong University]
 */

#include <iostream>
#include <comm/machine.hpp>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <arpa/inet.h>
#include <assert.h>
#include <pthread.h>
#include <omp.h>
#include<queue>
#include "ddgraph_basic_includes.hpp"
using namespace xGraph;
using namespace ddgraph;

typedef float EdgeDataType;
struct topvalue{
    unsigned int vid;
    float value;
};
struct topvalue_int{
    unsigned int vid;
    int value;
};
int main(int argc, const char **argv) {
	int M;
	task *T,t1;
	char *buf;

	metrics m("Application");

	m.start_time("Iters");
	std::vector<std::pair<int, std::string> > hosts;
	std::vector<task> Tc;
	std::cout<<"Today,begin...,2014/3/24"<<std::endl;

	ddgraph_init(argc, argv);
	/* determine the number of workers */
    	int machines = get_option_int("machines", 20); //get_option_int
	M = determine_machines("./conf/hosts.slave", machines, &hosts);
	if(M < 1){
		std::cout<<"There is not a worker. "<<std::endl;
		assert(false);
	}
	for(unsigned int i=0;i<=hosts.size()-1;i++)
            std::cout<<"worker"<<i<<": "<<hosts[i].second.c_str()<<std::endl;

    	/* Process input file - if not already preprocessed */
    	std::string filename = get_option_string("file"); // Base filename  src/util/cmdopts.hpp
    	int niters = get_option_int("niters", 4); //get_option_int
    	int crucial_option = get_option_int("crucial_option", 1);
        if(crucial_option<0 || crucial_option>3){
            std::cout<<"0<= CO <= 3" <<std::endl;
            assert(false);
        }
        int ntop = get_option_int("top", 20);
	if(crucial_option == -1) std::cout<<"crucial_option to be required: 0 adj only; 1 send inedge; 2 send outedge; 3 send inedge and outedge"<<std::endl;
	assert(crucial_option != -1);
    	int nshards = get_option_int("nshards", -1);
        if(nshards == -1){
            std::cout<<"Please input nshards......" <<std::endl;
            assert(nshards != -1);
        }

	if(nshards < M) M = nshards; //in-memory computing

	/* Assign tasks for workers*/
	T = (task*)malloc(sizeof(task) * nshards * niters);
	int k = 0;
	for(int i = 0; i < niters; i++) {
		for(int j = 0; j < nshards; j++) {
			if((k % M) == 0) k = 0;
			t1.machine = k;
			t1.iter = i;
			t1.interval = nshards * i + j;
			T[i * nshards + j] = t1;
			k ++;
		}
	}

	for(int i = 0; i < M; i++) {
		Tc.clear();
		for(int j = 0; j < niters * nshards; j ++) {
			if(T[j].machine == i) Tc.push_back(T[j]);
		}
		int len = 6 + 6*sizeof(int) + sizeof(task)*Tc.size()+15*M;
		buf = (char*)malloc(len);
		memcpy(buf,"XGREQS",6); //package 'S'=schedule
		
		int c;
		c = Tc.size();
		memcpy(buf + 6,(void*)&c,sizeof(int)); //package 'S'=schedule
		memcpy(buf + 6 + sizeof(int),(void*)&M,sizeof(int)); 
		memcpy(buf + 6 + 2*sizeof(int),(void*)&niters,sizeof(int));
		memcpy(buf + 6 + 3*sizeof(int),(void*)&nshards,sizeof(int));
		memcpy(buf + 6 + 4*sizeof(int),(void*)&crucial_option,sizeof(int));
		memcpy(buf + 6 + 5*sizeof(int),(void*)&ntop,sizeof(int));
		
		char ip[15];
		for(int k = 0; k < M; k++) {
			memset(ip, 0, 15);
			memcpy(ip, hosts[k].second.c_str(), 15);	
			memcpy((void*)(buf+6+6*sizeof(int)+k*15),ip,15);	
		}
		for(unsigned int j = 0; j < Tc.size(); j++) memcpy((void*)(buf+6+6*sizeof(int)+15*M+j*sizeof(task)),(void*)&(Tc[j]),sizeof(task));
		for(unsigned int j = 0; j < Tc.size(); j++) std::cout<<((task*)(buf+6+6*sizeof(int)+15*M+j*sizeof(task)))->machine<<" " \
                                                        <<((task*)(buf+6+6*sizeof(int)+15*M+j*sizeof(task)))->iter<<" " \
							<<((task*)(buf+6+6*sizeof(int)+15*M+j*sizeof(task)))->interval<<" " \
							<<(((task*)(buf+6+6*sizeof(int)+15*M+j*sizeof(task)))->interval)%nshards<<std::endl;
		int wret=write(hosts[i].first ,buf, len);
		if(wret != len){
			std::cout<<"Send Schedule Message to "<<hosts[i].second<< "Fail!"<<std::endl;
		}
		int recvbytes;
                if ((recvbytes = recv(hosts[i].first, buf, 1024,0)) == -1) {
                        std::cout<<"recv error!"<<std::endl;
                }
		assert(recvbytes != -1);
		char Phead[7];
		memcpy(Phead, buf, 6);
        	if(strcmp(Phead,"XGACKS")!=0){
                	std::cout<<"Recv Package S Fail:"<<hosts[i].second<<std::endl;
			assert(false);
        	}
		std::cout<<*((int*)(buf+6))<<Phead<<std::endl;
	}

	/*
	 * Section: Begin to Schedule.
     	 */
        std::cout<<"Tasks have been launched!"<<std::endl;

	k = niters * nshards;
	int ret;
	if(crucial_option != 0){
		for(int i = 0; i < k; i ++) {

			while((ret = recv(hosts[T[i].machine].first, buf, 6,0))<6){
				std::cout<< " Waiting for Worker Signal..."<<ret<<std::endl;
                exit(0);
				if(ret < 0){
					std::cout<<"Recv fail, exit..."<<std::endl;
					exit(0);
				}
			}

			if(i < k -1){
				while((ret = write(hosts[T[i+1].machine].first ,(char*)"XGREQN", 6)) < 6){

					std::cout<< " Sending signal to worker..."<<std::endl;
					if(ret < 0){
						std::cout<<"Recv fail, exit..."<<std::endl;
						exit(0);
					}
				}

			}
			std::cout<<"Finished...."<<i<<" / "<<k-1<<std::endl;
		}
	}else {
        if(M < nshards){
		    for(int i = 0; i < k / M; i ++){
			    for(int j = 0; j < M; j++){
				    while(write(hosts[j].first ,(char*)"XGREQN", 6)<6){};
			    }	
			    volatile int cc = 0;
			    omp_set_num_threads(M);
        	    #pragma omp parallel for
			    for(int t = 0; t < M; t ++){
				    while(recv(hosts[t].first, buf, 1024,0) < 0){
				    }
					    std::cout<<"Finished...."<< i*M + t<<" / "<<k-1<<std::endl;
				    __sync_add_and_fetch(&cc, 1);
			    }
			    while(cc < M){};
		    }
        }else if(M == nshards){
		    volatile int cc = 0;
		    omp_set_num_threads(M);
            #pragma omp parallel for
            for(int y = 0; y < M; y++){
                for(int t = 0; t < k/M; t++){
                    write(hosts[y].first ,(char*)"XGREQN", 6); 
                    ret = recv(hosts[y].first, buf, 1024,0);
                    if(ret > 0) std::cout<<"Finished...."<< P*t + y<<" / "<<k-1<<std::endl;
                }
                 __sync_add_and_fetch(&cc, 1);
            } 
            while(cc<M) {};
        }
	}


	/* Merge Results */
    std::queue<topvalue> q[50];
    std::queue<topvalue_int> q1[50];
    int vt;
	volatile int done = 0;
	    omp_set_num_threads(M);
        #pragma omp parallel for
		for(int i = 0; i < M; i ++) {
            char *topbuf;
            topbuf = (char*)malloc(1024*1024*10+1);
			recv(hosts[i].first, topbuf,1024*1024*10, 0);
            topvalue tv;
            topvalue_int tv_int;
            memcpy(&vt,topbuf,sizeof(int));
            for(int j = 0; j < ntop; j++){
                if(vt==0){
                    memcpy(&tv.vid, topbuf+(sizeof(int)+sizeof(float))*j+sizeof(int), sizeof(int));
                    memcpy(&tv.value, topbuf+(sizeof(int)+sizeof(float))*j+2*sizeof(int), sizeof(float));
                    q[i].push(tv);
                }else if(vt==1){
                    memcpy(&tv_int.vid, topbuf+(sizeof(int)+sizeof(int))*j+sizeof(int), sizeof(int));
                    memcpy(&tv_int.value, topbuf+(sizeof(int)+sizeof(int))*j+2*sizeof(int), sizeof(int));
                    //std::cout<<"vid: "<<tv_int.vid<<"  value:  "<<tv_int.value<<std::endl;
                    q1[i].push(tv_int);
                }
            }
            __sync_add_and_fetch(&done, 1); 
		}
	while(done < M){}	
    /*
    for(int v=0;v<M;v++){
        std::cout<<"============================== M:"<<v<<std::endl;
        for(int w=0;w<ntop;w++){
            topvalue_int tt;
            tt=q1[v].front();
            std::cout<<" vid: " <<tt.vid<<" value: "<<tt.value<<std::endl;
            q1[v].pop();
        } 
    }
    */
    int pos,f;
    unsigned int oldvid;
    float value;
    int   val;
    std::stringstream fname;

    std::cout << "Print top " << ntop << " vertices: "<< std::endl;
    fname << filename <<".map";
    f = open(((std::string)fname.str()).c_str(), O_RDONLY);
    if(f < 0){
        std::cout<<"Can not open the file : "<<fname.str()<<"   toplist fail! "<< strerror(errno)<<std::endl; 
    }else{
        if(vt==0){
            for(int t=0; t<ntop; t++){
                pos=0;
                value = (float)0.0;
                for(int p=0; p<M; p++){
                    if(q[p].front().value>value){
                        value = q[p].front().value;
                        pos = p;
                    } 
                }
                if(pread(f, &oldvid, sizeof(int),q[pos].front().vid*sizeof(int))<0){
                    std::cout<<"read oldvid fail!  " << strerror(errno)<<std::endl;;
                    break;
                }
                std::cout<<t+1<<".  "<<oldvid<<"  "<<q[pos].front().value<<std::endl;
                q[pos].pop();
            }
        }else if(vt==1){
            for(int t=0; t<ntop; t++){
                pos=0;
                val = 0;
                for(int p=0; p<M; p++){
                    if(q1[p].front().value>val){
                        val = q1[p].front().value;
                        pos = p;
                    } 
                }
                if(pread(f, &oldvid, sizeof(int),q1[pos].front().vid*sizeof(int))<0){
                    std::cout<<"read oldvid fail!  " << strerror(errno)<<std::endl;;
                    break;
                }
                std::cout<<t+1<<".  "<<oldvid<<"  "<<q1[pos].front().value<<std::endl;
                q1[pos].pop();
            }
        }
        close(f);
    }
	std::cout<<"Finished"<<std::endl;
	m.stop_time("Iters");

	/* inform workers to finish. */
       for(int i = 0; i < M; i ++) {
            send(hosts[i].first, (char*)"XGREQF", 6, 0);
        }


	/* close clients connect.*/
        for(unsigned int i=0;i<=hosts.size()-1;i++)
            close(hosts[i].first);

	/* Report Result */
	metrics_report(m); 
	return(0);
}
