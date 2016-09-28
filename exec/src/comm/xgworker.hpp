/* 
 * This file was created by YongLI Cheng 2014/3/16 \
 */ 
#include "io/stripedio.hpp"
#include <iostream>
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

using namespace ddgraph;

#define BACKLOG 10
#define PORT1 3333
#define PORT2 7777
#define MAXSIZE 1024 * 1024 * 50
#define COMMTHREADS 20

static int exec_interval;
static int interval;
static int ll_main_to_main;
static pthread_mutex_t tlock,olock;
static volatile int fd_master = -1;
static volatile int inbuf_loaded = 0;
static volatile int obuf_loaded = 0;
static int M, P, N,CO,sockfd1,sockfd2;
static volatile int interval_dirt_disk[1024];
static volatile int interval_dirt_memory[1024];
static volatile int interval_dirt_main[1024];
static volatile int exec_interval_buf = -1;
static char  **obuf;
static volatile int olength[1024];
static char **inbuf;
static volatile int inlength[1024];
static std::string bname;
static char *in_to_out_buf = NULL;
static char *out_to_in_buf = NULL;
static char *main_to_main_buf = NULL;
static int cout_g=0;
static dense_bitset *prev_bitset = NULL;
/* for test */
static int cc_prev = 0;

char *wbuf;
int read_all(int fd, void *buf, int n);
int write_all(int fd, void *buf, int n);
int write_to_worker(char *ip, int port,  void *buf, int from , int to, int flage, int n);

void *reply(void *client_fd){ //because this function is static ,its vars is static
	
	int fd,rval;
	char msg[1024];
	int Bid, offset, length;
	
	fd = *(int*)client_fd;
	rval = read_all(fd, msg, sizeof(int)*3);
	if(rval == -1)  std::cout<<"Recv Initial data...fail..."<<std::endl;
	memcpy(&Bid, msg, sizeof(int));
	memcpy(&offset, msg + sizeof(int), sizeof(int));
	memcpy(&length, msg + 2*sizeof(int), sizeof(int));
	//std::cout<<"Bid: "<<Bid<<"=================== offset: "<<offset<<" length: "<<length<<" fd: "<<fd<<std::endl;

        rval = read_all(fd, wbuf + offset, length);
        send(fd,(char*)"XGACKT", 6, 0);
        close(fd);
}
	 void *server1(void *info) {
		int ret, len;
		char msg[1024];
		char buf[1024 + 1];
		struct sockaddr_in my_addr;
		struct sockaddr_in remote_addr;
		int from, to, flag, length;
        char *tb = NULL;

		if ((sockfd1 = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
			std::cout<<"socket create failed!"<<std::endl;
		}
		assert(sockfd1!=-1);

		my_addr.sin_family = AF_INET;
		my_addr.sin_port = htons(PORT1);
		my_addr.sin_addr.s_addr = INADDR_ANY;
		bzero(&(my_addr.sin_zero), 8);
		if (bind(sockfd1, (struct sockaddr*) &my_addr, sizeof(struct sockaddr))==-1){
			std::cout<<"bind error!"<<std::endl;
			assert(false);
		}

		if (listen(sockfd1, BACKLOG) == -1) {
			std::cout<<"listen error!"<<std::endl;
			assert(false);
		}
		int isFirst = 1;
		int cc =0;
		int client_fd;
		
		while (1) {
			cc++;
			socklen_t sin_size = sizeof(struct sockaddr_in);
			if ((client_fd = accept(sockfd1, (struct sockaddr*) &remote_addr, &sin_size)) == -1) {
				std::cout<<"accept error!"<<std::endl;
				continue;
			}
			if(isFirst == 1) {//for xgMaster.
				isFirst = 0;
				std::cout<<"Received a connection from"<<inet_ntoa(remote_addr.sin_addr)<<std::endl;
				read(client_fd, buf, 1024);
            			char* msg1 = (char*)"XGACKT";
            			if (send(client_fd, msg1, strlen(msg1), 0) == -1)
                    			std::cout<<"send back XGACKS fail!"<<std::endl;
				fd_master = client_fd;
				continue;
			}

			/* todo for otherworker */

			/* for no_multi_threads */
			ret = read_all(client_fd, msg, sizeof(int)*4);
			if(ret == -1)  std::cout<<"Recv Initial data...fail..."<<std::endl;
			memcpy(&from, msg, sizeof(int));
			memcpy(&to, msg + sizeof(int), sizeof(int));
			memcpy(&flag, msg + 2*sizeof(int), sizeof(int));
			memcpy(&length, msg + 3*sizeof(int), sizeof(int));
			
            //std::cout<<"recv: length "<<length<<" from "<<from<<" to "<<to<<std::endl;
            int f = from % P;
            int t = to % P;
            std::string block_filename;
            if(tb != NULL) free(tb);
            tb = (char*)malloc(length+1);
            //std::cout<<"recv:.....1  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
            ret = read_all(client_fd, tb, length);
            //std::cout<<"recv:.....2  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
			send(client_fd,(char*)"XGACKT", 6, 0);
			if(client_fd > 0) close(client_fd);
            //std::cout<<"recv:.....3  length "<<length<<" from "<<from<<" to "<<to<<std::endl;

            while(exec_interval_buf == -1){} //for first interval
			if(to - from < M && t == exec_interval_buf){
                if(flag == 0){
				    if(to > P * N -1){ 
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 0);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        //std::cout<<"recv:..process...1  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                        writeedatablock(bname, tb, t, f , P, 0,length);
                        std::cout<<"recv:..process...1 end  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                    }else{
                        //std::cout<<"recv:..process...2  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                        __sync_add_and_fetch(&cout_g, 1);
                        pthread_mutex_lock(&tlock);
                        if(inlength[f]==-1){
                            if(inbuf[f] != NULL) free(inbuf[f]);
                            inbuf[f] = (char*)malloc(length+1);
                            memcpy(inbuf[f], tb, length);
                            inlength[f] = length;
                        }else if(inlength[f]>0){
                            assert(inlength[f]==length);
                            memcpy(inbuf[f], tb, length);
                        }
                        pthread_mutex_unlock(&tlock);
                        //std::cout<<"recv:..process...2 end length "<<length<<" from "<<from<<" to "<<to<<std::endl;

                    }
                }else if(flag == 1){
				    if(to > P * N -1){
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 1);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        writeedatablock(bname, tb, t, f , P, 1,length);
                    }else{
                        __sync_add_and_fetch(&cout_g, 1);
                        pthread_mutex_lock(&olock);
                        if(olength[f]==-1){
                            if(obuf[f] != NULL) free(obuf[f]);
                            obuf[f] = (char*)malloc(length+1);
                            memcpy(obuf[f],tb,length);
                            olength[f] = length;
                        }else if(olength[f]>0){
                            assert(olength[f]==length);
                            memcpy(obuf[f], tb, length);
                        }
                        pthread_mutex_unlock(&olock);
                    }
                }
				//interval_dirt_memory[to % P] --;
				__sync_add_and_fetch(&interval_dirt_memory[to % P], -1);
			} else if(to - from > M || t != exec_interval_buf){
                  if(flag == 0){
                        //check data length
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 0);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        //std::cout<<"recv:..process...3  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                        writeedatablock(bname, tb, t, f , P, 0,length);
                        //std::cout<<"recv:..process...3 end length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                  }else if(flag == 1){
                        //check data length
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 1);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        std::cout<<"recv:..process...4  length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                        __sync_add_and_fetch(&cout_g, 1);
                        writeedatablock(bname, tb, t, f , P, 1,length);
                        std::cout<<"recv:..process...4 end length "<<length<<" from "<<from<<" to "<<to<<std::endl;
                  }
                   __sync_add_and_fetch(&interval_dirt_disk[to % P], -1);
			}	
		}
	}

	void *server2(void *info) {
		int sockfd;
		int ret, len;
		char msg[1024];
		char buf[1024 + 1];
        char *tb = NULL;
		struct sockaddr_in my_addr;
		struct sockaddr_in remote_addr;
		int from, to, flag, length;

		if ((sockfd2 = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
			std::cout<<"socket create failed!"<<std::endl;
		}
		assert(sockfd2!=-1);

		my_addr.sin_family = AF_INET;
		my_addr.sin_port = htons(PORT2);
		my_addr.sin_addr.s_addr = INADDR_ANY;
		bzero(&(my_addr.sin_zero), 8);
		if (bind(sockfd2, (struct sockaddr*) &my_addr, sizeof(struct sockaddr))==-1){
			std::cout<<"bind error!"<<std::endl;
			assert(false);
		}

		if (listen(sockfd2, BACKLOG) == -1) {
			std::cout<<"listen error!"<<std::endl;
			assert(false);
		}
		int client_fd;
		
		while (1) {
			 socklen_t sin_size = sizeof(struct sockaddr_in);
			 if ((client_fd = accept(sockfd2, (struct sockaddr*) &remote_addr, &sin_size)) == -1) {
				std::cout<<"accept error!"<<std::endl;
				continue;
			 }


            /* for single_thread */
            ret = read_all(client_fd, msg, sizeof(int)*4);
            if(ret == -1)  std::cout<<"Recv Initial data...fail..."<<std::endl;
            memcpy(&from, msg, sizeof(int));
            memcpy(&to, msg + sizeof(int), sizeof(int));
            memcpy(&flag, msg + 2*sizeof(int), sizeof(int));
            memcpy(&length, msg + 3*sizeof(int), sizeof(int));

            int f = from % P;
            int t = to % P;
            std::string block_filename;
            if(tb != NULL) free(tb);
            tb = (char*)malloc(length+1);
            ret = read_all(client_fd, tb, length);
			send(client_fd,(char*)"XGACKT", 6, 0);
			if(client_fd > 0) close(client_fd);

            while(exec_interval_buf == -1){} //for first interval
			if(to - from < M && t == exec_interval_buf){
                if(flag == 0){
				    if(to > P * N -1){ 
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 0);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        writeedatablock(bname, tb, t, f , P, 0,length);
                    }else{
                        __sync_add_and_fetch(&cout_g, 1);
                        pthread_mutex_lock(&tlock);
                        if(inlength[f]==-1){
                            if(inbuf[f] != NULL) free(inbuf[f]);
                            inbuf[f] = (char*)malloc(length+1);
                            memcpy(inbuf[f], tb, length);
                            inlength[f] = length;
                        }else if(inlength[f]>0){
                            assert(inlength[f]==length);
                            memcpy(inbuf[f], tb, length);
                        }
                        pthread_mutex_unlock(&tlock);
                    }
                }else if(flag == 1){
				    if(to > P * N -1){
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 1);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        writeedatablock(bname, tb, t, f , P, 1,length);
                    }else{
                        __sync_add_and_fetch(&cout_g, 1);
                        pthread_mutex_lock(&olock);
                        if(olength[f]==-1){
                            if(obuf[f] != NULL) free(obuf[f]);
                            obuf[f] = (char*)malloc(length+1);
                            memcpy(obuf[f],tb,length);
                            olength[f] = length;
                        }else if(olength[f]>0){
                            assert(olength[f]==length);
                            memcpy(obuf[f], tb, length);
                        }
                        pthread_mutex_unlock(&olock);
                    }
                }
				//interval_dirt_memory[to % P] --;
				__sync_add_and_fetch(&interval_dirt_memory[to % P], -1);
			} else if(to - from > M || t != exec_interval_buf){
                  if(flag == 0){
                        //check data length
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 0);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        writeedatablock(bname, tb, t, f , P, 0,length);
                  }else if(flag == 1){
                        //check data length
                        /*
                        block_filename = filename_block_edata<float>(bname, to%P, from%P, P, 1);
                        if (file_exists(block_filename)){                                                                                                                                 
                             size_t fsize = get_filesize(block_filename);
                             if(fsize != length){
                                std::cout<<"edata is bad!"<<std::endl;
                             }    
                             assert(fsize==length);
                         }else{
                             std::cout<<"can not find file: " << block_filename <<std::endl;
                             assert(false);
                        }
                        */
                        __sync_add_and_fetch(&cout_g, 1);
                        writeedatablock(bname, tb, t, f , P, 1,length);
                  }
                   __sync_add_and_fetch(&interval_dirt_disk[to % P], -1);
             }
		}
	}
int connect_worker(char *ip, int i) {
        struct sockaddr_in serv_addr;
	
	int sockfd[COMMTHREADS];
	int c = 0;
        while (( sockfd[i] = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
               std::cout<<"socket error!"<<std::endl;
               return(-1);
        }

        bzero(&serv_addr,sizeof(serv_addr));
        serv_addr.sin_family    = AF_INET;
        serv_addr.sin_port      = htons(PORT1);
        serv_addr.sin_addr.s_addr= inet_addr(ip);

        if(connect(sockfd[i], (struct sockaddr *)&serv_addr,sizeof(struct sockaddr)) == -1) {
              std::cout<<"connect error!"<<errno<<std::endl;
              return(-1);
        }
	return(sockfd[i]);
}

int read_all(int fd, void *buf, int n) {
	int nleft;
	int nbytes;
	char *ptr;
	ptr = (char*)buf;
	nleft = n;
	for(; nleft > 0;){
		nbytes = read(fd, ptr, nleft);
		if(nbytes < 0){
			if(errno == EINTR) nbytes = 0;
			else return(-1);
		} else if(nbytes == 0) break;
		nleft -= nbytes;
		ptr += nbytes;
	}
	return(n - nleft);
}

int write_all(int fd, void *buf, int n) {
	int nleft, nbytes;
	char *ptr;
	nleft = n;
	ptr = (char*)buf;
	for(; nleft > 0;){
		nbytes = write(fd, ptr, nleft);
		if(nbytes <= 0){
			if(errno == EINTR) nbytes = 0;
			else return(-1);
		}
		nleft -= nbytes;
		ptr += nbytes;
  	}
	return(n);
}


int write_to_worker_multi(char *ip, int port, void *buf, int n) {
	int fd[COMMTHREADS],cc;
	char msg[COMMTHREADS][1024];
	int ret[COMMTHREADS];
        struct sockaddr_in serv_addr;

	int Bid[COMMTHREADS];   	//the id of sliding block
	int offset[COMMTHREADS]; 	//offset of the sliding block
	int length[COMMTHREADS]; 	//sent length chars

        timeval start_time,end;
        double lasttime;
        gettimeofday(&start_time, NULL);

	cc = n / COMMTHREADS;
        volatile int done = 0;
	omp_set_num_threads(COMMTHREADS);
	#pragma omp parallel for
	for(int i = 0; i < COMMTHREADS; i ++) {
		Bid[i] = 0;
		offset[i] = cc * i;
		length[i] = cc;
		if((i == COMMTHREADS -1) && ((n % COMMTHREADS)>0)) length[i] = n % COMMTHREADS;

	//std::cout<<Bid[i]<<"||===New=======||"<<offset[i]<<"||"<<length[i]<<"   threadid : "<<pthread_self()<<std::endl;
      	/*fd[i] = connect_worker(ip, i);

	if(fd[i] == -1){
           	std::cout<<"Connect to " << ip << " Fail!"<<std::endl;
		//if(fd[i] > 0) close(fd[i]);
           	assert(false);
      	}*/

	/* start:connect to worker. */

        if((fd[i] = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
               std::cout<<"socket error!"<<std::endl;
	       assert(false);
        }

        bzero(&serv_addr,sizeof(serv_addr));
        serv_addr.sin_family    = AF_INET;
        serv_addr.sin_port      = htons(PORT1);
        serv_addr.sin_addr.s_addr= inet_addr(ip);

        if(connect(fd[i], (struct sockaddr *)&serv_addr,sizeof(struct sockaddr)) == -1) {
              std::cout<<"connect error!"<<errno<<std::endl;
	      assert(false);
        }
	/* end:connect to worker. */

	//std::cout<<"|"<<fd[i]<<"|"<<std::endl;
	//std::cout<<"i : "<<i<<"|||"<<Bid[i]<<"||||"<<offset[i]<<"||"<<length[i]<<std::endl;
	memcpy(msg[i], (void*)&Bid[i], sizeof(int));
	memcpy(msg[i] + sizeof(int) , (void*)&offset[i], sizeof(int));
	memcpy(msg[i] + 2*sizeof(int), (void*)&length[i], sizeof(int));
        memcpy(&Bid[i], msg[i], sizeof(int));
        memcpy(&offset[i], msg[i] + sizeof(int), sizeof(int));
        memcpy(&length[i], msg[i] + 2*sizeof(int), sizeof(int));
	//std::cout<<"i : "<<i<<"|||"<<Bid[i]<<"||||"<<offset[i]<<"||"<<length[i]<<std::endl;
	ret[i] = write_all(fd[i], msg[i], sizeof(int)*3);
	if(ret[i] == -1) {
                 std::cout<<"send Intial data fail!"<<std::endl;
                 if(fd[i] > 0) close(fd[i]);
                 assert(false);
	}

        /*if (send(fd[i], msg[i], sizeof(int)*3,0) == -1){
                 std::cout<<"send Intial data fail!"<<std::endl;
		 if(fd[i] > 0) close(fd[i]);
		 assert(false);
	}*/

	ret[i] = write_all(fd[i], (void*)((char*)buf + offset[i]), length[i]);
	if(ret[i] == -1){
		std::cout<<"send to "<<ip<<" fail!"<<std::endl;
                if(fd[i] > 0) close(fd[i]);
		assert(false);
	} else {
		recv(fd[i], msg[i], 1024,0);
        	gettimeofday(&end, NULL);
                lasttime = end.tv_sec - start_time.tv_sec + ((double)(end.tv_usec - start_time.tv_usec)) / 1.0E6;
                //std::cout<<"send time:  "<<lasttime<<std::endl;
		__sync_add_and_fetch(&done, 1);
		close(fd[i]);
	}

	}
	while(done < COMMTHREADS){}
	return(0);
}


int write_to_worker(char *ip, int port, void *buf, int from, int to, int flag, int n) { //flag =0 inedge, flag=1 outedge
	int fd;
	char msg[1024];
	int ret;
        struct sockaddr_in serv_addr;

	int Bid;   	//the id of sliding block
	int offset; 	//offset of the sliding block
	int length; 	//sent length chars

        timeval start_time,end;
        double lasttime;
        gettimeofday(&start_time, NULL);

	Bid = 0;
	offset = 0;
	length = n;


	/* start:connect to worker. */

        if((fd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
               std::cout<<"socket error!"<<std::endl;
	       assert(false);
        }

        bzero(&serv_addr,sizeof(serv_addr));
        serv_addr.sin_family    = AF_INET;
        serv_addr.sin_port      = htons(port);
        serv_addr.sin_addr.s_addr= inet_addr(ip);

        if(connect(fd, (struct sockaddr *)&serv_addr,sizeof(struct sockaddr)) == -1) {
              std::cout<<"connect error!"<<errno<<std::endl;
	      assert(false);
        }
	/* end:connect to worker. */

	//std::cout<<"|"<<fd<<"|"<<std::endl;
	//std::cout<<"i : "<<i<<"|||"<<Bid[i]<<"||||"<<offset[i]<<"||"<<length[i]<<std::endl;
	memcpy(msg, (void*)&from, sizeof(int));
	memcpy(msg + sizeof(int) , (void*)&to, sizeof(int));
	memcpy(msg + 2*sizeof(int) , (void*)&flag, sizeof(int));
	memcpy(msg + 3*sizeof(int), (void*)&length, sizeof(int));
        //memcpy(&Bid, msg, sizeof(int));
        //memcpy(&offset, msg + sizeof(int), sizeof(int));
        //memcpy(&length, msg + 2*sizeof(int), sizeof(int));
	//std::cout<<"i : "<<i<<"|||"<<Bid[i]<<"||||"<<offset[i]<<"||"<<length[i]<<std::endl;
	ret = write_all(fd, msg, sizeof(int)*4);
	if(ret == -1) {
                 std::cout<<"send Intial data fail!"<<std::endl;
                 if(fd > 0) close(fd);
                 assert(false);
	}

	ret = write_all(fd, (void*)((char*)buf), length);
	if(ret == -1){
		std::cout<<"send to "<<ip<<" fail!"<<std::endl;
                if(fd > 0) close(fd);
		assert(false);
	} else {
		recv(fd, msg, 1024,0);
        	gettimeofday(&end, NULL);
                lasttime = end.tv_sec - start_time.tv_sec + ((double)(end.tv_usec - start_time.tv_usec)) / 1.0E6;
                //std::cout<<"send time:  "<<lasttime<<std::endl;
		close(fd);
	}

	return(0);
}
