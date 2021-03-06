/* 
 * This file was created by YongLI Cheng 2014/3/16 \
 */ 
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
#include <schedule/task.hpp>

#define SERVPORT 3333
#define BACKLOG 10
#define MAXSIZE 1024*1024*sizeof(int)+6

namespace xGraph {

struct RC {
	runcontrol *R;
	int *fd;
};
static void *reply(void *rc){
	int rval;
	int *li;
	int fd=(*(RC*)rc).fd;
	int I=(*(RC*)rc).intervals;
	task T=(*(RC*)rc).T;

	char buf[MAXSIZE];
	char Phead[6];
	char Pack_t[1];
	int offset=0;
	unsigned int nbyte;

	if ((rval = read(fd, buf, MAXSIZE)) < 0) {
		std::cout<<"reading stream error!"<<std::endl;
		close(fd);
		pthread_cancel(pthread_self());
	}
	nbyte=rval;
	li=(int*)(buf+6);
	memcpy(Phead,buf,5);
	if(strcmp(Phead,"XGREQ")!=0){
		std::cout<<"Bad Package Formationt|"<<Phead<<std::endl;
		close(fd);
		pthread_cancel(pthread_self());
	}

	Pack_t[0]=buf[5];
	if(Pack_t[0]=='A'){
		std::cout<<"Aaaaa"<<Pack_t<<"|"<<Phead<<"|"<<li[0]<<"|"<<li[1]<<std::endl;

		char* msg = (char*)"Hello,Mr hqlong, you are connected!\n";
		if (send(fd, msg, strlen(msg), 0) == -1)
			std::cout<<"send error!"<<std::endl;
		close(fd);
	}else if(Pack_t[0]=='B'){//'B':to Begin
		int i;
		while(nbyte<MAXSIZE){
			offset+=rval;
			rval=read(fd, buf+offset, MAXSIZE);
			nbyte+=rval;
			std::cout<<rval<<";;;;;;;;;;;;;;;;;;;;"<<std::endl;			
		}
		for(i=0;i<1024*1024;i++) {
			std::cout<<li[i]<<"|"<<std::endl;
		}
		//std::cout<<"Aaaaa"<<Pack_t<<"|"<<Phead<<"|"<<li[1048570]<<"|"<<i-1<<std::endl;
		char* msg = (char*)"Hello,Mr hqlong, you are connected!\n";
		if (send(fd, msg, strlen(msg), 0) == -1)
			std::cout<<"send error!"<<std::endl;
		close(fd);
/*	}else if(Pack_t[0]=='U'){//'F':A interval have finished.
		if(run.windowpos < run.intervals){
			send(fd, (char*)"XGACKN", 6, 0)	
		}
		runc.windowpos += 1;	
	}
*/
	}else{
		std::cout<<"Unkown Package Formation"<<"|"<<Phead<<std::endl;
		//close(fd);
		//pthread_cancel(pthread_self());
	}
	return(NULL);
}

	static void *server(void *info) {
		int sockfd, client_fd;
		RC rc;
		rc.R = (runcontrol*)info;
		struct sockaddr_in my_addr;
		struct sockaddr_in remote_addr;

		if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
			std::cout<<"socket create failed!"<<std::endl;
		}
		assert(sockfd!=-1);

		my_addr.sin_family = AF_INET;
		my_addr.sin_port = htons(SERVPORT);
		my_addr.sin_addr.s_addr = INADDR_ANY;
		bzero(&(my_addr.sin_zero), 8);
		if (bind(sockfd, (struct sockaddr*) &my_addr, sizeof(struct sockaddr))==-1){
			std::cout<<"bind error!"<<std::endl;
			assert(false);
		}

		if (listen(sockfd, BACKLOG) == -1) {
			std::cout<<"listen error!"<<std::endl;
			assert(false);
		}

		while (1) {
			socklen_t sin_size = sizeof(struct sockaddr_in);
			if ((client_fd = accept(sockfd, (struct sockaddr*) &remote_addr, &sin_size)) == -1) {
				std::cout<<"accept error!"<<std::endl;
				continue;
			}
			std::cout<<"Received a connection from"<<inet_ntoa(remote_addr.sin_addr)<<std::endl;

                	pthread_t t;
			rc.fd = client_fd;
                	int ret = pthread_create(&t, NULL, &reply, &rc);
			assert(ret>=0);
		}
	}
};
