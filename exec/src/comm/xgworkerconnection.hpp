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
#define BACKLOG 10
#define PORT 3333
#define MAXSIZE 1024 * 1024 * 15

static volatile int fd_master = -1;
char *wbuf;
int read_all(int fd, void *buf, int n);
int write_all(int fd, void *buf, int n);
int write_to_worker(char *ip, void *buf, int n);

static void *reply(void *client_fd){
	int fd = *(int*)client_fd;
        int rval;
	char msg[2024];
	int Bid; //which sliding Block to be writed.
	int offset; 
	int length;

	//while(read(fd, msg, 3*sizeof(int)) < 0){
	while(recv(fd, msg, 1024,0) < 0){
		std::cout<<"Recv Initial data..."<<std::endl;
	}
	memcpy(&Bid, msg, sizeof(int));
	memcpy(&offset, msg + sizeof(int), sizeof(int));
	memcpy(&length, msg + 2*sizeof(int), sizeof(int));
	std::cout<<"Bid: "<<Bid<<" offset: "<<offset<<" length: "<<length<<std::endl;

        rval = read_all(fd, wbuf + offset, length);
        std::cout<<"{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{"<<rval<<std::endl;
        send(fd,(char*)"XGACKT", 6, 0);
        if(fd > 0) close(fd);
}
	static void *server(void *info) {
		int sockfd, client_fd;
		int ret, len;
		char *buf1;
		char buf[1024 + 1];
		struct sockaddr_in my_addr;
		struct sockaddr_in remote_addr;

		if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
			std::cout<<"socket create failed!"<<std::endl;
		}
		assert(sockfd!=-1);

		my_addr.sin_family = AF_INET;
		my_addr.sin_port = htons(PORT);
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
		buf1 = (char*)malloc(MAXSIZE +1);
		int isFirst = 1;
		while (1) {
			socklen_t sin_size = sizeof(struct sockaddr_in);
			if ((client_fd = accept(sockfd, (struct sockaddr*) &remote_addr, &sin_size)) == -1) {
				std::cout<<"accept error!"<<std::endl;
				continue;
			}
			if(isFirst == 1) {//for xgMaster.
				isFirst = 0;
				std::cout<<"Received a connection from"<<inet_ntoa(remote_addr.sin_addr)<<std::endl;
				read(client_fd, buf, 1024);
            			char* msg = (char*)"XGACKT";
            			if (send(client_fd, msg, strlen(msg), 0) == -1)
                    			std::cout<<"send back XGACKS fail!"<<std::endl;
				fd_master = client_fd;
				continue;
			}

			/* todo for otherworker */

			std::cout<<"Received a connection from"<<inet_ntoa(remote_addr.sin_addr)<<client_fd<<std::endl;
/*
			int nbyte = 0;
			int offset = 0;
			int rval;
			rval = read_all(client_fd, buf1, MAXSIZE);
			std::cout<<"{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{"<<rval<<std::endl;
			send(client_fd,(char*)"XGACKT", 6, 0);
			if(client_fd > 0) close(client_fd);
*/

                	pthread_t t;
			wbuf = buf1;
                	int ret = pthread_create(&t, NULL, &reply, &client_fd);
			assert(ret>=0);
		}
	}

int connect_worker(char *ip) {
	int sockfd;
        struct sockaddr_in serv_addr;

        if (( sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
               std::cout<<"socket error!"<<std::endl;
               return(-1);
        }

        bzero(&serv_addr,sizeof(serv_addr));
        serv_addr.sin_family    = AF_INET;
        serv_addr.sin_port      = htons(PORT);
        serv_addr.sin_addr.s_addr= inet_addr(ip);

        if(connect(sockfd, (struct sockaddr *)&serv_addr,sizeof(struct sockaddr)) == -1) {
              std::cout<<"connect error!"<<std::endl;
              return(-1);
        }
	return(sockfd);
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
			if(errno = EINTR) nbytes = 0;
			else return(-1);
		}
		nleft -= nbytes;
		ptr += nbytes;
  	}
	return(n);
}

int write_to_worker(char *ip, void *buf, int n) {
	int fd,cc;
	char msg[1024];
	int ret;
	cc = 0;

	int Bid;   	//the id of sliding block
	int offset; 	//offset of the sliding block
	int length; 	//sent length chars

        timeval start_time,end;
        double lasttime;
        gettimeofday(&start_time, NULL);
      	while((fd = connect_worker(ip)) == -1){
           	cc ++;
           	std::cout<<"Connect to " << ip << " Fail!"<<std::endl;
		if(fd > 0) close(fd);
           	if(cc == 3) assert(false);
      	}

	Bid = 0;
	offset = 0;
	length = n;
	memcpy(msg, (void*)&Bid, sizeof(int));
	memcpy(msg + sizeof(int) , (void*)&offset, sizeof(int));
	memcpy(msg + 2*sizeof(int), (void*)&length, sizeof(int));
        memcpy(&Bid, msg, sizeof(int));
        memcpy(&offset, msg + sizeof(int), sizeof(int));
        memcpy(&length, msg + 2*sizeof(int), sizeof(int));
	std::cout<<Bid<<"||||"<<offset<<"||"<<length<<std::endl;
        //if (write(fd, msg, sizeof(int)*4) == -1){
        if (send(fd, msg, sizeof(int)*3,0) == -1){
                 std::cout<<"send Intial data fail!"<<std::endl;
		 if(fd > 0) close(fd);
		 assert(false);
	}

	ret = write_all(fd, buf, n);
	if(ret == -1){
		std::cout<<"send to "<<ip<<" fail!"<<std::endl;
		assert(false);
	} else {
		recv(fd, msg, 1024,0);
        	gettimeofday(&end, NULL);
                lasttime = end.tv_sec - start_time.tv_sec + ((double)(end.tv_usec - start_time.tv_usec)) / 1.0E6;
                std::cout<<"send time:  "<<lasttime<<std::endl;
		close(fd);
	 	return(ret);
	}
}
