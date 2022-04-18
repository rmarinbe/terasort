#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "terarec.h"
#include <math.h>

#define BUCKET_MIN 32032
#define BUCKET_MAX 80126
#define BUCKET_RANGE (BUCKET_MAX - BUCKET_MIN)

void terasort_orig(terarec_t *local_data, int  local_len, 
			  terarec_t **sorted_data, int* sorted_counts, long* sorted_displs){
	//The code below gathers the data together at a single root node, 
	//this root node sorts the data and then scatters it back to 
	//the others.  Replace this naive strategy with a better one in 
	//in the code below.  Your code will be evaluated on its correctness
	//and its performance.

	int rank, P;
	MPI_Comm_size (MPI_COMM_WORLD, &P);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	terarec_t* all_data = NULL;

	if(rank == 0){
		all_data = (terarec_t*) malloc(sizeof(terarec_t) * P * local_len);
	}

	MPI_Gather(local_data, local_len, mpi_tera_type,
			   all_data, local_len, mpi_tera_type, 0, MPI_COMM_WORLD);

	if(rank == 0){
		qsort(all_data, P * local_len, sizeof(terarec_t), teraCompare);
	}

	*sorted_data = malloc(sizeof(terarec_t) * local_len);

	MPI_Scatter(all_data, local_len, mpi_tera_type, 
				*sorted_data, local_len, mpi_tera_type, 
				0, MPI_COMM_WORLD);

	for(int i = 0; i < P; i++)
		sorted_counts[i] = local_len;
	sorted_displs[0] = 0;
	for(int i = 1; i < P; i++)
		sorted_displs[i] = sorted_displs[i-1] + sorted_counts[i-1];

	if(rank == 0)
		free(all_data);
}

int get_bucket(int bucket_len, char key[11]) {
    int key_val =key[0]*1000 + key[1];
    return ((key_val-BUCKET_MIN)/bucket_len);
}

void local_scan(terarec_t *local_data, int  local_len, int *bucket_ids, int *bucket_counts, int *send_displs, int procs) {
    int bkt_len = BUCKET_RANGE/procs;
    for (int i = 0; i < local_len; i++)
    {
        int bkt_num = get_bucket(bkt_len, local_data[i].key);
        bucket_ids[i] = bkt_num;
        bucket_counts[bkt_num] +=1;
    }
    
    send_displs[0] = 0;
    for(int i = 1; i < procs; i++)
    	send_displs[i] = send_displs[i-1] + bucket_counts[i-1];
    
}

void prepare_xdata(terarec_t *local_data, int  local_len, int *bucket_ids, int *send_displs, terarec_t *xdata, int procs) {

    int *elem_idx = (int*) calloc(procs, sizeof(int));
    
    for (int i = 0; i < local_len; i++)
    {
        int bkt = bucket_ids[i];
        xdata[send_displs[bkt] + elem_idx[bkt]]= local_data[i];
        ++elem_idx[bkt];
    }
    free(elem_idx);
}

void calc_counts_and_displs(int *bucket_counts, int *recv_counts, int *recv_displs, int *sorted_counts, long *sorted_displs, int rank, int procs) {

    sorted_displs[0] = 0;
    recv_displs[0] = 0;

    for (int i = 0; i < procs; i++)
    {
        //calculate recv_counts and recv_displs
        recv_counts[i] = bucket_counts[i*procs+rank];
        if(i >0 )
            recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            
        //calculate sorted_counts and sorted displs
        for (int j = 0; j < procs; j++) {
            if(j==0)
                sorted_counts[i] = 0;
            sorted_counts[i] += bucket_counts[j*procs+i];
        }
        if(i >0 )
            sorted_displs[i] = sorted_displs[i-1] + sorted_counts[i-1];

    }
}

void terasort(terarec_t *local_data, int  local_len, 
			  terarec_t **sorted_data, int* sorted_counts, long* sorted_displs){
	//The code below gathers the data together at a single root node, 
	//this root node sorts the data and then scatters it back to 
	//the others.  Replace this naive strategy with a better one in 
	//in the code below.  Your code will be evaluated on its correctness
	//and its performance.

    int rank, P;
    MPI_Request request;
    MPI_Comm_size (MPI_COMM_WORLD, &P);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    terarec_t *xdata = (terarec_t*) malloc(local_len*sizeof(terarec_t));
    int *bucket_ids = (int*) malloc(local_len*sizeof(int));
    int *bucket_counts = (int*) calloc( P*P, sizeof(int));
    int *recv_counts = (int*) calloc(P, sizeof(int));
    int *send_displs = (int*) malloc(P *sizeof(int));
    int *recv_displs = (int*) malloc(P *sizeof(int));
    
    local_scan(local_data, local_len, bucket_ids, bucket_counts+(P*rank), send_displs, P);

    MPI_Iallgather(MPI_IN_PLACE, 0, MPI_INT, bucket_counts, P, MPI_INT, MPI_COMM_WORLD, &request);
    
    prepare_xdata(local_data, local_len, bucket_ids, send_displs, xdata, P);
    
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    
    calc_counts_and_displs(bucket_counts, recv_counts, recv_displs, sorted_counts, sorted_displs, rank, P);

    *sorted_data = malloc(sizeof(terarec_t) * sorted_counts[rank]);
    
    MPI_Alltoallv(xdata, bucket_counts+(P*rank), send_displs, mpi_tera_type, *sorted_data, recv_counts, recv_displs, mpi_tera_type, MPI_COMM_WORLD);

	qsort(*sorted_data, sorted_counts[rank], sizeof(terarec_t), teraCompare);

/*    
    printf("my rank: %d\n", rank);
    for (int i = 0; i < P; i++)
    {
        printf("rank: %d | ", i);
        for (int j = 0; j < P; j++)
        {
            printf("%d, ", bucket_counts[i*P+j]);
        }
        printf("\n");
        
    }
    for (size_t i = 0; i < P; i++)
    {
        printf(" %d, %d;", send_displs[i], bucket_counts[P*rank+i]);
    }
    
    printf("\n");
*/

/*
    sleep(0.5*rank);
    for (int i = 0; i < sorted_counts[rank]; i++)
    {
        printf("rank: %d , key %s\n", rank, (*sorted_data)[i].key);
    }
*/

    free(xdata);
    free(bucket_ids);
    free(bucket_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
}