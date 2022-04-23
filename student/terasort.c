#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "terarec.h"
#include <math.h>
#include <stdbool.h>

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

void populate_sort_info(int *all_counts, int *recv_counts, int *recv_displs, int *sorted_counts, long *sorted_displs, int my_rank, int procs) {

    sorted_displs[0] = 0;
    recv_displs[0] = 0;

    for (int r = 0; r < procs; r++)
    {
        //calculate recv_counts and recv_displs, contribution from all 'r' ranks to current rank (me)
        recv_counts[r] = all_counts[r*procs + my_rank];
        if(r >0 )
            recv_displs[r] = recv_displs[r-1] + recv_counts[r-1];
            
        //calculate sorted_counts and sorted displs, contributin from all ranks to any rank
        for (int j = 0; j < procs; j++) {
            if(j==0)
                sorted_counts[r] = 0;
            sorted_counts[r] += all_counts[j*procs + r];
        }
        if(r >0 )
            sorted_displs[r] = sorted_displs[r-1] + sorted_counts[r-1];

    }
}

void terasort_bucket(terarec_t *local_data, int  local_len, 
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
    
    populate_sort_info(bucket_counts, recv_counts, recv_displs, sorted_counts, sorted_displs, rank, P);

    *sorted_data = malloc(sizeof(terarec_t) * sorted_counts[rank]);
    
    MPI_Alltoallv(xdata, bucket_counts+(P*rank), send_displs, mpi_tera_type, *sorted_data, recv_counts, recv_displs, mpi_tera_type, MPI_COMM_WORLD);

	qsort(*sorted_data, sorted_counts[rank], sizeof(terarec_t), teraCompare);

    free(xdata);
    free(bucket_ids);
    free(bucket_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
}


/*
void partition_with_splitters(terarec_t *data, int  dat_len, terarec_t *splitters, int spl_len, int *bkt_counts, int *bkt_displs, int procs) {

    int k = 0; // track the splitter num
    bkt_displs[0] = 0;

    for (int i = 0; i < dat_len; i++) {
        while(teraCompare(data+i, splitters+k) >= 0 && k < spl_len) {
            ++k;
            bkt_displs[k] = i;
        }
        ++bkt_counts[k];
    }
}
*/







//adpted from source: https://www.geeksforgeeks.org/binary-search/
int splitter_bsearch(terarec_t *data, terarec_t *splitter, int l, int r)
{
    if (r >= l) {
        int mid = l + (r - l) / 2;
 

        bool condition_1 = (teraCompare(data+mid, splitter) >= 0); 
        bool condition_2 = (mid ==0 || teraCompare(data+mid-1, splitter) < 0); 
        
        //splitter found
        if(condition_1 && condition_2)
            return mid;
        //splitter could only be present in left sub data
        if (condition_1)
            return splitter_bsearch(data, splitter, l, mid - 1);
        // Else the splitter could only be present in right sub data
        return splitter_bsearch(data, splitter, mid + 1, r);
    }
    // splitter not present
    return -1;
}



void partition_with_splitters(terarec_t *data, int  dat_len, terarec_t *splitters, int spl_len, int *bkt_counts, int *bkt_displs, int procs) {

    int k = 0; // track the splitter num
    int found_index;
    int next_displ = dat_len;
    bkt_displs[0] = 0;

    for (int i = 0; i < spl_len; i++) {
        found_index =  splitter_bsearch(data, &splitters[i], 0, dat_len);
        ++k;
        bkt_displs[k] = found_index ;
    }

    for (int i = spl_len; i >= 0; i--)
    {
        if(bkt_displs[i] == -1)
            continue;
        bkt_counts[i] = next_displ - bkt_displs[i];
        next_displ = bkt_displs[i];
    }
    
}

void terasort(terarec_t *local_data, int  local_len, 
			  terarec_t **sorted_data, int* sorted_counts, long* sorted_displs){
    
    int rank, P;
    double start;
    double end;
    
    MPI_Comm_size (MPI_COMM_WORLD, &P);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    int s = P-1; // sample, splitters
    int d = local_len/P;
    //corner cases:
    if(s > local_len) s = local_len-1;
    if(d==0) d = 1;
    
    terarec_t *local_samples = (terarec_t*) malloc(s*sizeof(terarec_t));
    terarec_t *splitters = (terarec_t*) malloc(s*sizeof(terarec_t));
    terarec_t *root_samples = NULL;
    int *all_counts = (int*) calloc( P * P, sizeof(int));
    int *send_displs = (int*) malloc(P * sizeof(int));
    int *recv_counts = (int*) calloc(P, sizeof(int));
    int *recv_displs = (int*) malloc(P * sizeof(int));

    if(rank == 0) root_samples = (terarec_t*) malloc(P*s*sizeof(terarec_t));
    // 1. Local Sort
    start = MPI_Wtime();
	qsort(local_data, local_len, sizeof(terarec_t), teraCompare);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 1.  Local Sort\n", end - start);

    // 2. Select P-1 samples (s)
    start = MPI_Wtime();
    for (int i = 0; i < s; i++)
        local_samples[i] = local_data[(i+1)*d];
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 2.  Select P-1 samples\n", end - start);
    // 3. Gather samples at root 
    start = MPI_Wtime();
	MPI_Gather(local_samples, s, mpi_tera_type, root_samples, s, mpi_tera_type, 0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 3.  Gather Samples at root\n", end - start);
    // 4. Sort samples at root and select splitters
    if(rank == 0) { 
        start = MPI_Wtime();
        qsort(root_samples, P*s, sizeof(terarec_t), teraCompare);
        end = MPI_Wtime();
        if(rank ==0) printf("%.6fs.| 4.  Sort samples on root\n", end - start);
        start = MPI_Wtime();
        for (int i = 0, x = (P*s)/P; i < s; i++) 
            splitters[i] = root_samples[(i+1)*x];
        end = MPI_Wtime();
        if(rank ==0) printf("%.6fs.| 5.  Select P-1 splitters\n", end - start);
    }
    // 5. Broadcast splitters
    start = MPI_Wtime();
    MPI_Bcast(splitters, s, mpi_tera_type, 0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 6.  Broadcast splitters\n", end - start);
    // 6. Partition using splitters
    start = MPI_Wtime();
    partition_with_splitters(local_data, local_len, splitters, s, all_counts+(P*rank), send_displs, P);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 7.  Partition using splitters bsearch\n", end - start);
    
    /*
    start = MPI_Wtime();
    partition_with_splitters(local_data, local_len, splitters, s, all_counts+(P*rank), send_displs, P);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 7.  Partition using splitters old\n", end - start);
    
    if(rank==0) {
    printf("\nPartiions old\n");
    for (size_t i = 0; i < P; i++)
    {
        printf("rank %d displ %d counts %d\n", rank, send_displs[i], all_counts[P*rank+i]); 
    }
    }
    */


    
    // 7. Gather all counts
    start = MPI_Wtime();
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, all_counts, P, MPI_INT, MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 8.  Gather all counts\n", end - start);
    // 8. populate sort info
    populate_sort_info(all_counts, recv_counts, recv_displs, sorted_counts, sorted_displs, rank, P);
    // 9. Exchange Values
    *sorted_data = malloc(sizeof(terarec_t) * sorted_counts[rank]);

    start = MPI_Wtime();
    MPI_Alltoallv(local_data, all_counts+(P*rank), send_displs, mpi_tera_type, *sorted_data, recv_counts, recv_displs, mpi_tera_type, MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 9.  Exchange data with all nodes\n", end - start);
    // 10. final local sort
    start = MPI_Wtime();
	qsort(*sorted_data, sorted_counts[rank], sizeof(terarec_t), teraCompare);
    end = MPI_Wtime();
    if(rank ==0) printf("%.6fs.| 11. Final Sort\n", end - start);

/*
    sleep(0.5*rank);
    printf("\nLocal Keys \n");
    for (int i = 0; i < local_len; i++)
    {
        printf("rank %d key %s \n", rank, local_data[i].key);
    }
    
    if(rank==0) {
    printf("\nRoot Samples\n");
    for (size_t i = 0; i < P*s; i++)
    {
        printf("root %s\n", root_samples[i].key); 
    }
    }
    
    printf("\nAll Counts\n");
    for (int i = 0; i < P*P; i++)
    {
       printf("rank %d counts %d\n", rank, all_counts[i]);
    }
    
    printf("\nSplitters\n");
    for (size_t i = 0; i < s; i++)
    {
        printf("rank %d split %s\n", rank, splitters[i].key); 
    }
    printf("\nTo Send\n");
    for (size_t i = 0; i < P; i++)
    {
        printf("rank %d displ %d counts %d\n", rank, send_displs[i], all_counts[P*rank+i]); 
    }
    
    printf("\nSorted keys disp %ld counts %d \n", sorted_displs[rank], sorted_counts[rank]);
    for (int i = 0; i < sorted_counts[rank]; i++)
    {
        printf("rank %d sorted key %s \n", rank, (*sorted_data)[i].key);
    }
    printf("\n--------------\n");
*/
    
    free(local_samples);
    free(splitters);
    if(rank ==0) {
        free(root_samples);
    }
    free(all_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);

    
    
}