#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "terarec.h"
#include <math.h>
#include <stdbool.h>

typedef struct{
	char key[TERA_KEY_LEN];
} terakey_t;

MPI_Datatype mpi_terakey_type;

void terakeyMPICommitType(){
	static int committed = 0;

	if(committed)
		return;

    int blocklengths[1] = {TERA_KEY_LEN};
    MPI_Datatype types[1] = {MPI_CHAR};

    MPI_Aint     offsets[1];
    offsets[0] = offsetof(terakey_t, key);

    MPI_Type_create_struct(1, blocklengths, offsets, types, &mpi_terakey_type);
    MPI_Type_commit(&mpi_terakey_type);

    committed = 1;
}

int teraKeyCompare(const void *a, const void *b){
	return memcmp( ((terakey_t*) a)->key, ((terakey_t*) b)->key, TERA_KEY_LEN);
}

int teraRecKeyCompare(const void *a, const void *b){
	return memcmp( ((terarec_t*) a)->key, ((terakey_t*) b)->key, TERA_KEY_LEN);
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

//adpted from source: https://www.geeksforgeeks.org/binary-search/
int splitter_bsearch(terarec_t *data, terakey_t *splitter, int l, int r)
{
    if (r >= l) {
        int mid = l + (r - l) / 2;
 

        bool condition_1 = (teraRecKeyCompare(data+mid, splitter) >= 0); 
        bool condition_2 = (mid ==0 || teraRecKeyCompare(data+mid-1, splitter) < 0); 
        
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



void partition_with_splitters(terarec_t *data, int  dat_len, terakey_t *splitters, int spl_len, int *bkt_counts, int *bkt_displs, int procs) {

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

void swap(terarec_t** a, terarec_t** b) {
    terarec_t* temp = *a;
    *a = *b;
    *b = temp;
}


// adapted from Algorithms: 24-part Lecture Series Robert Sedgewick / Kevin Wayne
void merge_terarec(terarec_t *data, terarec_t *aux, int lo, int mid, int hi) {

    int i = lo, j= mid+1;
    
    for (int k = lo; k <= hi; k++)    {
        if( i > mid) {
            memcpy(data+k, aux+j, sizeof(terarec_t)*(hi-k+1));
            return;
            //data[k] = aux[j++];
        }
        else if(j > hi) {
            memcpy(data+k, aux+i, sizeof(terarec_t)*(hi-k+1));
            return;
            //data[k] = aux[i++];
        }
        else if(teraCompare(aux+j, aux+i) < 0)
            data[k] = aux[j++];  
        else
            data[k] = aux[i++];
    }   
}

void merge(terarec_t **data, terarec_t **aux, int data_len, int *displs, int displs_len) {
    for (int sz = 1; sz < displs_len; sz=sz+sz)
    {
        swap(data, aux);
        for (int i = 1; i <= displs_len; i += sz+sz) {
            int lo = displs[i-1];
            int mid_idx = i+sz-1;
            int mid = (mid_idx < displs_len)? displs[mid_idx]-1 : data_len-1;
            int hi =  (mid_idx+sz < displs_len)? displs[mid_idx+sz]-1 : data_len-1;
            merge_terarec(*data, *aux, lo, mid, hi);
        }
    }

}



void exch(terarec_t *data, int i, int j) {
    terarec_t temp = data[i];
    data[i] = data[j];
    data[j] = temp;
}

void insertion_sort(terarec_t *data, int N) {
    for (int i = 0; i < N; i++)
        for (int j = i; j > 0; j--)
            if(teraCompare(data+j, data+j-1) < 0)
                exch(data, j, j-1);
            else
                break;
}

void swap_merge(terarec_t *data, terarec_t *aux, int lo, int mid, int hi) {
    
    //memcpy(aux+lo, data+lo, sizeof(terarec_t)*(hi-lo+1));

    int i = lo, j= mid+1;
    
    for (int k = lo; k <= hi; k++)    {
        if( i > mid) {
            memcpy(aux+k, data+j, sizeof(terarec_t)*(hi-k+1));
            return;
        }
        else if(j > hi) {
            memcpy(aux+k, data+i, sizeof(terarec_t)*(hi-k+1));
            return;
        }
        else if(teraCompare(data+j, data+i) < 0)
            aux[k] = data[j++];  
        else
            aux[k] = data[i++];
    }   
}



void sort_terarec(terarec_t *orig, terarec_t *data, terarec_t *aux, int lo, int hi, int cutoff) {
    if(hi <= lo+cutoff-1) {
	    //qsort(data+lo, hi-lo+1, sizeof(terarec_t), teraCompare);
        if(orig == data)
	        qsort(data+lo, hi-lo+1, sizeof(terarec_t), teraCompare);
            //insertion_sort(data+lo, hi-lo+1);
        else
	        qsort(aux+lo, hi-lo+1, sizeof(terarec_t), teraCompare);
            //insertion_sort(aux+lo, hi-lo+1);
        return;
    }
    int mid = lo + (hi-lo)/2;
    sort_terarec(orig, aux, data, lo, mid, cutoff);
    sort_terarec(orig, aux, data, mid+1, hi, cutoff);
    swap_merge(aux, data, lo, mid, hi);
}

//void sort_terarec_2(terarec_t *data, terarec_t *aux,  int *displs, int *counts, int len, int lo, int hi) {
//    if(hi <= lo) {
//        return;
//    }
//    int mid = lo + (hi-lo)/2;
//    sort_terarec_2(data, aux, displs, counts, len, lo, mid);
//    sort_terarec_2(data, aux, displs, counts, len, mid+1, hi);
//    copy_merge(data, aux, displs[lo], displs[mid]-1, displs[hi]+counts[hi]);
//}


//https://stackoverflow.com/questions/3437404/min-and-max-in-c
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// source https://stackoverflow.com/questions/8377412/ceil-function-how-can-we-implement-it-ourselves
int my_ceil(float input) {
    int output = (int)input;
    if (input == (float)output)
        return output;
    return output + 1;
}

void merge_sort(terarec_t **data, terarec_t **aux, int N, int chunk) {
    
    //double start;
    //double end;
    //int rank;
    //MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    //start = MPI_Wtime();
    int displs_len = my_ceil(N/(double)chunk); 
    int *displs = (int*) malloc(displs_len * sizeof(int));
    int offset =0;
    for (int i = 0, k=0; i < N; i+= chunk, k++)
    {
        int len = MIN(chunk, N-i);
	    //qsort(*data+i, len, sizeof(terarec_t), teraCompare);
        insertion_sort(*data+i, len);
        displs[k] = offset;
        offset += len;
    }
    //end = MPI_Wtime();
    //printf("%.6fs.| 0.  Sub Sort local data, rank %d \n", end - start, rank);

    merge(data, aux, N, displs, displs_len);
    free(displs);
    
}

void terasort(terarec_t *local_data, int  local_len, 
			  terarec_t **sorted_data, int* sorted_counts, long* sorted_displs){
    
    int rank, P;
    double start;
    double end;
    
    start = MPI_Wtime();
    MPI_Comm_size (MPI_COMM_WORLD, &P);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    terakeyMPICommitType();
    
    int s = P-1; // sample, splitters
    int d = local_len/P;
    //corner cases:
    if(s > local_len) s = local_len-1;
    if(d==0) d = 1;
    
    //terakey_t *local_samples = (terakey_t*) malloc(s*sizeof(terakey_t));
    terarec_t *original_local_data = local_data;
    terarec_t *aux1 = (terarec_t*) malloc(sizeof(terarec_t) * local_len);
    terakey_t *splitters = (terakey_t*) malloc(s*sizeof(terakey_t));
    terakey_t *key_samples = (terakey_t*) malloc(s*sizeof(terakey_t));
    terakey_t *root_samples = NULL;
    
    int *all_counts = (int*) calloc( P * P, sizeof(int));
    int *send_displs = (int*) malloc(P * sizeof(int));
    int *recv_counts = (int*) calloc(P, sizeof(int));
    int *recv_displs = (int*) malloc(P * sizeof(int));

    if(rank == 0) root_samples = (terakey_t*) malloc(P*s*sizeof(terakey_t));
    // 1. Local Sort
	//qsort(local_data, local_len, sizeof(terarec_t), teraCompare);
    //sort_terarec(original_local_data, local_data, aux1, 0, local_len-1, 10);
    merge_sort(&local_data, &aux1, local_len, 7);
    end = MPI_Wtime();
    printf("%.6fs.| 1.  Local Sort, rank %d \n", end - start, rank);

    // 2. Select P-1 samples (s)
    start = MPI_Wtime();
    //for (int i = 0; i < s; i++)
    //    local_samples[i] = local_data[(i+1)*d];
    for (int i = 0; i < s; i++) {
        memcpy(key_samples[i].key, local_data[(i+1)*d].key, TERA_KEY_LEN);
    }
    end = MPI_Wtime();
    printf("%.6fs.| 2.  Select P-1 samples, rank %d \n", end - start, rank);
    // 3. Gather samples at root 
    start = MPI_Wtime();
	MPI_Gather(key_samples, s, mpi_terakey_type, root_samples, s, mpi_terakey_type, 0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf("%.6fs.| 3.  Gather Samples at root, rank %d \n", end - start, rank);
    // 4. Sort samples at root and select splitters
    if(rank == 0) { 
        start = MPI_Wtime();
        qsort(root_samples, P*s, sizeof(terakey_t), teraKeyCompare);
        end = MPI_Wtime();
        printf("%.6fs.| 4.  Sort samples on root, rank %d\n", end - start, rank);
        start = MPI_Wtime();
        for (int i = 0, x = (P*s)/P; i < s; i++) 
            splitters[i] = root_samples[(i+1)*x];
        end = MPI_Wtime();
        printf("%.6fs.| 5.  Select P-1 splitters, rank %d\n", end - start, rank);
    }
    // 5. Broadcast splitters
    start = MPI_Wtime();
    MPI_Bcast(splitters, s, mpi_terakey_type, 0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf("%.6fs.| 6.  Broadcast splitters, rank %d\n", end - start, rank);
    // 6. Partition using splitters
    start = MPI_Wtime();
    partition_with_splitters(local_data, local_len, splitters, s, all_counts+(P*rank), send_displs, P);
    end = MPI_Wtime();
    printf("%.6fs.| 7.  Partition using splitters bsearch, rank %d\n", end - start, rank);
       
    // 7. Gather all counts
    start = MPI_Wtime();
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, all_counts, P, MPI_INT, MPI_COMM_WORLD);
    // 8. populate sort info
    populate_sort_info(all_counts, recv_counts, recv_displs, sorted_counts, sorted_displs, rank, P);
    end = MPI_Wtime();
    printf("%.6fs.| 8.  Gather all counts, rank %d\n", end - start, rank);
    
    // 9. Exchange Values
    start = MPI_Wtime();
    *sorted_data = malloc(sizeof(terarec_t) * sorted_counts[rank]);

    MPI_Alltoallv(local_data, all_counts+(P*rank), send_displs, mpi_tera_type, *sorted_data, recv_counts, recv_displs, mpi_tera_type, MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf("%.6fs.| 9.  Exchange data with all nodes, rank %d\n", end - start, rank);
    // 10. final local sort
    start = MPI_Wtime();
    terarec_t *aux2 = (terarec_t*) malloc(sizeof(terarec_t) * sorted_counts[rank]);
	//qsort(*sorted_data, sorted_counts[rank], sizeof(terarec_t), teraCompare);
    merge(sorted_data, &aux2, sorted_counts[rank], recv_displs, P);
    
    free(key_samples);
    free(splitters);
    if(rank ==0) {
        free(root_samples);
    }
    free(all_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
    free(aux2);
    if(original_local_data != local_data)
        swap(&local_data, &aux1);
    free(aux1);
    
    end = MPI_Wtime();
    printf("%.6fs.| 11. Final Sort, rank %d\n", end - start, rank);

       
}