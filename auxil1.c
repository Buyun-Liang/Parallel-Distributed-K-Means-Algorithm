#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define MAX_LINE 210

/*-------------------- begin prototyping */
int saxpy_(int *n, float *a, float *x, int *incx, float *y, int *incy);
int sscal_(int *n, float *alpha, float *x, int *inc);
void scopy_(int *nfeat, float *x, int *incx, float *y, int *incy);
void get_rand_ftr(float *ctr, float *fdata, int m, int nfeat);
/*-------------------- end prototyping */

void get_rand_ftr(float *ctr, float *fdata, int m, int nfeat){
// gets a random convex  combination of all samples
  float tot, t; 
  int j, one=1;
  /*-------------------- initialize to zero */
  for (j=0; j<nfeat; j++)
    ctr[j] = 0.0;
  tot = 0.0;
  /*-------------------- loop over all samples*/ 
  for (j=0; j<m; j++) {
    t = (float)(rand() / (float) RAND_MAX);
    t = t*t;
    //if (t < 0.5){
    //    for (k=0; k<nfeat; k++)      ctr[k] += t*fdata[j*nfeat+k];
    saxpy_(&nfeat,&t,&fdata[j*nfeat],&one,ctr,&one);
    tot +=t;    
  }
  tot = 1.0/tot;
  sscal_(&nfeat,&tot,ctr,&one);
}

int assign_ctrs(float *dist, int k){
  float min;
  min = dist[0];
  int i, ctr = 0;
  for(i=1; i<k; i++) {
    if(min > dist[i]) {
      ctr = i;
      min = dist[i];
    }
  }
  return ctr;
}

/*-------------------- reading data */
int read_csv_matrix(float *mtrx, char file_name[], int* nrow, int *nfeat){
/* -------------------- reads data from a csv file to mtrx */
  FILE *finputs; 
  char line[MAX_LINE], subline[MAX_LINE];
  
  if (NULL == (finputs = fopen(file_name, "r" )))
    exit(1);
  memset(line,0,MAX_LINE);
  //  
  int k, j, start, rlen, lrow=0, lfeat=0, jcol=0, jcol0=0, len, first= 1;
  char *delim;
  delim =",";
  /*-------------------- big while loop */
  while(fgets(line,MAX_LINE,finputs)) {
    if(first) {
//--------------------ignore first line of csv file 
      first = 0;
      continue;
    }
    len = strlen(line);
    lrow++;
    start = 0;
/*-------------------- go through the line */
    for (j=0; j<len; j++){
      if (line[j] == *delim || j==len-1) {
	k = j-start;
	memcpy(subline,&line[start],k*sizeof(char));
	//-------------------- select items to drop here --*/
	if (start >0){     //  SKIPPING THE FIRST RECORD 
	  subline[k] = '\0';
	  mtrx[jcol++] = atof(subline);
	}
	start = j+1;
      }
    }
/*-------------------- next row */
    rlen  = jcol -jcol0;
    jcol0 = jcol;
    if (lrow == 1) lfeat = rlen;
/*-------------------- inconsistent rows */
    if (rlen != lfeat)  return(1);
  }
/*-------------------- done */
  fclose(finputs);
  //  for (j=0; j<jcol; j++)    printf(" %e  \n",mtrx[j]);
  *nrow = lrow;
  *nfeat = lfeat;
  return(0);
}
/*-------------------- assign a center to an item */

float dist2(float *x, float *y, int len){
  int i;
  float dist = 0;
  for(i=0; i<len; i++){
    dist += (x[i] - y[i]) * (x[i] - y[i]);
  }
  return dist/len;
}

/*=======================================================================*/

int MyKmeans_p(float *fdata, int *clustId, int *counter, int *params,
	       float tol, MPI_Comm comm) {
/*==================================================
  IN: 
    fdata     = float* (nfeat*m)   = input data (local to this process).
    params[ ]  = int* contains 
    params[0] = Nc     = number of clusters
    params[1] = m      =  number of samples
    params[2] = nfeat  = number of features
    params[3] = maxits = max number of Kmeans iterations
    tol       = tolerance for determining if centers have converged.
    comm      = communicator
  OUT:
   clustId[i] = cluster of sample i for i=1...m
   counter[j] = size of cluster j for j=1:Nc
   ===================================================*/
  // some declarations
  int nprocs, myid;
  /*-------------------- START */
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&myid);
  //-------------------- unpack params.  
  int Nc = params[0];
  int m  = params[1]; 
  int nfeat = params[2]; 
  int maxits = params[3];
  int NcNf =Nc*nfeat;
  /*-------------------- replace these by your function*/   
  
  /* ------------------- initialize centers - random means*/
  
  float ctrs[NcNf];
  for (int i = 0; i < Nc; i++) {
	get_rand_ftr( &ctrs[i*nfeat], fdata, m, nfeat);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, ctrs, NcNf, MPI_FLOAT,MPI_SUM, MPI_COMM_WORLD);
  
  // obtain the average location
  for (int ii = 0; ii < Nc; ii++) {
  	for (int jj = 0; jj < nfeat; jj++) {
		ctrs[ii * nfeat + jj] /= nprocs;
	}
  }	

  /* ------------------- main loop  */

  // clustId[i] = cluster of sample i for i=1...m
  int i;
  //set tol be negative to ensure main loop ends at max iterations; 
  //can be change to positive value
  tol = -0.1;
  
  for (i = 0; i < maxits; i++) {
	/* ------- Get Closest Centers Ctrs */
	  
	//counter[j] = size of cluster j for j=1:Nc 
  	for (int ii = 0; ii < Nc; ii++) {
		counter[ii] = 0;
	}
	float pt_sums[NcNf];

	// initialize pt_sums for every iteration
	for (int ii = 0; ii < Nc; ii++) {
        	for (int jj = 0; jj < nfeat; jj++) {
                	pt_sums[ii * nfeat + jj] = 0;
        	}
  	}

	for (int j = 0; j < m; j++) {
		float dists[Nc];
		float vt[nfeat];
		for (int jj = 0; jj < nfeat; jj++) {
			vt[jj] = fdata[j*nfeat + jj];
		}
		for (int k = 0; k < Nc; k++) {
			dists[k] = dist2( &ctrs[k*nfeat], vt, nfeat);
		}
		/* -------------------- closest ctr to V(i,:) */
		int kk = assign_ctrs(dists, Nc);

		clustId[j] = kk;
		counter[kk] = counter[kk] + 1;
		
		for (int l = 0; l < nfeat; l++) {
			pt_sums[kk * nfeat + l] += fdata[j*nfeat + l]; 	
		}
	}

	// communication with other PEs
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, counter, Nc, MPI_INT,MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, pt_sums, NcNf, MPI_FLOAT,MPI_SUM, MPI_COMM_WORLD);
	/* -------------------- second half: reset centroids */
 
	/* if Not_conv indicates not converged yet */
	int not_conv = 0;
	if (myid == 0) {
		for (int ii = 0; ii < Nc; ii++) {
			if (counter[ii] > 0) {
				float new_c[nfeat];
				for (int jj = 0; jj < nfeat; jj++) {
					new_c[jj] = pt_sums[ii * nfeat + jj]/counter[ii];
				}
				/* -------------------- this is a convergence test */
				//printf(" myid %d dist = %.10f \n",myid,dist2(&ctrs[ii*nfeat], new_c, nfeat));

				if ( dist2(&ctrs[ii*nfeat], new_c, nfeat) > tol  ) {
					for (int jj = 0; jj < nfeat; jj++) {
						ctrs[ii*nfeat + jj] = new_c[jj];
					}
					not_conv = 1;
				}
			}
			else {
				/* -------------------- cluster is empty - reset centriod  */
				get_rand_ftr( &ctrs[ii*nfeat], fdata, m, nfeat);
				not_conv = 1;	
			}
		}
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&not_conv, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(ctrs, NcNf, MPI_FLOAT, 0, MPI_COMM_WORLD);

	if (not_conv == 0) {
		printf(" break loop: myid %d iteration = %d  \n",myid,i);
		return(0);
	}

  }
  printf(" myid %d iteration = %d  \n",myid,i);
  return(0);
}

