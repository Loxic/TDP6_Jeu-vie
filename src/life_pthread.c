#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
#include <semaphore.h>

#define NBR_THREADS 8

int BS;

int nb_threads=NBR_THREADS;
int maxloop;

int *_board;
int *_nbngb;
sem_t *nbdone;

pthread_cond_t barrier_cond;
pthread_mutex_t barrier_mut;
int barrier_id=0;
int barrier=0;

#define cell( _i_, _j_ ) board[ ldboard * (_j_) + (_i_) ]
#define ngb( _i_, _j_ )  nbngb[ ldnbngb * ((_j_) - 1) + ((_i_) - 1 ) ]

void barrier(){
  pthread_mutex_lock(&barrier_mut);
  int id = barrier_id;
  barrier++;
  if(barrier == nb_threads){
    barrier = 0;
    barrier_id++;
    pthread_cond_broadcast(&barrier_cond);
  }
  while(id ==barrier_id){
    pthread_cond_wait(&barrier_cond, &barrier_mut);
  }
  pthread_mutex_unlock(&barrier_mut);
}

double mytimer(void)
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

void output_board(int N, int *board, int ldboard, int loop)
{
    int i,j;
    printf("loop %d\n", loop);
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    if ( cell( i, j ) == 1)
		printf("X");
	    else
		printf(".");
	}
	printf("\n");
    }
}

/**
 * This function generates the iniatl board with one row and one
 * column of living cells in the middle of the board
 */
int generate_initial_board(int N, int *board, int ldboard)
{
    int i, j, num_alive = 0;

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    if (i == N/2 || j == N/2) {
		cell(i, j) = 1;
		num_alive ++;
	    }
	    else {
		cell(i, j) = 0;
	    }
	}
    }

    return num_alive;
}

void * thread_compute(void *arg){

  int tid = (int)arg;

  int subBSi=BS/nb_threads;
  if(tid==nb_threads-1)
    subBSi+= BS%nb_threads;
  int subBSj=BS;

  int ldboard = BS+2;
  int ldnbngb = BS;
  int *board = _board + tid*subBSi;
  int *nbngb = _nbngb + tid*subBSi;

    for (loop = 1; loop <= maxloop; loop++) {
      
	for(int j=1; j<=subBSj; j++){
	  ngb(1, j) = cell(0, j-1) + cell(1, j-1) + cell(2, j-1) +
	              cell(0,   j) +                cell(2,   j) +
	              cell(0, j+1) + cell(1, j+1) + cell(2, j+1);
	  ngb(subBSi, j) = cell(subBSi-1, j-1) + cell(subBSi, j-1) + cell(subBSi+1, j-1) +
	                   cell(subBSi-1,   j)                     + cell(subBSi+1,   j) +
	                   cell(subBSi-1, j+1) + cell(subBSi, j+1) + cell(subBSi+1, j+1);
	}
	sem_post(nbdone+(tid-1)%nb_threads);
	sem_post(nbdone+(tid+1)%nb_threads);

	for (j = 1; j <= subBSj; j++) {
	  for (i = 2; i <= subBSi-1; i++) {
		ngb( i, j ) =
		    cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
		    cell( i-1, j   ) +                  cell( i+1, j   ) +
		    cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
	    }
	}

	num_alive = 0;

	for (j = 1; j <= subBSj; j++) {
	    for (i = 2; i <= subBSi-1; i++) {
		if ( (ngb( i, j ) < 2) ||
		     (ngb( i, j ) > 3) ) {
		    cell(i, j) = 0;
		}
		else {
		    if ((ngb( i, j )) == 3)
			cell(i, j) = 1;
		}
		if (cell(i, j) == 1) {
		    num_alive ++;
		}
	    }
	}

	sem_wait(nbdone+tid);
	sem_wait(nbdone+tid);

	for(int j=1; j<=subBSj; j++){
	  switch (ngb(1, j)){
	  case 3:
	    cell(1, j) = 1;
	  case 2:
	    break;
	  default:
	    cell(1, j) = 0;
	  }
	  if(cell(1,j)==1)
	    num_alive++;

	  switch (ngb(subBSi, j)){
	  case 3:
	    cell(subBSi, j) = 1;
	  case 2:
	    break;
	  default:
	    cell(subBSi, j) = 0;
	  }
	  if(cell(subBSi,j)==1)
	    num_alive++; 
	}

	for(int i=1; i<=subBSi; i++){
	  cell(i, 0) = cell(i, subBSj);
	  cell(i, subBSj+1) = cell(i, 1);
	}

	barrier();

	if(tid == 0){
	  cell(   0, 0   ) = cell(BS, BS);
	  cell(   0, BS+1) = cell(BS,  1);
	  cell(BS+1, 0   ) = cell( 1, BS);
	  cell(BS+1, BS+1) = cell( 1,  1);	  
	  for(int j=1; j<=subBSj; j++){
	    cell(   0, j) = cell(BS, j);
	    cell(BS+1, j) = cell(1, j);
	  }
	}
	
	barrier();
    }
}

int main(int argc, char* argv[])
{
    int i, j, loop, num_alive;
    int ldboard, ldnbngb;
    double t1, t2;
    double temps;

    if (argc < 3) {
	printf("Usage: %s nb_iterations size [nb_threads]\n", argv[0]);
	return EXIT_SUCCESS;
    } else {
	maxloop = atoi(argv[1]);
	BS = atoi(argv[2]);
	//printf("Running sequential version, grid of size %d, %d iterations\n", BS, maxloop);
    }
    if(argc > 3)
      nb_threads = atoi(argv[3]);
    num_alive = 0;

    /* Leading dimension of the board array */
    ldboard = BS + 2;
    /* Leading dimension of the neigbour counters array */
    ldnbngb = BS;

    _board = malloc( ldboard * ldboard * sizeof(int) );
    _nbngb = malloc( ldnbngb * ldnbngb * sizeof(int) );

    num_alive = generate_initial_board( BS, &(cell(1, 1)), ldboard );

    printf("Starting number of living cells = %d\n", num_alive);
    t1 = mytimer();

    nbdone = malloc(nb_threads*sizeof(*nbdone));


    t2 = mytimer();
    temps = t2 - t1;
    printf("Final number of living cells = %d\n", num_alive);
    printf("%.2lf\n",(double)temps * 1.e3);

    free(_board);
    free(_nbngb);
    return EXIT_SUCCESS;
}

