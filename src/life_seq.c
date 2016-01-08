#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

int BS;

#define cell( _i_, _j_ ) board[ ldboard * (_j_) + (_i_) ]
#define ngb( _i_, _j_ )  nbngb[ ldnbngb * ((_j_) - 1) + ((_i_) - 1 ) ]

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

int main(int argc, char* argv[])
{
    int i, j, loop, num_alive, maxloop;
    int ldboard, ldnbngb;
    double t1, t2;
    double temps;

    int *board;
    int *nbngb;

    if (argc < 3) {
	printf("Usage: %s nb_iterations size\n", argv[0]);
	return EXIT_SUCCESS;
    } else {
	maxloop = atoi(argv[1]);
	BS = atoi(argv[2]);
	//printf("Running sequential version, grid of size %d, %d iterations\n", BS, maxloop);
    }
    num_alive = 0;

    /* Leading dimension of the board array */
    ldboard = BS + 2;
    /* Leading dimension of the neigbour counters array */
    ldnbngb = BS;

    board = malloc( ldboard * ldboard * sizeof(int) );
    nbngb = malloc( ldnbngb * ldnbngb * sizeof(int) );

    num_alive = generate_initial_board( BS, &(cell(1, 1)), ldboard );

    printf("Starting number of living cells = %d\n", num_alive);
    t1 = mytimer();

    for (loop = 1; loop <= maxloop; loop++) {

	cell(   0, 0   ) = cell(BS, BS);
	cell(   0, BS+1) = cell(BS,  1);
	cell(BS+1, 0   ) = cell( 1, BS);
	cell(BS+1, BS+1) = cell( 1,  1);

	for (i = 1; i <= BS; i++) {
	    cell(   i,    0) = cell( i, BS);
	    cell(   i, BS+1) = cell( i,  1);
	    cell(   0,    i) = cell(BS,  i);
	    cell(BS+1,    i) = cell( 1,  i);
	}


	for (j = 1; j <= BS; j++) {
	    for (i = 1; i <= BS; i++) {
		ngb( i, j ) =
		    cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
		    cell( i-1, j   ) +                  cell( i+1, j   ) +
		    cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
	    }
	}

	num_alive = 0;
	for (j = 1; j <= BS; j++) {
	    for (i = 1; i <= BS; i++) {
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

        /* Avec les celluls sur les bords (utile pour vérifier les comm MPI) */
        /* output_board( BS+2, &(cell(0, 0)), ldboard, loop ); */

        /* Avec juste les "vraies" cellules: on commence à l'élément (1,1) */
        output_board( BS, &(cell(1, 1)), ldboard, loop);

	printf("%d cells are alive\n", num_alive);
    }

    t2 = mytimer();
    temps = t2 - t1;
    printf("Final number of living cells = %d\n", num_alive);
    printf("%.2lf\n",(double)temps * 1.e3);

    free(board);
    free(nbngb);
    return EXIT_SUCCESS;
}

