#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

int BS;

#define cell( _i_, _j_ ) board[ ldboard * (_j_) + (_i_) ]
#define ngb( _i_, _j_ )  nbngb[ ldnbngb * ((_j_) - 1) + ((_i_) - 1 ) ]

enum direction {LEFT = 0, UP = 1, RIGHT = 2, DOWN = 3,UPLEFT = 4, UPRIGHT = 5, DOWNLEFT = 6, DOWNRIGHT = 7};


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

void make_neighbours_table(int * neighbours, MPI_Comm comm_cart) {
  int displ = 1;
  int index = 1;
  int coords[2];
  int myrank;

  MPI_Comm_rank(comm_cart, &myrank);

  MPI_Cart_shift(comm_cart, index, displ, &neighbours[LEFT], &neighbours[RIGHT]);
  index = 0;
  MPI_Cart_shift(comm_cart, index, displ, &neighbours[UP], &neighbours[DOWN]);  
  MPI_Cart_coords(comm_cart, myrank, 2, coords);
  coords[0]--;
  coords[1]--;
  MPI_Cart_rank(comm_cart, coords, &neighbours[UPLEFT]);
  coords[0]+=2;
  MPI_Cart_rank(comm_cart, coords, &neighbours[UPRIGHT]);
  coords[1]+=2;
  MPI_Cart_rank(comm_cart, coords, &neighbours[DOWNRIGHT]);
  coords[0]-=2;
  MPI_Cart_rank(comm_cart, coords, &neighbours[DOWNLEFT]);
}

void make_communications(MPI_Request * req, MPI_Comm comm_cart, int * neighbours, int block_size, int * board, int ldboard, MPI_Datatype block_line) {
  MPI_Send_init(&cell(1, 1),block_size, MPI_INT, neighbours[LEFT], 0, comm_cart, &req[0]);
  MPI_Recv_init(&cell(1, block_size + 1),block_size, MPI_INT, neighbours[RIGHT], 0, comm_cart, &req[1]); //GAUCHE

  MPI_Send_init(&cell(1, 1), 1, block_line, neighbours[UP], 0, comm_cart, &req[2]);
  MPI_Recv_init(&cell(block_size + 1, 1), 1, block_line, neighbours[DOWN], 0, comm_cart, &req[3]); //HAUT

  MPI_Send_init(&cell(1, block_size),block_size, MPI_INT, neighbours[RIGHT], 0, comm_cart, &req[4]);
  MPI_Recv_init(&cell(1, 0),block_size, MPI_INT, neighbours[LEFT], 0, comm_cart, &req[5]); //DROITE

  MPI_Send_init(&cell(block_size, 1), 1, block_line, neighbours[DOWN], 0, comm_cart, &req[6]);
  MPI_Recv_init(&cell(0, 1), 1, block_line, neighbours[UP], 0, comm_cart, &req[7]); //BAS

  MPI_Send_init(&cell(1, 1), 1, MPI_INT, neighbours[UPLEFT], 0, comm_cart, &req[8]);
  MPI_Recv_init(&cell(block_size + 1, block_size + 1), 1, MPI_INT, neighbours[DOWNRIGHT], 0, comm_cart, &req[9]);  //HAUT-GAUCHE
  
  MPI_Send_init(&cell(1, block_size), 1, MPI_INT, neighbours[UPRIGHT], 0, comm_cart, &req[10]);
  MPI_Recv_init(&cell(block_size + 1, 0), 1, MPI_INT, neighbours[DOWNLEFT], 0, comm_cart, &req[11]);  //HAUT-DROITE

  MPI_Send_init(&cell(block_size, 1), 1, MPI_INT, neighbours[DOWNLEFT], 0, comm_cart, &req[12]);
  MPI_Recv_init(&cell(0, block_size + 1), 1, MPI_INT, neighbours[UPRIGHT], 0, comm_cart, &req[13]);  //BAS-GAUCHE

  MPI_Send_init(&cell(block_size, block_size), 1, MPI_INT, neighbours[DOWNRIGHT], 0, comm_cart, &req[14]);
  MPI_Recv_init(&cell(0, 0), 1, MPI_INT, neighbours[UPLEFT], 0, comm_cart, &req[15]);  //BAS-DROIT
  

}

int main(int argc, char* argv[])
{
    int i, j, loop, num_alive, maxloop;
    int ldgboard,ldboard, ldnbngb;
    double t1, t2;
    double temps;
    int *gboard;
    int *board;
    int *nbngb;

    int size;
    int coord[2], id;
    int procs_per_lines_col;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    procs_per_lines_col = sqrt(size);
    if(procs_per_lines_col * procs_per_lines_col != size) {
      fprintf(stderr, "Renseignez un nombre carré de processeurs siouplait !\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    int dims[2]; dims[0] = procs_per_lines_col; dims[1] = procs_per_lines_col;
    int periods[2]; periods[0] = 1; periods[1] = 1;
    MPI_Comm comm_cart;
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);
    MPI_Comm_rank(comm_cart, &id);
    MPI_Cart_coords(comm_cart, id, 2, coord);

    if (argc < 3) {
	printf("Usage: %s nb_iterations size\n", argv[0]);
	return EXIT_SUCCESS;
    } else {
	maxloop = atoi(argv[1]);
	BS = atoi(argv[2]);
	//printf("Running sequential version, grid of size %d, %d iterations\n", BS, maxloop);
    }
    num_alive = 0;


    //Generate the neighbours table
    
    /* Leading dimension of the global board array */
    ldgboard = BS + 2;
    /* Leading dimension of the board array */
    ldboard = BS/procs_per_lines_col + 2;
    /* Leading dimension of the neigbour counters array */
    ldnbngb = BS/procs_per_lines_col;

    board = malloc( ldboard * ldboard * sizeof(int) );
    nbngb = malloc( ldnbngb * ldnbngb * sizeof(int) );
    
    if(id == 0) {
      gboard = malloc(ldgboard * ldgboard * sizeof(int));
      num_alive = generate_initial_board( BS, &gboard[1+ldgboard], ldgboard );
      //fprintf(stderr,"Starting number of living cells = %d\n", num_alive);
    }

    MPI_Datatype block;
    MPI_Type_vector(ldboard-2, ldboard-2, ldgboard, MPI_INT, &block);
    MPI_Type_create_resized(block, 0, sizeof(int), &block);
    MPI_Type_commit(&block);

    MPI_Datatype subblock;
    MPI_Type_vector(ldboard-2, ldboard-2, ldboard, MPI_INT, &subblock);
    MPI_Type_create_resized(subblock, 0, sizeof(int), &subblock);
    MPI_Type_commit(&subblock);
    
    int * counts = (int*) malloc(size*sizeof(int));
    int * displs = (int*) malloc(size*sizeof(int));
    // Définition des déplacements pour chaque proc
    for (int i = 0; i < procs_per_lines_col; ++i)
      {
	for (int j = 0; j < procs_per_lines_col; ++j)
	  {
	    counts[i+j*procs_per_lines_col]= 1;
	    displs[i+j*procs_per_lines_col]= i*ldgboard*(ldboard-2)+j*(ldboard-2);
	  }
      }
    MPI_Scatterv(&gboard[1+ldgboard], counts, displs, block, &board[ldboard+1], 1,
				subblock,0, comm_cart);
    

    int neighbours[8];
    make_neighbours_table(neighbours, comm_cart);    
    MPI_Request req[16];

    int block_size = ldboard - 2;
    MPI_Datatype block_line;
    MPI_Type_vector(block_size+2, 1, ldboard,MPI_INT, &block_line);
    MPI_Type_commit(&block_line);

    t1 = mytimer();
    make_communications(req, comm_cart, neighbours, block_size, board, ldboard, block_line);
    for (loop = 1; loop <= maxloop; loop++) {
      MPI_Startall(16, req);
	  
	  /*	cell(   0, 0   ) = cell(BS, BS);
	cell(   0, BS+1) = cell(BS,  1);
	cell(BS+1, 0   ) = cell( 1, BS);
	cell(BS+1, BS+1) = cell( 1,  1);

	for (i = 1; i <= BS; i++) {
	    cell(   i,    0) = cell( i, BS);
	    cell(   i, BS+1) = cell( i,  1);
	    cell(   0,    i) = cell(BS,  i);
	    cell(BS+1,    i) = cell( 1,  i);
	}
	  */

      //Inner cells 
	for (j = 2; j <= block_size; j++) {
	    for (i = 2; i <= block_size; i++) {
		ngb( i, j ) =
		    cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
		    cell( i-1, j   ) +                  cell( i+1, j   ) +
		    cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
	    }
	}

	//On LEFT
	MPI_Waitall(2, &req[0], MPI_STATUS_IGNORE);
	MPI_Waitall(2, &req[8], MPI_STATUS_IGNORE);
	MPI_Waitall(2, &req[12], MPI_STATUS_IGNORE);
	//CALCUL LIGNE GAUCHE
	for(j = 1; j <= block_size; j++) {
	  ngb( 1, j ) =
	    cell( 0, j-1 ) + cell( 1, j-1 ) + cell( 2, j-1 ) +
	    cell( 0, j   ) +                  cell( 2, j   ) +
	    cell( 0, j+1 ) + cell( 1, j+1 ) + cell( 2, j+1 );
	}
	
	//On TOP
	MPI_Waitall(2, &req[2], MPI_STATUS_IGNORE);
	MPI_Waitall(2, &req[10], MPI_STATUS_IGNORE);
	//CALCUL LIGNE DESSUS
	for(i = 1; i <= block_size; i++) {
	  ngb( i, 1 ) =
	    cell( i - 1, 0) + cell( i, 0 ) + cell( i + 1, 0 ) +
	    cell( i - 1, 1) +                cell( i + 1, 1 ) +
	    cell( i - 1, 2) + cell( i, 2 ) + cell( i + 1, 2 );
	}


	//On RIGHT
	MPI_Waitall(2, &req[4], MPI_STATUS_IGNORE);
	MPI_Waitall(2, &req[14], MPI_STATUS_IGNORE);
	//CALCULER A DROITE
	for(j = 1; j <= block_size; j++) {
	  ngb( block_size, j ) =
	    cell( block_size - 1, j-1 ) + cell( block_size , j-1 ) + cell( block_size + 1, j-1 ) +
	    cell( block_size - 1, j   ) +                            cell( block_size + 1, j   ) +
	    cell( block_size - 1, j+1 ) + cell( block_size, j+1 ) + cell(  block_size + 1, j+1 );
	}
	

	
	//ON BOT
	MPI_Waitall(2,&req[6], MPI_STATUS_IGNORE);
	//CALCULER EN BAS
	for(i = 1; i <= block_size; i++) {
	  ngb( i, block_size ) =
	    cell( i - 1, block_size - 1) + cell( i, block_size - 1 ) + cell( i + 1, block_size - 1 ) +
	    cell( i - 1, block_size ) +                cell( i + 1, block_size ) +
	    cell( i - 1, block_size + 1 ) + cell( i, block_size + 1 ) + cell( i + 1, block_size + 1 );
	}

	
	// Met a jour la partie
	num_alive = 0;
	for (j = 1; j <= block_size; j++) {
	    for (i = 1; i <= block_size; i++) {
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
	//output_board( BS, &(cell(1, 1)), ldboard, loop);

	//printf("%d cells are alive\n", num_alive);
    }
    MPI_Gatherv(&board[ldboard+1], 1, subblock,&gboard[ldgboard+1], counts,displs, block, 0, comm_cart);

    t2 = mytimer();

    temps = t2 - t1;
    MPI_Allreduce(MPI_IN_PLACE,&temps, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
    MPI_Allreduce(MPI_IN_PLACE,&num_alive, 1, MPI_INT, MPI_SUM, comm_cart);
    if(id == 0) {
      //printf("Final number of living cells = %d\n", num_alive);
      printf("%.2lf\n",(double)temps * 1.e3);
    }
    free(board);
    free(nbngb);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

