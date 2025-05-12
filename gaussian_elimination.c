/*
by: Ira Garrett
March 13 2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define  N  4

int      thread_index[N];
double   A[N][N+1];
int pivotComplted[N] = {0};

pthread_barrier_t    barrier;


/********************************************************
** function: print_matrix() - print the matrix A       **
********************************************************/
int  print_matrix()
{
    int i, j;

    printf("--------------------------------------\n");

    for ( i = 0; i < N; i++ ) {
        printf( "|" );

        for ( j = 0; j < N + 1; j++ )
            printf( "%6.2f ", A[i][j] );

        printf( " |\n" );
    }
    printf("--------------------------------------\n");
}


/**************************************************
** threads function: ge() - Gaussian elimination **
**************************************************/
void  *ge( void *arg )
{
    int      i, j, k, prow;
    int      myid = *( (int *) arg );
    double   temp,  factor;

    for ( i = 0; i < N - 1; i++ ) {
        if ( i == myid ) {
            printf( "\nPartial pivoting executed by thread %d on row %d:\n", myid, i + 1 );
            temp = 0.0;
            prow = i;

            for ( j = i; j <= N; j++ ) {
                if ( fabs( A[j][i] ) > temp ) {
                    temp = fabs( A[j][i] );
                    prow = j;
                }                                // end IF branch 'fabs > temp'
            }                                  // end FOR loop 'j'

            printf( "   Pivot_row = %d,  pivot value = %6.2f\n", prow + 1, A[prow][i] );

            if ( prow != i ) {                 // test for pivot

                //  swap rows

                for ( j = i; j < N + 1; j++ ) {
                    temp = A[i][j];
                    A[i][j] = A[prow][j];
                    A[prow][j] = temp;
                }                                // end FOR loop 'j'
            }                                  // end IF branch 'prow != i'
        }                                    // end IF branch 'i == myid'

        // wait for partial pivoting done

        printf( "Thread %d stops and waits for the other threads to perform pivots\n", myid );

        pthread_barrier_wait( &barrier );

        for ( j = i + 1; j < N; j++ ) {
            if ( j == myid ) {
                printf( "Thread %d performs factorization on row %d\n", myid, j + 1 );
                factor = A[j][i] / A[i][i];

                for ( k = i + 1; k <= N; k++ )
                    A[j][k] -= A[i][k] * factor;

                A[j][i] = 0.0;
            }                                  // end IF branch 'j == myid'
        }                                    // end FOR loop 'j'

        // wait for current row reductions to finish

        printf( "Thread %d stops and waits for the other threads to factor their rows\n", myid );

        pthread_barrier_wait( &barrier );

        if ( i == myid )
            print_matrix();
    }                                      // end FOR loop 'i'
}                                        // end function ge


/****************************************************************
** function: main - initialize matrix for Gaussian elimination **
****************************************************************/
int  main( int argc, char *argv[] )
{

    int     i, j, NTHREADS;
    double  sum;

    //////////////////////////
    /*my code additions start*////         //all the references to "N" are replaced with NTHREADS////
    /////////////////////////
    printf("Number of chars :\"argc\" = %d\n",argc);


    /*
     * logic to determine the number of threads to create
     */

    //exit if there is a problem
    if(argc <0 ) exit(-1);

    //set NTHREADS to be N if there are no arguments or more than 1
    if(argc == 1 || argc > 2)
    {
        printf("# NTHREADS = %d because you inputted %d parameters instead of 1\n", (N), argc-1);
        NTHREADS = N;
    }

    if(argc == 2)
    {
        int numberOfThreads = atoi(argv[1]);
        if(!(numberOfThreads > N)) {
            printf("#NTHREADS = %d\n", numberOfThreads);
            NTHREADS = numberOfThreads;
        }

        else
        {
            printf("the number you inputted \"%d\" is larger than the allowed size of %d.\n#NTHREADS  = %d\n",atoi(argv[1]),N,N);
            NTHREADS = N;
        }
    }



    pthread_t   threads[NTHREADS];

    printf( "main: initialize matrix A[N][N+1] as [A|B]\n" );

    for ( i = 0; i < N; i++ )
        for ( j = 0; j < N; j++ )
            A[i][j] = 1.0;

    for ( i = 0; i < N; i++ )
        A[i][N-i-1] = 1.0 * N;

    for ( i = 0; i < N; i++ )
        A[i][N] = 2.0 * N - 1;

    print_matrix();                                 // show initial matrix [A|B]

    pthread_barrier_init( &barrier, NULL, NTHREADS );      // set up barrier

    printf( "main: create N = %d working threads\n", N );

    //allocate and re-allocate created threads to go through each row and solve
    for(int j = 0; (pivotComplted[j] == 0); j++) {
        printf("PivotCompleted [%d] = %d\n",j, pivotComplted[i]);

        for (int i = 0; i < NTHREADS; i++) {
            //if (pivotComplted[i] == 0) {
            thread_index[i] = i;
            pthread_create(&threads[i], NULL, ge, (void *) &thread_index[i]);
            pivotComplted[i] = 1;
            //}
        }
        if(i == NTHREADS) i = 0;
    }


    printf( "main: wait for all %d working threads to join\n", N );

    for ( i = 0; i < NTHREADS; i++ ) {
        pthread_join( threads[i], NULL );
    }
    printf( "main: back substitution : " );

    for ( i = NTHREADS - 1; i >= 0; i-- ) {
        sum = 0.0;

        for ( j = i + 1; j < NTHREADS; j++ )
            sum += A[i][j] * A[j][NTHREADS];

        A[i][NTHREADS] = (A[i][NTHREADS]- sum) / A[i][i];
    }
    // print solution

    printf( "The solution is :\n" );

    //iterate through all the elements in the matrix (N)
    for ( i = 0; i < N; i++ ) {
        printf( "%6.2f ", A[i][NTHREADS] );
    }
    printf("\n");
}
