#ifndef SOLVER_DATA_H
#define SOLVER_DATA_H

#include <mpi.h>

#ifdef USE_SHAN
#include <SHAN_comm.h>
#include <SHAN_type.h>
#endif

#define NGRAD 7
#define NFLUX 3

typedef struct 
{
  int global  __attribute__((aligned(64)));
} counter_t;

typedef struct RangeList_t
{
  // next slice - linked list
  struct RangeList_t *succ;

  // meta data
  int  start;
  int  stop;
  int  ftype; //  face type   
  
  // points of color 
  int  nfirst_points_of_color; // first touch  
  int  *first_points_of_color;
  int  nlast_points_of_color; // last touch
  int  *last_points_of_color;
  
  // comm vars, color local, send 
  int nsendcount;
  int *sendpartner;
  int *sendcount;
  int **sendindex;
  int **sendoffset;

  // comm vars, color local, recv
  int nrecvcount;
  int *recvpartner;
  int *recvcount;
  int **recvindex;
  int **recvoffset;

} RangeList;


typedef struct 
{
  int     nfaces;
  int     nallfaces;
  int     nownpoints;
  int     nallpoints;
  int     ncolors;
  int     (*fpoint)[2];
  double  (*fnormal)[3];
  double  *pvolume;
  double  (*var)[NGRAD];
  double  (*grad)[NGRAD][3];
  double  (*psd_flux)[NFLUX];
  RangeList *fcolor;
  int     niter;

  // points of color 
  int  *first_points_of_color;
  int  *last_points_of_color;

#ifdef USE_SHAN
    shan_segment_t dataSegment;
#endif
  
} solver_data ;

typedef struct 
{

  int nProc;
  int iProc;

  int ndomains;
  int ncommdomains;
  int nownpoints;
  int naddpoints;

  /* general comm */
  int *addpoint_owner;
  int *addpoint_id;
  int *commpartner;
  int *sendcount;
  int *recvcount;
  int **recvindex;
  int **sendindex;

  /* comm vars mpi */
  int nreq;
  MPI_Request *req;
  MPI_Status *stat;
  double **recvbuf;
  double **sendbuf;

  /* global stage counter */
  volatile counter_t *recv_counter;
  volatile counter_t *send_counter;

#ifdef USE_SHAN
    shan_neighborhood_t neighborhood_id;
#endif

} comm_data ;


void init_solver_data(solver_data *sd, int NITER);
void read_solver_data(int ncid, solver_data *sd);
void free_solver_data(solver_data *sd);

#endif
