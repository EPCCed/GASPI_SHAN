#include <mpi.h>
#include <algorithm>
#include <cassert>

#include "common/heat.hpp"


MPI_Request *topBorderReqs = nullptr;
MPI_Request *bottomBorderReqs = nullptr;

inline void solveBlock(block_t *matrix, int nbx, int nby, int bx, int by)
{
	block_t &targetBlock = matrix[bx*nby + by];
	const block_t &centerBlock = matrix[bx*nby + by];
	const block_t &topBlock    = matrix[(bx-1)*nby + by];
	const block_t &leftBlock   = matrix[bx*nby + (by-1)];
	const block_t &rightBlock  = matrix[bx*nby + (by+1)];
	const block_t &bottomBlock = matrix[(bx+1)*nby + by];
	
	double sum = 0.0;
	for (int x = 0; x < BSX; ++x) {
		const row_t &topRow = (x > 0) ? centerBlock[x-1] : topBlock[BSX-1];
		const row_t &bottomRow = (x < BSX-1) ? centerBlock[x+1] : bottomBlock[0];
		
		for (int y = 0; y < BSY; ++y) {
			double left = (y > 0) ? centerBlock[x][y-1] : leftBlock[x][BSY-1];
			double right = (y < BSY-1) ? centerBlock[x][y+1] : rightBlock[x][0];
			
			double value = 0.25 * (topRow[y] + bottomRow[y] + left + right);
			double diff = value - targetBlock[x][y];
			sum += diff * diff;
			targetBlock[x][y] = value;
		}
	}
}



inline void sendFirstComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
    for (int by = 1; by < nby-1; ++by) {
	MPI_Isend(&matrix[nby+by][0]
		  , BSY
		  , MPI_DOUBLE
		  , rank - 1
		  , by
		  , MPI_COMM_WORLD
		  , &(topBorderReqs[by])
	    );
    }
}


inline void wait4sendFirstComputeRow(int nreq)
{
    MPI_Waitall(nreq, (MPI_Request *)&topBorderReqs[1], MPI_STATUSES_IGNORE);
}

inline void sendLastComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
    for (int by = 1; by < nby-1; ++by) {
	MPI_Isend(&matrix[(nbx-2)*nby + by][BSX-1]
		  , BSY
		  , MPI_DOUBLE
		  , rank + 1
		  , by
		  , MPI_COMM_WORLD
		  , &(bottomBorderReqs[by])
	    );	
    }
}

inline void wait4sendLastComputeRow(int nreq)
{
    MPI_Waitall(nreq, (MPI_Request *)&bottomBorderReqs[1], MPI_STATUSES_IGNORE);
}

inline void receiveUpperBorder(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
    for (int by = 1; by < nby-1; ++by) {
	MPI_Recv(&matrix[by][BSX-1], BSY, MPI_DOUBLE, rank - 1, by, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

inline void receiveLowerBorder(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
    for (int by = 1; by < nby-1; ++by) {
	MPI_Recv(&matrix[(nbx-1)*nby + by][0], BSY, MPI_DOUBLE, rank + 1, by, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

inline void solveGaussSeidel(block_t *matrix, int nbx, int nby, int rank, int rank_size, int timesteps)
{
    for (int t = 0; t < timesteps; ++t) 
    {
	if (rank == 0) 
	{
	    receiveLowerBorder(matrix, nbx, nby, rank, rank_size);
	    for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
		    solveBlock(matrix, nbx, nby, bx, by);
		}
	    }
	    sendLastComputeRow(matrix, nbx, nby, rank, rank_size);
	    wait4sendLastComputeRow(nby - 2);
	} 
	else if (rank == rank_size - 1)
	{
	    sendFirstComputeRow(matrix, nbx, nby, rank, rank_size);
	    wait4sendFirstComputeRow(nby - 2);
	    receiveUpperBorder(matrix, nbx, nby, rank, rank_size);
	    for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
		    solveBlock(matrix, nbx, nby, bx, by);
		}
	    }
	}
	else
	{
	    sendFirstComputeRow(matrix, nbx, nby, rank, rank_size);
	    wait4sendFirstComputeRow(nby - 2);
	    receiveUpperBorder(matrix, nbx, nby, rank, rank_size);
	    receiveLowerBorder(matrix, nbx, nby, rank, rank_size);
	    for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
		    solveBlock(matrix, nbx, nby, bx, by);
		}
	    }
	    sendLastComputeRow(matrix, nbx, nby, rank, rank_size);
	    wait4sendLastComputeRow(nby - 2);
	}
    }
}

double solve(block_t *matrix, int rowBlocks, int colBlocks, int timesteps)
{
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	topBorderReqs = (MPI_Request *) malloc(colBlocks * sizeof(MPI_Request));
	bottomBorderReqs = (MPI_Request *) malloc(colBlocks * sizeof(MPI_Request));
	assert(topBorderReqs != 0);
	assert(bottomBorderReqs != 0);

	for (int by = 0; by < colBlocks; ++by) {
	    topBorderReqs[by] = MPI_REQUEST_NULL;
	    bottomBorderReqs[by] = MPI_REQUEST_NULL;
	}

	solveGaussSeidel(matrix, rowBlocks, colBlocks, rank, rank_size, timesteps);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	return 0.0;
}

