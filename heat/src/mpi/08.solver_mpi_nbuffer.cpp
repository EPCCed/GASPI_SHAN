#include <mpi.h>
#include <algorithm>
#include <cassert>

#include "common/heat.hpp"

struct border_reqs_t {
	MPI_Request send;
	MPI_Request recv;
};

border_reqs_t *topBorderReqs = nullptr;
border_reqs_t *bottomBorderReqs = nullptr;


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

inline void startSendingFirstComputeRow(block_t *matrix, int nbx, int nby, int by, int rank)
{
	MPI_Isend(&matrix[nby+by][0], BSY, MPI_DOUBLE, rank - 1, by, MPI_COMM_WORLD, &(topBorderReqs[by].send));
}

inline void startSendingLastComputeRow(block_t *matrix, int nbx, int nby, int by, int rank)
{
	MPI_Isend(&matrix[(nbx-2)*nby + by][BSX-1], BSY, MPI_DOUBLE, rank + 1, by, MPI_COMM_WORLD, &(bottomBorderReqs[by].send));
}

inline void startReceivingUpperBorder(block_t *matrix, int nbx, int nby, int by, int rank)
{
	MPI_Irecv(&matrix[by][BSX-1], BSY, MPI_DOUBLE, rank - 1, by, MPI_COMM_WORLD, &(topBorderReqs[by].recv));
}

inline void startReceivingLowerBorder(block_t *matrix, int nbx, int nby, int by, int rank)
{
	MPI_Irecv(&matrix[(nbx-1)*nby + by][0], BSY, MPI_DOUBLE, rank + 1, by, MPI_COMM_WORLD, &(bottomBorderReqs[by].recv));
}

inline void solveBorderBlock(block_t *matrix, int nbx, int nby, int bx, int by, int rank, int rank_size, bool lastTimestep)
{
	block_t &targetBlock = matrix[bx*nby + by];
	const block_t &centerBlock = matrix[bx*nby + by];
	const block_t &topBlock    = matrix[(bx-1)*nby + by];
	const block_t &leftBlock   = matrix[bx*nby + (by-1)];
	const block_t &rightBlock  = matrix[bx*nby + (by+1)];
	const block_t &bottomBlock = matrix[(bx+1)*nby + by];
	
	const bool topCommunication = (rank > 0 && bx == 1);
	const bool bottomCommunication = (rank < rank_size-1 && bx == nbx-2);
	
	if (topCommunication) {
		MPI_Waitall(2, (MPI_Request *)&topBorderReqs[by], MPI_STATUSES_IGNORE);
	}
	
	double sum = 0.0;
	for (int x = 0; x < BSX; ++x) {
		const bool topBlockRow = (x == 0);
		const bool bottomBlockRow = (x == BSX-1);
		
		if (bottomCommunication && bottomBlockRow) {
			MPI_Waitall(2, (MPI_Request *)&bottomBorderReqs[by], MPI_STATUSES_IGNORE);
		}
		
		const row_t &topRow = (topBlockRow) ? topBlock[BSX-1] : centerBlock[x-1];
		const row_t &bottomRow = (bottomBlockRow) ? bottomBlock[0] : centerBlock[x+1];
		
		for (int y = 0; y < BSY; ++y) {
			double left = (y > 0) ? centerBlock[x][y-1] : leftBlock[x][BSY-1];
			double right = (y < BSY-1) ? centerBlock[x][y+1] : rightBlock[x][0];
			
			double value = 0.25 * (topRow[y] + bottomRow[y] + left + right);
			double diff = value - targetBlock[x][y];
			sum += diff * diff;
			targetBlock[x][y] = value;
		}
		
		if (topCommunication && topBlockRow) {
			startSendingFirstComputeRow(matrix, nbx, nby, by, rank);
			if (!lastTimestep) {
				startReceivingUpperBorder(matrix, nbx, nby, by, rank);
			}
		}
	}
	
	if (bottomCommunication) {
		startSendingLastComputeRow(matrix, nbx, nby, by, rank);
		startReceivingLowerBorder(matrix, nbx, nby, by, rank);
	}
}

inline void solveGaussSeidel(block_t *matrix, int nbx, int nby, int rank, int rank_size, bool lastTimestep)
{
	for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
			if (bx == 1 || bx == nbx-2) {
				solveBorderBlock(matrix, nbx, nby, bx, by, rank, rank_size, lastTimestep);
			} else {
				solveBlock(matrix, nbx, nby, bx, by);
			}
		}
	}
}

double solve(block_t *matrix, int rowBlocks, int colBlocks, int timesteps)
{
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);

	topBorderReqs = (border_reqs_t *) malloc(colBlocks * sizeof(border_reqs_t));
	bottomBorderReqs = (border_reqs_t *) malloc(colBlocks * sizeof(border_reqs_t));
	assert(topBorderReqs != 0);
	assert(bottomBorderReqs != 0);

	for (int by = 0; by < colBlocks; ++by) {
		topBorderReqs[by].send = MPI_REQUEST_NULL;
		topBorderReqs[by].recv = MPI_REQUEST_NULL;
		bottomBorderReqs[by].send = MPI_REQUEST_NULL;
		bottomBorderReqs[by].recv = MPI_REQUEST_NULL;
	}
	
	if (rank > 0) {
		for (int by = 1; by < colBlocks-1; ++by) {
			startReceivingUpperBorder(matrix, rowBlocks, colBlocks, by, rank);
		}
	}
	
	for (int t = 0; t < timesteps-1; ++t) {
		solveGaussSeidel(matrix, rowBlocks, colBlocks, rank, rank_size, false);
	}
	solveGaussSeidel(matrix, rowBlocks, colBlocks, rank, rank_size, true);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	free(topBorderReqs);
	free(bottomBorderReqs);
	return 0.0;
}

