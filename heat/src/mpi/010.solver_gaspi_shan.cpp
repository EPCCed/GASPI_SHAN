#include <mpi.h>
#include <algorithm>

#include "common/heat.hpp"
#include "common/assert.h"

inline void solveBlock(HeatConfiguration &conf, int nbx, int nby, int bx, int by)
{
    block_t *matrix            = conf.matrix;
    block_t &targetBlock       = matrix[bx*nby + by];
    const block_t &centerBlock = matrix[bx*nby + by];
    const block_t &topBlock    = matrix[(bx-1)*nby + by];
    const block_t &leftBlock   = matrix[bx*nby + (by-1)];
    const block_t &rightBlock  = matrix[bx*nby + (by+1)];
    const block_t &bottomBlock = matrix[(bx+1)*nby + by];
	
    double sum = 0.0;
    for (int x = 0; x < BSX; ++x) 
    {
	const row_t &topRow = (x > 0) ? centerBlock[x-1] : topBlock[BSX-1];
	const row_t &bottomRow = (x < BSX-1) ? centerBlock[x+1] : bottomBlock[0];
		
	for (int y = 0; y < BSY; ++y) 
	{
	    double left = (y > 0) ? centerBlock[x][y-1] : leftBlock[x][BSY-1];
	    double right = (y < BSY-1) ? centerBlock[x][y+1] : rightBlock[x][0];
			
	    double value = 0.25 * (topRow[y] + bottomRow[y] + left + right);
	    double diff = value - targetBlock[x][y];
	    sum += diff * diff;
	    targetBlock[x][y] = value;
	}
    }
}


inline void TrySendFirstComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, bool *output)
{
    int res;
    if (input && ! *output) 
    {
	int const idx = conf.rank_to_idx[rank - 1];
	if ((res = shan_comm_notify_or_write(&conf.neighborhood_id
					     , &conf.dataSegment
					     , type_id
					     , idx
		 )) == SHAN_SUCCESS)
	{
	    *output = true;
	}
    }
}

inline void TryReceiveUpperBorder(HeatConfiguration &conf, int rank, int type_id, bool input, bool *output)
{
    int res;
    if (input && ! *output) 
    {
	int const idx = conf.rank_to_idx[rank - 1];
	if ((res = shan_comm_test4Recv(&conf.neighborhood_id
					   , &conf.dataSegment
					   , type_id
					   , idx
		 )) == SHAN_SUCCESS)
	{
	    *output = true;
	}
    }
}

inline void TryReceiveLowerBorder(HeatConfiguration &conf, int rank, int type_id, bool input, bool *output)
{
    int res;
    if (input && ! *output) 
    {
	int const idx = conf.rank_to_idx[rank + 1];
	if ((res = shan_comm_test4Recv(&conf.neighborhood_id
					   , &conf.dataSegment
					   , type_id
					   , idx
		 )) == SHAN_SUCCESS)
	{
	    *output = true;
	}
    }
}

inline void TryWaitSendFirstComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, bool *output)
{
    int res;
    if (input && !*output) 
    {
	int const idx = conf.rank_to_idx[rank - 1];
	if ((res = shan_comm_test4Send(&conf.neighborhood_id
					   , type_id
					   , idx
		 )) == SHAN_SUCCESS)
	{
	    *output = true;
	}
    }
}

inline void TrySolveBlock(HeatConfiguration &conf, int nbx, int nby, int by, bool input, bool *output)
{
    if (input && !*output) 
    {
	for (int bx = 1; bx < nbx-1; ++bx) 
	{
	    solveBlock(conf, nbx, nby, bx, by);
	}
	*output = true;
    }
}

inline void TrySendLastComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, bool *output)
{
    int res;
    if (input && !*output) 
    {
	int const idx = conf.rank_to_idx[rank + 1];
	if ((res = shan_comm_notify_or_write(&conf.neighborhood_id
					     , &conf.dataSegment
					     , type_id
					     , idx
		 )) == SHAN_SUCCESS)
	{
	    *output = true;
	}
    }
}

inline void TryWaitSendLastComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, bool *output)
{
    int res;
    if (input && !*output) 
    {
	int const idx = conf.rank_to_idx[rank + 1];
	if ((res = shan_comm_test4Send(&conf.neighborhood_id
					   , type_id
					   , idx
		     )) == SHAN_SUCCESS)
	{
	    *output = true;
	}
    }
}



inline void solveGaussSeidel(HeatConfiguration &conf
				 , int nbx
				 , int nby
				 , int rank
				 , int rank_size
    )
{
    int res, num_blocks = 0;
    stage_t *stage = conf.stage;

    for (int i = 0; i < nby - 2; ++i) 
    {
	stage[i].sendFirst = false;
	stage[i].recvUpper = false;
	stage[i].recvLower = false;
	stage[i].testFirst = false;
	stage[i].solve     = false;
	stage[i].sendLast  = false;
	stage[i].testLast  = false;
	stage[i].complete  = false;
    }

    while (num_blocks < nby - 2)
    {
	for (int by = 1; by < nby - 1; ++by) 
	{
	    int type_id = by - 1;
	    if (!stage[type_id].complete)
	    {
		if (rank != 0)
		{
		    TrySendFirstComputeRow(conf, rank, type_id
					   , true
					   , &stage[type_id].sendFirst
			);
		    TryReceiveUpperBorder(conf, rank, type_id
					  , true
					  , &stage[type_id].recvUpper
			);
		}

		if (rank != rank_size - 1)
		{
		    TryReceiveLowerBorder(conf, rank, type_id
					  , true
					  , &stage[type_id].recvLower
			);
		}

		if (rank != 0)
		{
		    TryWaitSendFirstComputeRow(conf, rank, type_id
					       , stage[type_id].sendFirst
					       , &stage[type_id].testFirst
			);
		}
		     
		bool block_dependency;
		if (rank == 0)
		{
		    block_dependency = stage[type_id].recvLower;
		}
		else if (rank == rank_size - 1)
		{
		    block_dependency = stage[type_id].recvUpper
			&& stage[type_id].testFirst;
		}
		else
		{
		    block_dependency = stage[type_id].recvUpper
			&& stage[type_id].recvLower
			&& stage[type_id].testFirst;
		}
		
		TrySolveBlock(conf, nbx, nby, by
			      , block_dependency
			      , &stage[type_id].solve
		    );

		if (rank != rank_size - 1)
		{

		    TrySendLastComputeRow(conf, rank,type_id
					  , stage[type_id].solve
					  , &stage[type_id].sendLast
			);
		    
		    TryWaitSendLastComputeRow(conf, rank, type_id
					      , stage[type_id].sendLast
					      , &stage[type_id].testLast
			);
		}

		    
		if ((rank != rank_size - 1 && stage[type_id].testLast)
		    || (rank ==  rank_size - 1 && stage[type_id].solve))
		{		    
		    stage[type_id].complete = true;
		    num_blocks++;
		}
	    }
	}
    }


}


double solve(HeatConfiguration &conf, int rowBlocks, int colBlocks, int timesteps)
{
    int rank, rank_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
    
    MPI_Barrier(MPI_COMM_WORLD);

    for (int t = 0; t < timesteps; ++t) 
    {
	solveGaussSeidel(conf, rowBlocks, colBlocks, rank, rank_size);
    }

	
    MPI_Barrier(MPI_COMM_WORLD);
	
    return 0.0;
}

