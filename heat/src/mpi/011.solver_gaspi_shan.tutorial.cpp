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


inline void TrySendFirstComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, int *output)
{
    if (input) 
    {
	int const idx = conf.rank_to_idx[rank - 1];
	ASSERT(shan_comm_notify_or_write(&conf.neighborhood_id
					     , &conf.dataSegment
					     , type_id
					     , idx
		    ) == SHAN_SUCCESS);
	(*output)++;
    }
}

inline void TryReceiveUpperBorder(HeatConfiguration &conf, int rank, int type_id, bool input, int *output)
{
    if (input) 
    {
	int const idx = conf.rank_to_idx[rank - 1];
	if (shan_comm_test4Recv(&conf.neighborhood_id
				       , &conf.dataSegment
				       , type_id
				       , idx
		 ) == SHAN_SUCCESS)
	{
	    (*output)++;
	}
    }
}

inline void TryReceiveLowerBorder(HeatConfiguration &conf, int rank, int type_id, bool input, int *output)
{
    if (input) 
    {
	int const idx = conf.rank_to_idx[rank + 1];
	if (shan_comm_test4Recv(&conf.neighborhood_id
					   , &conf.dataSegment
					   , type_id
					   , idx
		 ) == SHAN_SUCCESS)
	{
	    (*output)++;
	}
    }
}

inline void TryWaitSendFirstComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, int *output)
{
    if (input) 
    {
	int const idx = conf.rank_to_idx[rank - 1];
	if (shan_comm_test4Send(&conf.neighborhood_id
					   , type_id
					   , idx
		 ) == SHAN_SUCCESS)
	{
	    (*output)++;
	}
    }
}

inline void TrySolveBlock(HeatConfiguration &conf, int nbx, int nby, int by, bool input, int *output)
{
    if (input) 
    {
	for (int bx = 1; bx < nbx-1; ++bx) 
	{
	    solveBlock(conf, nbx, nby, bx, by);
	}
	(*output)++;
    }
}

inline void TrySendLastComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, int *output)
{
    if (input) 
    {
	int const idx = conf.rank_to_idx[rank + 1];
	ASSERT(shan_comm_notify_or_write(&conf.neighborhood_id
					     , &conf.dataSegment
					     , type_id
					     , idx
		    ) == SHAN_SUCCESS);
	(*output)++;
    }
}

inline void TryWaitSendLastComputeRow(HeatConfiguration &conf, int rank, int type_id, bool input, int *output)
{
    if (input) 
    {
	int const idx = conf.rank_to_idx[rank + 1];
	if (shan_comm_test4Send(&conf.neighborhood_id
					   , type_id
					   , idx
		     ) == SHAN_SUCCESS)
	{
	    (*output)++;
	}
    }
}



inline void solveGaussSeidel(HeatConfiguration &conf
			     , int nbx
			     , int nby
			     , int rank
			     , int rank_size
			     , int timesteps
    )
{
    stage_t *stage = conf.stage;
    bool dep = true;
    bool complete = false;

    while (!complete)
    {
	complete = true;
	for (int by = 1; by < nby - 1; ++by) 
	{
	    int type_id = by - 1;
	    if (stage[type_id].counter < timesteps)
	    {
		complete = false;
		if (rank == 0)
		{
		    dep = (stage[type_id].recvLower == stage[type_id].counter - 1);
		    TryReceiveLowerBorder(conf, rank, type_id
					  , dep
					  , &stage[type_id].recvLower
			);
		    
		    // dep = ?
		    TrySolveBlock(conf, nbx, nby, by
				  , dep
				  , &stage[type_id].solve
			);

		    // dep = ?
		    TrySendLastComputeRow(conf, rank, type_id
					  , dep
					  , &stage[type_id].sendLast
			);

		    // dep = ?
		    TryWaitSendLastComputeRow(conf, rank, type_id
					      , dep
					      , &stage[type_id].testLast
			);

		    if (stage[type_id].testLast == stage[type_id].counter)
		    {
			stage[type_id].counter++;
		    }
		} 
		else if (rank == rank_size - 1)
		{
		    dep = (stage[type_id].recvUpper == stage[type_id].counter - 1);
		    TryReceiveUpperBorder(conf, rank, type_id
					  , dep
					  , &stage[type_id].recvUpper
			);

		    // dep = ?
		    TrySendFirstComputeRow(conf, rank, type_id
					   , dep
					   , &stage[type_id].sendFirst
			);

		    // dep = ?
		    TryWaitSendFirstComputeRow(conf, rank, type_id
					       , dep
					       , &stage[type_id].testFirst
			);

		    // dep = ?
		    TrySolveBlock(conf, nbx, nby, by
				  , dep
				  , &stage[type_id].solve
			);

		    if (stage[type_id].solve == stage[type_id].counter)
		    {
			stage[type_id].counter++;
		    }
		}
		else
		{
		    // dep = ?
		    TryReceiveLowerBorder(conf, rank, type_id
					  , dep
					  , &stage[type_id].recvLower
			);


		    // dep = ?
		    TryReceiveUpperBorder(conf, rank, type_id
					  , dep
					  , &stage[type_id].recvUpper
			);


		    // dep = ?
		    TrySendFirstComputeRow(conf, rank, type_id
					   , dep
					   , &stage[type_id].sendFirst
			);


		    // dep = ?
		    TryWaitSendFirstComputeRow(conf, rank, type_id
					       , dep
					       , &stage[type_id].testFirst
			);


		    dep = (type_id == 0 || stage[type_id].solve == stage[type_id - 1].solve - 1) &&
			(stage[type_id].solve == stage[type_id].testFirst - 1) &&
			(stage[type_id].solve == stage[type_id].recvUpper - 1) &&
			(stage[type_id].solve == stage[type_id].recvLower - 1);
		    TrySolveBlock(conf, nbx, nby, by
				  , dep
				  , &stage[type_id].solve
			);

		    // dep = ?
		    TrySendLastComputeRow(conf, rank, type_id
					  , dep
					  , &stage[type_id].sendLast
			);
		    
		    // dep = ?
		    TryWaitSendLastComputeRow(conf, rank, type_id
					      , dep
					      , &stage[type_id].testLast
			);


		    if (stage[type_id].testLast == stage[type_id].counter)
		    {			
			stage[type_id].counter++;
		    }
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

    ASSERT(colBlocks - 2 > rank_size);

    stage_t *stage = conf.stage;
    for (int i = 0; i < colBlocks - 2; ++i) 
    {
	stage[i].sendFirst = -1;
	stage[i].recvUpper = -1;
	stage[i].recvLower = -1;
	stage[i].testFirst = -1;
	stage[i].solve     = -1;
	stage[i].sendLast  = -1;
	stage[i].testLast  = -1;
	stage[i].counter   =  0 ;
    }

    solveGaussSeidel(conf, rowBlocks, colBlocks, rank, rank_size, timesteps);

    MPI_Barrier(MPI_COMM_WORLD);
	
    return 0.0;
}

