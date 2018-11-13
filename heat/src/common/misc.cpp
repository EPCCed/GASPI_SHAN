#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sys/time.h>
#include <string.h>


#include <mpi.h>

#ifdef USE_SHAN
#include "SHAN_comm.h"
#include "SHAN_type.h"
#include "common/assert.h"
#endif

#include "common/matrix.hpp"
#include "common/heat.hpp"


#ifdef USE_SHAN


void check_free(void *ptr)
{
    if(ptr != NULL)
    {
	free(ptr);
    }

    return;
}

void *check_malloc(long bytes)
{
    void *tmp;
    ASSERT(bytes > 0);
    tmp = malloc(bytes);
    ASSERT(tmp != NULL);

    return tmp;
}

int initSharedNotification(HeatConfiguration &conf, int nbx, int nby)
{
    MPI_Comm MPI_COMM_SHM;
    MPI_Comm_split_type (MPI_COMM_WORLD
			 , MPI_COMM_TYPE_SHARED
			 , 0
			 , MPI_INFO_NULL
			 , &MPI_COMM_SHM
	);


    int rank, rank_size, iProcLocal;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
    MPI_Comm_rank(MPI_COMM_SHM, &iProcLocal);

    int  const segment_id = 0;
    long const dataSz = nbx * nby * sizeof(block_t);
    int res = shan_alloc_shared(&conf.dataSegment
				, segment_id
				, SHAN_DATA
				, dataSz			  
				, MPI_COMM_SHM 
	);
    ASSERT(res == SHAN_SUCCESS);
  
    void *shm_ptr = NULL;    
    res = shan_get_shared_ptr(&conf.dataSegment
			      , iProcLocal
			      , &shm_ptr
	);
    ASSERT(res == SHAN_SUCCESS);

    conf.matrix = (block_t *) shm_ptr;

    conf.rank_to_idx = (int *) check_malloc(rank_size * sizeof(int));
    for (int i = 0; i < rank_size; ++i) 
    {
	conf.rank_to_idx[i] = -1;
    }

    conf.stage = (stage_t *) malloc(sizeof(stage_t) * (nby - 2));
    for (int i = 0; i < nby - 2; ++i) 
    {
	conf.stage[i].counter     = 0;
	conf.stage[i].sendFirst   = 0;
	conf.stage[i].recvUpper   = 0;
	conf.stage[i].recvLower   = 0;
	conf.stage[i].testFirst   = 0;
	conf.stage[i].solve       = 0;
	conf.stage[i].sendLast    = 0;
	conf.stage[i].testLast    = 0;
    }

    int num_neighbors = (rank == 0 || rank == rank_size - 1) ? 1 : 2;
    int *neighbors = (int *) check_malloc(num_neighbors * sizeof(int));

    int k = 0;
    if (rank != 0)
    {
	conf.rank_to_idx[rank - 1] = k;
	neighbors[k++] = rank - 1;
    }
    if (rank != rank_size - 1)
    {
	conf.rank_to_idx[rank + 1] = k;
	neighbors[k++] = rank + 1;
    }

    int num_type   = nby - 2;	
    long *maxSendSz = (long *) check_malloc(num_type * sizeof(long));
    long *maxRecvSz = (long *) check_malloc(num_type * sizeof(long));
    int *max_nelem_send = (int *) check_malloc(num_type * sizeof(int));
    int *max_nelem_recv = (int *) check_malloc(num_type * sizeof(int));
  
    for (int i = 0; i < num_type; ++i) 
    {
	maxSendSz[i] = BSY * sizeof(double);
	maxRecvSz[i] = BSY * sizeof(double);
	max_nelem_send[i] = 1;
	max_nelem_recv[i] = 1;
    }

    int const neighbor_hood_id = 0;
    shan_comm_init_comm(&conf.neighborhood_id
			, neighbor_hood_id
			, neighbors
			, num_neighbors
			, maxSendSz
			, maxRecvSz
			, max_nelem_send
			, max_nelem_recv
			, num_type 
			, MPI_COMM_SHM
			, MPI_COMM_WORLD
	);

    check_free(maxSendSz);
    check_free(maxRecvSz);
    check_free(max_nelem_send);
    check_free(max_nelem_recv);

    for (int i = 0; i < num_type; ++i) 
    {
	int *nelem_send;
	int *nelem_recv;
	int *send_sz;
	int *recv_sz;
	long *send_offset;
	long *recv_offset;

	shan_comm_type_offset(&conf.neighborhood_id
			      , i
			      , &nelem_send
			      , &nelem_recv
			      , &send_sz
			      , &recv_sz
			      , &send_offset
			      , &recv_offset
	    );
	
	int by = i + 1;
	for (int j = 0; j < num_neighbors; ++j) 
	{
	    if (neighbors[j] == rank - 1 && neighbors[j] != -1)
	    {		
		nelem_send[j] = 1;
		nelem_recv[j] = 1;
		send_sz[j] = BSY * sizeof(double);
		recv_sz[j] = BSY * sizeof(double);
		send_offset[j] = (nby + by) * sizeof(block_t) + 0 * sizeof(row_t); 
		recv_offset[j] =  by * sizeof(block_t) + (BSX-1) * sizeof(row_t);
	    }
	    if (neighbors[j] == rank + 1 && neighbors[j] != rank_size)
	    {		
		nelem_send[j] = 1;
		nelem_recv[j] = 1;
		send_sz[j] = BSY * sizeof(double);
		recv_sz[j] = BSY * sizeof(double);
		send_offset[j] = ((nbx-2)*nby + by) * sizeof(block_t) + (BSX-1) * sizeof(row_t);
		recv_offset[j] = ((nbx-1)*nby + by) * sizeof(block_t) + 0 * sizeof(row_t);
	    }
	}
    }

    check_free(neighbors);


}
#endif



int initialize(HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlockOffset)
{
#ifndef USE_SHAN
    conf.matrix = (block_t *) malloc(rowBlocks * colBlocks * sizeof(block_t));
    if (conf.matrix == NULL) {
	fprintf(stderr, "Error: Memory cannot be allocated!\n");
	exit(1);
    }
#else
    initSharedNotification(conf, rowBlocks, colBlocks);
#endif

    initializeMatrix(conf, conf.matrix, rowBlocks, colBlocks, rowBlockOffset);
	
    return 0;
}

int finalize(HeatConfiguration &conf)
{
    assert(conf.matrix != nullptr);

#ifndef USE_SHAN
    free(conf.matrix);
#else
    int res = shan_free_shared(&conf.dataSegment);
    ASSERT(res == SHAN_SUCCESS);

    res = shan_comm_free_comm(&conf.neighborhood_id);
    ASSERT(res == SHAN_SUCCESS);
#endif

    conf.matrix = nullptr;
    
    return 0;
}

int writeImage(std::string imageFileName, block_t *matrix, int rowBlocks, int colBlocks)
{
    // RGB table
    unsigned int r[1024], g[1024], b[1024];
	
    // Prepare the RGB table
    int n = 1023;
    for (int i = 0; i < 256; i++) {
	r[n] = 255; g[n] = i; b[n] = 0;
	n--;
    }
	
    for (int i = 0; i < 256; i++) {
	r[n] = 255 - i; g[n] = 255; b[n] = 0;
	n--;
    }
	
    for (int i = 0; i < 256; i++) {
	r[n] = 0; g[n] = 255; b[n] = i;
	n--;
    }
	
    for (int i = 0; i < 256; i++) {
	r[n] = 0; g[n] = 255 - i; b[n] = 255;
	n--;
    }
	
    // Find minimum and maximum
    double min = DBL_MAX;
    double max = -DBL_MAX;
    traverseByRows(matrix, rowBlocks, colBlocks,
		   [&](int x, int y, double value) {
		       if (value > max)
			   max = value;
		       if (value < min)
			   min = value;
		   }
	);
	
    int rows = (rowBlocks - 2) * BSX + 2;
    int cols = (colBlocks - 2) * BSY + 2;
    std::ofstream file;
    file.open(imageFileName.c_str());
	
    file << "P3" << std::endl;
    file << cols << " " << rows << std::endl;
    file << 255 << std::endl;
	
    traverseByRows(matrix, rowBlocks, colBlocks,
		   [&](int x, int y, double value) {
		       int k = 0;
		       if (max - min != 0) {
			   k = (int)(1023.0 * (value - min)/(max - min));
		       }
		       file << r[k] << " " << g[k] << " " << b[k] << "  ";
		       if (y == cols - 1) file << std::endl;
		   }
	);
	
    file.close();
	
    return 0;
}

void printUsage(int argc, char **argv)
{
	fprintf(stdout, "Usage: %s <-s size> | <-r rows -c cols> <-t timesteps> [OPTION]...\n", argv[0]);
	fprintf(stdout, "Parameters:\n");
	fprintf(stdout, "  -s, --size=SIZE\t\tuse SIZExSIZE matrix as the surface\n");
	fprintf(stdout, "  -r, --rows=ROWS\t\tuse ROWS as the number of rows of the surface\n");
	fprintf(stdout, "  -c, --cols=COLS\t\tuse COLS as the number of columns of the surface\n");
	fprintf(stdout, "  -t, --timesteps=TIMESTEPS\tuse TIMESTEPS as the number of timesteps\n\n");
	fprintf(stdout, "Optional parameters:\n");
	fprintf(stdout, "  -f, --sources-file=NAME\tget the heat sources from the NAME configuration file (default: heat.conf)\n");
	fprintf(stdout, "  -o, --output[=NAME]\t\tsave the computed matrix to a PPM file, being 'heat.ppm' the default name (disabled by default)\n");
	fprintf(stdout, "  -h, --help\t\t\tdisplay this help and exit\n\n");
}

void readParameters(int argc, char **argv, HeatConfiguration &conf)
{
	static struct option long_options[] = {
		{"size",         required_argument,  0, 's'},
		{"rows",         required_argument,  0, 'r'},
		{"cols",         required_argument,  0, 'c'},
		{"timesteps",    required_argument,  0, 't'},
		{"sources-file", required_argument,  0, 'f'},
		{"output",       optional_argument,  0, 'o'},
		{"help",         no_argument,        0, 'h'},
		{0, 0, 0, 0}
	};

	int c;
	int index;
	while ((c = getopt_long(argc, argv, "ho::fs:r:c:t:", long_options, &index)) != -1) {
		switch (c) {
			case 'h':
				printUsage(argc, argv);
				exit(0);
			case 'f':
				conf.confFileName = optarg;
				break;
			case 'o':
				conf.generateImage = true;
				if (optarg) {
					conf.imageFileName = optarg;
				}
				break;
			case 's':
				conf.rows = atoi(optarg);
				conf.cols = atoi(optarg);
				break;
			case 'r':
				conf.rows = atoi(optarg);
				break;
			case 'c':
				conf.cols = atoi(optarg);
				break;
			case 't':
				conf.timesteps = atoi(optarg);
				break;
			case '?':
				exit(1);
			default:
				abort();
		}
	}
	
	if (!conf.rows || !conf.cols || !conf.timesteps) {
		printUsage(argc, argv);
		exit(1);
	}
}

HeatConfiguration readConfiguration(int argc, char **argv)
{
	// Default configuration
	HeatConfiguration conf;
	
	// Read the execution parameters
	readParameters(argc, argv, conf);
	
	std::string line;
	std::ifstream file(conf.confFileName.c_str());
	if (!file.is_open()) {
		fprintf(stderr, "Error: Configuration file %s not found!\n", conf.confFileName.c_str());
		exit(1);
	}
	
	std::getline(file, line);
	int n = std::sscanf(line.c_str(), "%d", &(conf.numHeatSources));
	if (n != 1) {
		fprintf(stderr, "Error: Configuration file not correct!\n");
		exit(1);
	}
	
	conf.heatSources = (HeatSource *) malloc(sizeof(HeatSource) * conf.numHeatSources);
	assert(conf.heatSources != nullptr);
	
	for (int i = 0; i < conf.numHeatSources; i++) {
		std::getline(file, line);
		n = std::sscanf(line.c_str(), "%f %f %f %f",
			&(conf.heatSources[i].row),
			&(conf.heatSources[i].col),
			&(conf.heatSources[i].range),
			&(conf.heatSources[i].temperature));
		
		if (n != 4) {
			fprintf(stderr, "Error: Configuration file not correct!\n");
			exit(1);
		}
	}
	
	file.close();
	
	return conf;
}

void refineConfiguration(HeatConfiguration &conf, int rowValue, int colValue)
{
	assert(conf.rows > 0);
	assert(conf.cols > 0);
	assert(conf.timesteps > 0);
	
	if (conf.rows % rowValue) {
		// Make the number of rows divisible by the value
		fprintf(stderr, "Warning: The number of rows (%d) is not divisible by %d. Rounding it...\n", conf.rows, rowValue);
		conf.rows = round(conf.rows, rowValue);
	}
	if (conf.cols % colValue) {
		// Make the number of cols divisible by the value
		fprintf(stderr, "Warning: The number of cols (%d) is not divisible by %d. Rounding it...\n", conf.cols, colValue);
		conf.cols = round(conf.cols, colValue);
	}
}

void printConfiguration(const HeatConfiguration &conf)
{
	fprintf(stdout, "Rows x Cols       : %u x %u\n", conf.rows, conf.cols);
	fprintf(stdout, "Timesteps         : %u\n", conf.timesteps);
	fprintf(stdout, "Num. heat sources : %u\n", conf.numHeatSources);
	
	for (int i = 0; i < conf.numHeatSources; i++) {
		fprintf(stdout, "  %2d: (%2.2f, %2.2f) %2.2f %2.2f \n", i+1,
			conf.heatSources[i].row,
			conf.heatSources[i].col,
			conf.heatSources[i].range,
			conf.heatSources[i].temperature
		);
	}
}

void initializeMatrix(const HeatConfiguration &conf, block_t *matrix, int rowBlocks, int colBlocks, int rowBlockOffset)
{
	const int totalRowBlocks = conf.rowBlocks + 2;
	const int totalRows = conf.rows + 2;
	const int numRows = (rowBlocks - 2) * BSX + 2;
	const int numCols = (colBlocks - 2) * BSY + 2;
	const int rowOffset = rowBlockOffset * BSX;
	
	// Set all elements to zero
	traverseByRows(matrix, rowBlocks, colBlocks,
		[&](int x, int y, double &value) {
			value = 0;
		}
	);
	
	for (int i = 0; i < conf.numHeatSources; i++) {
		const HeatSource &src = conf.heatSources[i];
		
		// Initialize top row
		if (rowBlockOffset == 0) {
			traverseRow(matrix, rowBlocks, colBlocks, 0, 0, numCols,
				[&](int x, int y, double &value) {
					double dist = sqrt(pow((double)y / (double)numCols - src.col, 2) + pow(src.row, 2));
					if (dist <= src.range) {
						value += (src.range - dist) / src.range * src.temperature;
					}
				}
			);
		}
		
		// Initialize bottom row
		if (rowBlockOffset + rowBlocks == totalRowBlocks) {
			traverseRow(matrix, rowBlocks, colBlocks, numRows - 1, 0, numCols,
				[&](int x, int y, double &value) {
					double dist = sqrt(pow((double)y / (double)numCols - src.col, 2) + pow(1 - src.row, 2));
					if (dist <= src.range) {
						value += (src.range - dist) / src.range * src.temperature;
					}
				}
			);
		}
		
		// Initialize left column
		traverseCol(matrix, rowBlocks, colBlocks, 0, 1, numRows - 1,
			[&](int x, int y, double &value) {
				double dist = sqrt(pow(src.col, 2) + pow((double)(rowOffset + x)/(double)totalRows - src.row, 2));
				if (dist <= src.range) {
					value += (src.range - dist) / src.range * src.temperature;
				}
			}
		);
		
		// Initialize right column
		traverseCol(matrix, rowBlocks, colBlocks, numCols - 1, 1, numRows - 1,
			[&](int x, int y, double &value) {
				double dist = sqrt(pow(1 - src.col, 2) + pow((double)(rowOffset + x)/(double)totalRows - src.row, 2));
				if (dist <= src.range) {
					value += (src.range - dist) / src.range * src.temperature;
				}
			}
		);
	}
}

double get_time()
{
	struct timeval tv;
	gettimeofday(&tv, 0);
	
	return tv.tv_sec + 1e-6 * tv.tv_usec;
}

