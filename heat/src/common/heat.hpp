#ifndef HEAT_HPP
#define HEAT_HPP

#include <string>

#ifdef USE_SHAN
#include "SHAN_comm.h"
#endif

#include "common/matrix.hpp"

template<class T = int>
inline T round(T a, T b)
{
	return ((a + b - 1) / b) * b;
}

struct stage_t {
    int counter;
    int sendFirst;
    int recvUpper;
    int recvLower;
    int testFirst;
    int solve;
    int sendLast;
    int testLast;
};

struct HeatSource {
	float row;
	float col;
	float range;
	float temperature;
	
	HeatSource() :
		row(0.0),
		col(0.0),
		range(0.0),
		temperature(0.0)
	{
	}
};

struct HeatConfiguration {
    int timesteps;
    int rows;
    int cols;
    int rowBlocks;
    int colBlocks;
    block_t *matrix;
#ifdef USE_SHAN
    stage_t *stage;
#endif
    int numHeatSources;
    HeatSource *heatSources;
    std::string confFileName;
    std::string imageFileName;
    bool generateImage;

#ifdef USE_SHAN
    int* rank_to_idx;
    shan_neighborhood_t neighborhood_id;
    shan_segment_t dataSegment;
#endif
    HeatConfiguration() :
	timesteps(0),
	rows(0),
	cols(0),
	matrix(nullptr),
	numHeatSources(0),
	heatSources(nullptr),
	confFileName("heat.conf"),
	imageFileName("heat.ppm"),
	generateImage(false)
	{
	}
};

int initialize(HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlockOffset = 0);
int finalize(HeatConfiguration &conf);
int writeImage(std::string fileName, block_t *matrix, int rowBlocks, int colBlocks);
HeatConfiguration readConfiguration(int argc, char **argv);
void refineConfiguration(HeatConfiguration &conf, int rowValue, int colValue);
void printConfiguration(const HeatConfiguration &conf);
void initializeMatrix(const HeatConfiguration &conf, block_t *matrix, int rowBlocks, int colBlocks, int rowBlockOffset = 0);
double get_time();

#ifdef USE_SHAN
double solve(HeatConfiguration &conf, int rowBlocks, int colBlocks, int timesteps);
#else
double solve(block_t *matrix, int rowBlocks, int colBlocks, int timesteps);
#endif

#endif // HEAT_HPP
