#pragma once
#include "datatypes.cuh"
#include "cublas_v2.h"
#include <cmath>
#include <cub/cub.cuh>

__global__ void isGood(int *goodArr, int *tj, circle *circles, circle *s, int
		numCircs);

__global__ void createMatrix(int *goodArr, int *tj, float *A, float *b, circle
		*circles, int numCircs, int *m);



