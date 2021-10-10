#include "refine.cuh"
//extern circle refineEigen(matOut outTest); 
// PsuedoCode

//i <- 0
//foreach circle (c_j, phi_j, sigma_j) in the input
//    if (|arccos(c_j * s) - phi_j| < 3 * sigma_j) {  // if circle is good
//       A[i,1..3] <- c_j        * 1/sigma_j
//       b[i]      <- cos(phi_j) * 1/sigma_j
//       i++
//    }

// Set Boolean Array of which Circles are good 
__global__ void isGood(int *goodArr, circle *circles, circle *s, int numCircs) {
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int blockWidth = blockDim.x;
	int numBlocks = gridDim.x;
	int offset = bid * blockWidth;
	
	for (int j = offset; j < numCircs; j += blockWidth * numBlocks) {
		if (j + tid < numCircs) {
			goodArr[j + tid] = 0;
			if (abs(acosf(circles[j + tid].c * s[0].c) - circles[j + tid].phi) < (3 * circles[j + tid].dphi)) {
				goodArr[j + tid] = 1;					
			} 	
		}	
	}
}

// Create Matrix A and Vector b
__global__ void createMatrix(int *goodArr, int *tj, float *A, float *b, circle *circles, int numCircs, int *m) {
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int blockWidth = blockDim.x;
	int numBlocks = gridDim.x;
	int offset = bid * blockWidth;
	if (tid == 0) {
		m[0] = goodArr[numCircs - 1] + tj[numCircs - 1];
	}
	for (int j = offset; j < numCircs; j += blockWidth * numBlocks) {
		if (j + tid < numCircs) {
			if (goodArr[j + tid] == 1) {
				float weight = 1/circles[j + tid].dphi;
				A[tj[j + tid]] = circles[j + tid].c.x * weight;
				A[tj[j + tid] + m[0]] =  circles[j + tid].c.y * weight;
				A[tj[j + tid] + 2 * m[0]] =  circles[j + tid].c.z * weight;
				b[tj[j + tid]] = cosf(circles[j + tid].phi) * weight;
			}
		}
	}
}
