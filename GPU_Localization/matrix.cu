#include <cmath>
#include "argMax.cuh"

#define threadsPerBlock 128
// A^T b
__global__ void Sgemv(int *m, int n, float *A, float *b, float *g) { 		/// Assume Column-Major Still -- This Makes Readability Really Bad -- Write out the operation to make sense of this mess
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int blockWidth = blockDim.x;
	int numBlocks = gridDim.x;

	for (int i = 0; i < 3; i += numBlocks) {				/// Each Row of the Vector will be handled by a block
		if (i + bid < 3) {
			float ind = 0;
			for (int j = 0; j < m[0]; j += blockWidth) {		/// Sum required elements by each thread
				if (j + tid < m[0]) {
					ind += A[(bid * m[0]) + (j + tid)] * b[j + tid];
				}
			}
			float gTemp = BlockReduce<float, threadsPerBlock>::sum(ind, threadsPerBlock); /// Blockwise sum to find Total
			if (tid == 0) {	
				g[bid] = gTemp; 						      /// Set Respective Vector Element
			}	
				
		}
	}
	
}
// A^T A
__global__ void Ssyrk(int *m, int n, float *A, float *H) {			/// Assume Column-Major Still -- This Makes Readability Really Bad -- Write out the operation to make sense of this mess
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int blockWidth = blockDim.x;
	int numBlocks = gridDim.x;

	for (int i = 0; i < 6; i += numBlocks) {				/// Each Element in the Matrix will be Handled by a Block -- Only Populate Lower Triangular
		if (i + bid < 6) {
			float ind = 0;
			if (i + bid < 3) {					/// First Column of H
				for (int j = 0; j < m[0]; j += blockWidth) {	/// 128 operations at a time, go until m
					if (j + tid < m[0]) {
						ind += A[((i + bid) * m[0]) + (j + tid)] * A[j + tid];		/// Sum Required Elements on each Thread
					}
				}
				float hTemp1 = BlockReduce<float, threadsPerBlock>::sum(ind, threadsPerBlock);	/// Blockwise Sum to Find Total
				if (tid == 0) {
					H[i + bid] = hTemp1;							/// Set Respective Matrix Element
				}
			}

			if (i + bid > 2 && i + bid < 5) {			/// Second Column of H
				for (int j = 0; j < m[0]; j += blockWidth) {	/// 128 operations at a time, go until m
					if (j + tid < m[0]) {
						ind += A[((i + bid - 2) * m[0]) + (j + tid)] * A[j + tid + m[0]]; 	/// Sum Required Elements on each Thread
					}
				}
				float hTemp2 = BlockReduce<float, threadsPerBlock>::sum(ind, threadsPerBlock);  /// Blockwise Sum to Find Total
				if (tid == 0) {
					H[i + bid + 1] = hTemp2;						/// Set Respective Matrix Element
				}
			}

			if (i + bid == 5) {					/// Last Column of H
				for (int j = 0; j < m[0]; j += blockWidth) {	/// 128 operations at a time, go until m
					if (j + tid < m[0]) {
						ind += A[((i + bid - 3) * m[0]) + (j + tid)] * A[((i + bid - 3) * m[0]) + (j + tid)]; /// Sum Required Elements on each Thread
					}
				}
				float hTemp3 = BlockReduce<float, threadsPerBlock>::sum(ind, threadsPerBlock); 	/// Blockwise Sum to Find Total
				if (tid == 0) {
					H[i + bid + 3] = hTemp3;								/// Set Respective Matrix Element
				}
			}

		}
	}
}
