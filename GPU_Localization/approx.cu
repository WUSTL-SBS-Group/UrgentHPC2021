#include <iostream>
// Approximation CUDA Kernel
// PsuedoCode
// circle = <c, phi, dphi>
// Given Circles c1, c2, ..., cm
// For each k  =  1 .... m				- Parallelize By Block
// 	Compute Directions s1, ..., s720 from ck 	-- Some Math from Spherical Astrophysics
// 	for each sp, 1 <= p <= 720			- Parallelize by Thread (Thread width 128 *Use for Loop*)
//		Compute LL(Sp | DataSet)		- Very long for loop (I think)
//	Save max Sp, LL					- Max of each Ck to max later
#include <cmath>
#include <random>
#include "Vector.h"
#include "datatypes.cuh"
#include "argMax.cuh"
#define radial 720
#define blockwidth 128

constexpr float PI2 = 3.141592653589793238;

// low-precision approximation to arccos.  This is substantially
// faster than the C library's acos() on the Pi.  See
// https://stackoverflow.com/questions/3380628/fast-arc-cos-algorithm
__device__ inline float fast_acos(float a)
{
  constexpr float PI = 3.14159265f;
  constexpr float C  = 0.10501094f;
  float r, s, t, u;
  t = (a < 0) ? (-a) : a;  // handle negative arguments
  u = 1.0f - t;
  s = sqrt(u + u);
  r = C * u * s + s;  // or fmaf (C * u, s, s) if FMA support in hardware
  if (a < 0) r = PI - r;  // handle negative arguments
  return r;
}
// Compute a Direction given a circle and a radial position
// cir = Circle to base Direction Off of, j = radial index
__device__ Vector<float> computeDirs(circle cir, int j) {
	Vector<float> de;
	de.x = -cir.c.y;
	de.y = cir.c.x;
	de.z = 0.00000001;
	de = (1/norm(de)) * de;
	Vector<float> dn;
	dn = cir.c & de;
	Vector<float> CV = cosf(cir.phi) * cir.c;
	float SV = sinf(cir.phi);
	float bearing = j * (2 * PI2 / radial);
	Vector<float> d = cosf(bearing) * dn + sinf(bearing) * de;
	Vector<float> dir = CV + SV * d;
	return dir;
}
// Find LL (Sp | Dataset)
// cir = DataSet Circle, s = Sp
__device__ float LL(circle cir, Vector<float> s) {
	float q = s * cir.c;
	float UB = cosf(fmaxf(cir.phi - 3 * cir.dphi, 0.0f));
	float LB = cosf(fminf(cir.phi + 3 * cir.dphi, PI2));
	float back = 0;
	if (q > LB && q < UB) {
		float aco = fast_acos(q) - cir.phi;
		back = (-0.5f * logf(cir.dphi * cir.dphi)) + -aco * aco * (1/(2 * (cir.dphi * cir.dphi)));
	}
	return back;
}

// Find best Sp, LL, for each sampled circle.
// setLengh = Number of sampled circles, eventsTot = number of total circles
// circles = Sampled Circles, allCircs = All Circles
// out = Array of best sp, LL for each sampled circle 

__global__ void approx(int setLength,int eventsTot, int *circles, circle *allCircs, approxOut *out, int *bigSamp) {
	int bid = blockIdx.x;
	int tid = threadIdx.x;	
	int numBlocks = gridDim.x;
	auto &bcast = broadcast<int, blockwidth>;
								
	for (int k = 0; k < setLength; k += numBlocks) { 								// If setLength > numBlocks
		if (k + bid < setLength) {	
			out[k + bid].LL = -INFINITY;
			float ll = -INFINITY;										// Best Per Thread LL
			Vector<float> sp;										// Best Per Thread Sp
			for (int i = 0; i < (radial / blockwidth); i++) { 						//128 Threads at a time for all sp
				float llTemp = 0;
				Vector<float> spTemp = computeDirs(allCircs[circles[k + bid]], i * blockwidth + tid); 	// Compute a direction 
				for (int j = 0; j < eventsTot; j++) {							// find LL for a direction
					llTemp += LL(allCircs[bigSamp[j]], spTemp);
				}
				if (llTemp > ll) {
					ll = llTemp;
					sp = spTemp;
				}
			}
			int maxTid = BlockArgMax<int, float, blockwidth>::argmax(tid,ll);		// Find Max on block
			maxTid = bcast(maxTid, 0);							// Broadcast Index
	
			if (tid == maxTid) {								// If thread containing max 
				out[k + bid].LL = ll;
				out[k + bid].c = sp; 				// You now have max LL and best sp for each sampled Circle
			}
		}
	}
	
}

// Convert from Eta/dEtaSq to phi/sigma
__global__ void phiSigma(int eventsTot, circle *circles) {
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	const int blockWidth = blockDim.x;
	const int numBlocks = gridDim.x;
	int offset = bid * blockWidth;
	for (int j = offset; j < eventsTot; j += numBlocks * blockWidth) {
		if (j + tid < eventsTot) {
			float eta = circles[j + tid].phi;
			float detaSq = circles[j + tid].dphi;
			float phi = fast_acos(eta);
			float sigma = sqrtf(detaSq/(1 - eta * eta));
			circles[j + tid].phi = phi;
			circles[j + tid].dphi = sigma;
		}
	}
}
