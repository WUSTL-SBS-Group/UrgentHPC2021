#include "datatypes.cuh"
#include<cmath>
#include<math.h>
#include<iostream>
#include<fstream>
#include<sys/time.h>
#include<random>
#include <cub/cub.cuh>

#include "matrix.cu"
#include "util.h"

#define samples 20
#define samplesAll 1000
#define iters 20
#define threadsPerBlock 128

using namespace std::chrono;

std::random_device dev;
std::mt19937 prng(dev());

extern __global__ void approx(int setLength, int eventsTot, int *circles, circle *allCircs, approxOut *out, int *bigSamp);
extern __global__ void phiSigma(int eventsTot, circle *circles);

extern __global__ void isGood(int *goodArr, circle *circles, circle *s, int numCircs);
extern __global__ void createMatrix(int *goodArr, int *tj, float *A, float *b, circle *circles, int numCircs, int *m); 
extern circle refineEigen(matOut outTest);

using namespace std;
/// Parse Circles Method ///
int parseCircles(std::string filename, circle *Circles) {
	ifstream inFile(filename);
	string txt;
	int i = 0;
	while(!inFile.eof()) {
		getline(inFile,txt);
		if (inFile.eof())
			break;
		circle tempCirc;
		getline(inFile,txt);
		tempCirc.c.x = stof(txt);
		getline(inFile,txt);
		tempCirc.c.y = stof(txt);
		getline(inFile,txt);
		tempCirc.c.z = stof(txt);
		getline(inFile,txt);
		tempCirc.phi = stof(txt);
		getline(inFile,txt);
		tempCirc.dphi = stof(txt);
		if (tempCirc.phi < 1 && tempCirc.phi > -1) {
			Circles[i] = tempCirc;
			i++;
		}
	}
	return i;	
}
/// Method for randomly sampling indices ///
void sample(int *bigSample, int *smallSample, int n, int events, int smallSamples) {
	int *sampleSet = new int[events];
	for (int i = 0; i < events; i++) {
		sampleSet[i] = i;
	}
	for (int j = 0; j < n; j++) {							/// Sample large set
		uniform_int_distribution<int> D(j, events - 1);
		bigSample[j] = j;
		swap(bigSample[j], sampleSet[D(prng)]);
	}
	for (int j = 0; j < smallSamples; j++) {					/// Sample small set from large set
		smallSample[j] = bigSample[j];
	}
}

/// Method for comparison with Marion's Code /// 
void fixedCompare(int *bigSample, int *smallSample, int n, int events, int smallSamples) {
	for (int j = 0; j < n; j++){
		bigSample[j] = j;
	}
	for (int j = 0; j < smallSamples; j++) {
		smallSample[j] = j;
	}
}

int main(int argc, char** argv) {	
	/// Find Appropriate Number of Blocks per Kernel ///
	int nBlocksPer_approx, nBlocksPer_phiSig, nBlocksPer_isGood, nBlocksPer_createMatrix;
	int nBlocksPer_Sgemv, nBlocksPer_Ssyrk;
	cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nBlocksPer_approx, approx, threadsPerBlock, 0);
	cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nBlocksPer_phiSig, phiSigma, threadsPerBlock, 0);
	cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nBlocksPer_isGood, isGood, threadsPerBlock, 0);
	cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nBlocksPer_createMatrix, createMatrix, threadsPerBlock, 0); 
	cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nBlocksPer_Sgemv, Sgemv, threadsPerBlock, 0);
	cudaOccupancyMaxActiveBlocksPerMultiprocessor(&nBlocksPer_Ssyrk, Ssyrk, threadsPerBlock, 0);
	int device, numSMs;
	cudaGetDevice(&device);
	cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, device);
	int nBlocks_approx = numSMs * nBlocksPer_approx;
	int nBlocks_phiSig = numSMs * nBlocksPer_phiSig;
	int nBlocks_isGood = numSMs * nBlocksPer_isGood;
	int nBlocks_createMatrix = numSMs * nBlocksPer_createMatrix;
	int nBlocks_Sgemv = numSMs * nBlocksPer_Sgemv;
	int nBlocks_Ssyrk = numSMs * nBlocksPer_Ssyrk;
	
	/// Set up Timing ///
	struct timeval start_time;
	struct timeval approx_time;
	struct timeval end_time;
	struct timeval refineKernels_time;
	struct timeval cublas_time;
	struct timeval eigen_time;
	float refineSec = 0; float refineuSec = 0;
	float cublasSec = 0; float cublasuSec = 0;
	float eigenSec = 0; float eigenuSec = 0;

	/// Parse Data/Set up Circles ///
	circle *allCircs = new circle[100000];
	int eventsActual = parseCircles(argv[1], allCircs);	
	circle s;
	approxOut *bestCircs = new approxOut[samples];
	int numEvents = min(samplesAll, eventsActual);
	
	float *g = new float[3];
	float *H = new float[9];

	/// Find Random Indices ///
	int *sampled = new int[samples]; 		/// Array of sampled indices for sampled circles in approximation
	int *bigSample = new int[numEvents]; 		/// Array of sampled indices for all circles in approximation
	sample(bigSample, sampled, numEvents, eventsActual, samples);		/// Use for Random Sampling

//	fixedCompare(bigSample, sampled, samplesAll, eventsActual, samples);		/// Use when trying to get deterministic behavior

	/// Set up GPU Memory ///
	circle *allCircsGPU, *sGPU;
	int *sampledGPU, *goodArrGPU, *tjGPU, *m, *bigSampGPU;
	approxOut *outGPU;
	float *gpuA, *gpuB, *gpuG, *gpuH;
	void *d_temp_storage = NULL;
	size_t temp_storage_bytes = 767;	/// THIS MAY NOT BE SAFE

	cudaMalloc(&m, sizeof(int));
	cudaMalloc(&d_temp_storage, temp_storage_bytes);
	cudaMalloc(&sampledGPU, samples * sizeof(int));
	cudaMalloc(&outGPU, samples * sizeof(approxOut));
	cudaMalloc(&goodArrGPU, eventsActual * sizeof(int));
	cudaMalloc(&tjGPU, eventsActual * sizeof(int));
	cudaMalloc(&sGPU, sizeof(circle));
	cudaMalloc(&gpuG, 3 * sizeof(float));
	cudaMalloc(&gpuH, 9 * sizeof(float));
	cudaMalloc(&gpuA, 3 * eventsActual * sizeof(float));
	cudaMalloc(&gpuB, eventsActual * sizeof(float));
	cudaMalloc(&allCircsGPU, eventsActual * sizeof(circle));
	cudaMalloc(&bigSampGPU, numEvents * sizeof(int));

	/// Start Calculation /// 

	gettimeofday(&start_time,  NULL);			/// Start Timer

	cudaMemcpy(allCircsGPU, allCircs, eventsActual * sizeof(circle), cudaMemcpyHostToDevice);
	cudaMemcpy(sampledGPU, sampled, samples * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(bigSampGPU, bigSample, numEvents * sizeof(int), cudaMemcpyHostToDevice);

	/// Convert Circles to Proper Format ///
	phiSigma<<<nBlocks_phiSig, threadsPerBlock>>>(eventsActual, allCircsGPU);

	/// Approximate Best Circles ///
	approx<<<nBlocks_approx, threadsPerBlock>>>(samples, numEvents, sampledGPU, allCircsGPU, outGPU, bigSampGPU);
	cudaMemcpy(bestCircs, outGPU, samples * sizeof(approxOut), cudaMemcpyDeviceToHost);
	
	// Find Weighted Average Best s ///
	float Lmax = -100000;
	for (int i = 0; i < samples; i++) {
		Lmax = std::max(Lmax, bestCircs[i].LL);
	}
	s.c.x = 0;
	s.c.y = 0;
	s.c.z = 0;
	for (int i = 0; i < samples; i++) {
		s.c = s.c + (expf(bestCircs[i].LL - Lmax) * bestCircs[i].c);
	}	
	normalize(s.c);	
	
	gettimeofday(&approx_time, NULL);			/// Set Timer for Approximation End
	
	/// Refine ///
	circle bestGuess;
	bestGuess = s;
	for (int i = 0; i < iters; i++) {
		/// Find Good Circles/Indices to Write to ///
		gettimeofday(&refineKernels_time, NULL);		/// Set Timer for Refinement Kernels

		cudaMemcpy(sGPU, &bestGuess, sizeof(circle), cudaMemcpyHostToDevice);

		isGood<<<nBlocks_isGood, threadsPerBlock>>>(goodArrGPU, allCircsGPU, sGPU, eventsActual);
		
		cub::DeviceScan::ExclusiveSum(d_temp_storage,temp_storage_bytes,goodArrGPU,tjGPU, eventsActual);	

		/// Create A and B ///
		createMatrix<<<nBlocks_createMatrix, threadsPerBlock>>>(goodArrGPU,tjGPU,gpuA,gpuB, allCircsGPU, eventsActual, m); // m is Now Calculated in this Kernel
		cudaDeviceSynchronize();

		/// Linear Algebra Operations ///
		gettimeofday(&cublas_time, NULL);			/// Set Timer for Linear Algebra Operations

		Sgemv<<<nBlocks_Sgemv, threadsPerBlock>>>(m, 3, gpuA, gpuB, gpuG); // A^T b
		Ssyrk<<<nBlocks_Ssyrk, threadsPerBlock>>>(m, 3, gpuA, gpuH); // A^T A
	
		cudaMemcpy(g, gpuG, 3 * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(H, gpuH, 9 * sizeof(float), cudaMemcpyDeviceToHost);

		cudaDeviceSynchronize();

		/// Refine Eigen ///
		gettimeofday(&eigen_time, NULL);			/// Set Timer for Eigen Operations
		
		bestGuess = refineEigen(matOut(g,H));

		gettimeofday(&end_time, NULL);			/// End Timer

		/// Loop Timing ///
		refineSec += cublas_time.tv_sec - refineKernels_time.tv_sec;
		refineuSec += cublas_time.tv_usec - refineKernels_time.tv_usec;
		cublasSec += eigen_time.tv_sec - cublas_time.tv_sec;
		cublasuSec += eigen_time.tv_usec - cublas_time.tv_usec;
		eigenSec += end_time.tv_sec - eigen_time.tv_sec;
		eigenuSec += end_time.tv_usec - eigen_time.tv_usec;

	}
	
	/// PRINTS / RESULTS ///

	/// Num Events ///
	
	cout << eventsActual << " Events\n";

	/// Accuracy ///

	Vector<float> actual;
	actual.x = 0;
	actual.y = 0;
	actual.z = 1;
	float degreesOff = acosf((actual * bestGuess.c) / (norm(actual) * norm(bestGuess.c)));
	float degreesOffApprox = acosf((actual * s.c) / (norm(actual) * norm(s.c)));
	cout << "After Approximation : <" << s.c.x << " " << s.c.y << " " << s.c.z << ">\n";
	cout << degreesOffApprox * (180.0 / M_PI) << " Degrees Off\n";
	cout << "After Refinement : <" << bestGuess.c.x << " " << bestGuess.c.y << " " << bestGuess.c.z << ">\n";
	cout << degreesOff * (180.0/ M_PI) << " Degrees Off\n";
	
	/// Timing ///

	float timeSec = end_time.tv_sec - start_time.tv_sec;
	float timeuSec = end_time.tv_usec - start_time.tv_usec;
	float timeApprox = approx_time.tv_sec - start_time.tv_sec;
	float timeApproxU = approx_time.tv_usec - start_time.tv_usec;
	float timeRefine = end_time.tv_sec - approx_time.tv_sec;
	float timeRefineU = end_time.tv_usec - approx_time.tv_usec;

	cout << "Time Taken : \nTotal: " << 0.001 * (timeSec * 1000000 + timeuSec) << " Milliseconds \n";
	cout << "Approximation : " << 0.001 * (timeApprox * 1000000 + timeApproxU) << " Milliseconds \n";
	cout << "Refinement : " << 0.001 * (timeRefine * 1000000 + timeRefineU) << " Milliseconds \n";
	cout << "Refinement Kernels : " << 0.001 * (refineSec * 1000000 + refineuSec) << " Milliseconds \n";
	cout << "Linear Algebra : " << 0.001 * (cublasSec * 1000000 + cublasuSec) << " Milliseconds \n";
	cout << "Eigen : " << 0.001 * (eigenSec * 1000000 + eigenuSec) << " Milliseconds \n";

	std::ofstream approxFile;
	std::ofstream refineFile;

	approxFile.open("timings/approxTime.txt", std::ios_base::app); // append instead of overwrite
	refineFile.open("timings/refineTime.txt", std::ios_base::app);
	approxFile << (0.001 * (timeApprox * 1000000 + timeApproxU)) << endl;
	refineFile << (0.001 * (timeRefine * 1000000 + timeRefineU)) << endl;
	approxFile.close();
	refineFile.close();

	return 0;

	
}
