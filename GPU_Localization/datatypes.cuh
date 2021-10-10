#ifndef __DATATYPES_CUH
#define __DATATYPES_CUH

#include </usr/local/cuda/include/cuda_runtime.h>
#include "Vector.h"
struct Hit {
	Vector<float> p;
	float e;

	void set(float E, const Vector<float> &P)
	{
		p = P;
		e = E;
	}
};

typedef float FloatType;

struct Event {
	Hit hits[16];
	int numHits;
	int eventID;
};


struct hitList {
		int p0, p1, p0Orig, p1Orig, IntSet, n;

		FloatType Wprev, dWprevSq, chi2;

		__host__ __device__
		hitList() {}
					
		__host__ __device__
		hitList(int pZ, int pOn, int pZO, int pOnOrig, int set, int en, FloatType wP, FloatType dwP, FloatType chi)
			: p0(pZ),p1(pOn),p0Orig(pZO),p1Orig(pOnOrig),IntSet(set),n(en),Wprev(wP),dWprevSq(dwP),chi2(chi)
		{}
};



struct Triple {
		float eta,dEtaSq;

		__host__ __device__
		Triple() {}
				
		__host__ __device__
		Triple(FloatType et, FloatType det)
			: eta(et),dEtaSq(det)
		{}

};

struct hitOut {
		int p0;
		FloatType chi2;

		__host__ __device__
		hitOut() {}
					
		__host__ __device__
		hitOut(int pZ, FloatType chi)
			: p0(pZ),chi2(chi)
		{}
};

struct circle {
		Vector<float> c;
		float phi, dphi;

		__host__ __device__
		circle() {}

		__host__ __device__
		circle(Vector<float> cin, float Phi, float dPhi)
			: c(cin),phi(Phi),dphi(dPhi)
		{}
};

struct approxOut {
	Vector<float> c;
	float LL;

	__host__ __device__
	approxOut() {}

	__host__ __device__
	approxOut(Vector<float> cin, float ll)
		: c(cin),LL(ll)
	{}
};

struct matOut {
	float *g, *H;

	__host__ __device__
	matOut() {}

	__host__ __device__
	matOut(float *G, float *h)
		: g(G),H(h)
	{}
};

#endif
