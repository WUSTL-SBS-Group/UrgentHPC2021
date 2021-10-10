#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Eigenvalues"
#include "datatypes.cuh"

using namespace Eigen;

typedef class Matrix<double, 6, 6> Matrix6d;

circle refineEigen(matOut outTest) {
	Vector3d g;
	Matrix3d H;
	/// Convert from Column-Major Order ///
	g(0) = static_cast<double>(outTest.g[0]);
	g(1) = static_cast<double>(outTest.g[1]);
	g(2) = static_cast<double>(outTest.g[2]);
	H(0,0) = static_cast<double>(outTest.H[0]);
	H(1,0) = static_cast<double>(outTest.H[1]);
	H(2,0) = static_cast<double>(outTest.H[2]);
	H(1,1) = static_cast<double>(outTest.H[4]);
	H(2,1) = static_cast<double>(outTest.H[5]);
	H(2,2) = static_cast<double>(outTest.H[8]);
	H(0,1) = H(1,0);
	H(0,2) = H(2,0);
	H(1,2) = H(2,1);
	/// Follow Calculations from Marion's Code ///
	
	double scaleFactor = 2/(abs(H.mean()) + abs(g.mean()));
	H *= scaleFactor;
	g *= scaleFactor;

	Matrix3d A0 = H*H - g*g.transpose();

	Matrix6d LA;
	LA << Matrix3d::Zero(), A0, A0, 2*H;

	Matrix6d LB;
	LB << -A0, Matrix3d::Zero(), Matrix3d::Zero(), Matrix3d::Identity();

	GeneralizedEigenSolver<Matrix6d> eigsolver(LA, LB, false);
	assert(eigsolver.info() == Success);

	auto eigs = eigsolver.eigenvalues();


	double lambda = INFINITY;
	for (int j = 0; j < eigs.size(); j++) {
		if (abs(eigs(j).imag()) < 1e-6) 
			lambda = std::min(lambda, eigs(j).real());
	}
	assert(lambda != INFINITY);

	Vector3d xopt = (H - lambda * Matrix3d::Identity()).inverse() * g;

	xopt.stableNormalize();
	circle outCirc;
	outCirc.c.x = xopt(0);
	outCirc.c.y = xopt(1);
	outCirc.c.z = xopt(2);

	return outCirc;
}
