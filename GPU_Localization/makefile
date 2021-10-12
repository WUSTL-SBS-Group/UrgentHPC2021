refineEigen: refineEigen.cpp
	g++ -O3 -march=native -c -o refineEigen.o refineEigen.cpp
postCompton: refineEigen.cpp approx.cu refine.cu postCompton.cu matrix.cu
	nvcc -O3 -arch=compute_75 -code=sm_75 -o  postCompton postCompton.cu refineEigen.o approx.cu refine.cu 
	g++ -o shuffler shuffleCircles.cpp
