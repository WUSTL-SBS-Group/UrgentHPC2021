# UrgentHPC2021
Public code repository for "Supporting Multi-messenger Astrophysics with Fast Gamma-ray Burst Localization" submitted to the UrgentHPC 2021 workshop.

# Centroiding
This folder contains the code and data used to generate figure 4. fiber_values.txt contains the raw fiber values we centroided, raw_geant.txt contains the ground truth data we compared the centroiding results to, and fiberOut.txt and geantOut.txt hold the same information formatted to be fed into the program via the Python script formatData.py. Building and running main.cpp will output data to three files: output.txt, toMatlabX.txt and toMatlabY.txt. output.txt will contains a list of ambiguous events that were thought to be improperly centroided and the toMatlab files will contain the centroids that are thought to be reasonably correct. The raw geant data can then be compared to the centroids generated in order to find accuracy statistics. 

# GPU_Localization 
This folder contains the code and data used for all GPU localization discussed in the paper. Two bash scripts are provided, the first of which, generateCircles.sh, will generate data to be fed into the localization pipeline. This script creates 200 datasets for 24 fluences (0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, ... , 1.00) that can be found in their respective directories in the dataSets folder. The second script, runTrials.sh, will run the localization pipeline on all the data generated from the previous script and output timing data to the timings folder. Accuracy/Localization data will be output to the console for each run. 
