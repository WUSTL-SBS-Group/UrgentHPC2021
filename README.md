# UrgentHPC2021
Public code repository for "Supporting Multi-messenger Astrophysics with Fast Gamma-ray Burst Localization" submitted to the UrgentHPC 2021 workshop.

# Centroiding
This folder contains the code and data used to generate figure 4. fiber_values.txt contains the raw fiber values we centroided, raw_geant.txt contains the ground truth data we compared the centroiding results to, and fiberOut.txt and geantOut.txt hold the same information formatted to be fed into the program via the Python script formatData.py. Running building and running main.cpp will output data to three files: output.txt, toMatlabX.txt and toMatlabY.txt. output.txt will contains a list of ambiguous events that were thought to be improperly centroided and the toMatlab files will contain the centroids that are thought to be reasonably correct. The raw geant data can then be compared to the centroids generated in order to find accuracy statistics. 

# 
