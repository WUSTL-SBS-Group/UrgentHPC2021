#include "centroid.h"
// Driver File for Centroiding Testing, Requires Fiber Data and Geant Data for Comparison. Data Must be Passed into Included Python Script for Correct Formatting.

int main(void) {

	fiberMap fiberM, rawMap, averageMap;
	// Place Fiber Readings into a Map
	fiberM = parseFiber(fiberM, "data\\fiberOut.txt");
	cout << "Fiber Events Read: " << fiberM.size() << endl;
	// Place Geant Data into a Map
	rawMap = parseRaw(rawMap, "data\\geantOut.txt");
	cout << "Events in Geant Data: " << rawMap.size() << endl;
	// Take Weighted Average for Each Key
	averageMap = weightedAverage(fiberM, averageMap);
	cout << "Events after Averaging: " << averageMap.size() << endl;
	// Find Accuracy
	int ambiguous = compareMaps(rawMap, averageMap, fiberM);
	cout << "Ambiguous Events: " << ambiguous << endl;

	return 0;
}
