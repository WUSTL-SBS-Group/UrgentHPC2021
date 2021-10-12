#include "centroid.h"
// As the name suggests, parse initial Fiber Readings into A Hash Map
fiberMap parseFiber(fiberMap map, std::string filename) {
	ifstream inFile(filename);
	string txt;
	getline(inFile,txt);
	
	while (!inFile.eof()) {
		mapKey keyTemp;
		fiberVal valTemp;
		getline(inFile,txt);
		if (inFile.eof()) break;
		keyTemp.id = stoi(txt);
		getline(inFile,txt);
		keyTemp.layer = stoi(txt);
		getline(inFile,txt);
		keyTemp.signal = txt;
		getline(inFile,txt);
		keyTemp.axis = txt;
		getline(inFile,txt);
		valTemp.fiberID = stoi(txt);
		getline(inFile,txt);
		valTemp.x = stof(txt);
		getline(inFile,txt);
		valTemp.z = stof(txt);
		getline(inFile,txt);
		valTemp.energy = stof(txt);
		valTemp.flag = false;
		valTemp.cannotFix = false;
		std::pair<mapKey, fiberVal> tmp {keyTemp, valTemp};
		map.insert(tmp);
	}

	return map;
}
// Take the weighted average of positions sharing the same key and place into new map
fiberMap weightedAverage(fiberMap map, fiberMap averageMap) {
	// Iterate through all entries in the Hash Map
	for (auto x : map) {
		float num = 0, den = 0, localPeak = 0;
		int localMax = 0;
		bool flagged = false, threeMax = false;
		int currID, prevID;
		fiberVal prevVal;

		float z = x.second.z;
		mapKey key = x.first;

		auto range = map.equal_range(key);						// Iterate through every value with the key
		for (auto it = range.first; it != range.second; it++) {
			fiberVal val = it->second;
			num += val.x * val.energy;							// Numerator += relevant position value * energy
			den += val.energy;									// Denominator  += energy
			if (it == range.first) {							// If first value iterated through
				prevID = val.fiberID;
				prevVal = val;
			}
			currID = val.fiberID;
			localPeak = max(localPeak, val.energy);				// Set local Maxes as you go

			// If Non-Contiguous Fiber Values, flag
			if (!(currID + 1 == prevID || currID - 1 == prevID || currID == prevID) && flagged == false && averageMap.find(key) == averageMap.end()) {
			//	cout << "FLAG " << it->first.id  << " " << it->first.layer << " " << it->first.axis << " " << prevID << " " << currID << endl;
				flagged = true;
			}
			// Check/set local maxes
			if (threeMax == false && averageMap.find(key) == averageMap.end()) {
				if (val.energy < prevVal.energy && prevVal.energy == localPeak) localMax++;
				if (val.energy > prevVal.energy) localPeak = val.energy;
				if (localMax > 2) threeMax = true;
			}


			prevID = currID;
			prevVal = val;

		}

		if (averageMap.find(key) == averageMap.end()) {			// If the key is not in the map, insert
			fiberVal average;
			average.x = num / den;
			average.z = z;
			average.energy = den;
			average.fiberID = key.id;

			average.flag = flagged;
			average.cannotFix = threeMax;

			std::pair<mapKey, fiberVal> tmp {key, average};
			averageMap.insert(tmp);
		}
	}

	return averageMap;

}
// Parse Geant Data into a hash map
fiberMap parseRaw(fiberMap map, std::string filename) {
	ifstream inFile(filename);
	string txt;
	getline(inFile,txt);

	while (!inFile.eof() ) {
		int id,layer;
		float x,y,z,e;
		string mat;
		getline(inFile,txt);
		if (inFile.eof()) break;
		id = stoi(txt);
		getline(inFile, txt);
		x = stof(txt);
		getline(inFile,txt);
		y = stof(txt);
		getline(inFile,txt);
		z = stof(txt);
		getline(inFile,txt);
		e = stof(txt);
		getline(inFile,txt);
		getline(inFile,txt);
		layer = stoi(txt);
		getline(inFile,txt);
		if (txt == "GUNout")
			continue;
		else if (txt == "CsI")
			txt = "WLS_Slow";
		else if (txt == "WLS")
			txt = "WLS_Fast";

		mat = txt;
		mapKey xKey, yKey;
		fiberVal xVal, yVal;
		xKey.id = id; yKey.id = id;
		xKey.layer = layer; yKey.layer = layer;
		xKey.signal = mat; yKey.signal = mat;
		xKey.axis = "x"; yKey.axis = "y";
		xVal.x = x; xVal.z = z; xVal.energy = e;
		yVal.x = y; yVal.z = z; yVal.energy = e;
		std::pair<mapKey, fiberVal> tmp1 {xKey, xVal};
		std::pair<mapKey, fiberVal> tmp2 {yKey, yVal};
		map.insert(tmp1);
		map.insert(tmp2);
	}
	return map;
}
// Compare the Maps for accuracy
int compareMaps(fiberMap rawMap, fiberMap averageMap, fiberMap fibers) {
	std::ofstream outFile("data/output.txt");
	std::ofstream matFileX("data/toMatlabX.txt");
	std::ofstream matFileY("data/toMatlabY.txt");
	int zeros = 0, more = 0, flags = 0, flagsTrue = 0, noMap = 0, bad = 0;
	// Iterate through Averaged Map
	for (auto i : averageMap) {
		mapKey key = i.first;
		fiberVal val = i.second;

		if (val.flag == true)
			flags++;								// Increment Number of Flagged Events if Flagged

		int keyCount = rawMap.count(key);
		if (keyCount == 1) {						// If Geant data indicates that this is an unambiguous event, compare averaged value and geant data
			auto search = rawMap.find(key);
			mapKey key2 = search->first;
			fiberVal val2 = search->second;

			if (outFile.is_open()) {
				outFile << "Average: " << key.id << " " << key.layer << " " << key.signal << " " << key.axis << " " << val.x << "\n";
				outFile << "Raw: " << key2.id << " " << key2.layer << " " << key2.signal << " " << key2.axis << " " << val2.x << "\n";
				if (key.axis == "x") {
					matFileX << val.x << "\n" << (abs(val.x - val2.x)) << "\n";
				}
				else {
					matFileY << val.x << "\n" << (abs(val.x - val2.x)) << "\n";
				}
			}

		}
		else {										// If Geant Data Indicates That This is an Ambiguous Event
			if (keyCount > 1) {
				more++;								// Geant Data Has Multiple Hits For a Key
				auto rangeGeant = rawMap.equal_range(key);
				auto rangeFiber = fibers.equal_range(key);
				for (auto it = rangeGeant.first; it != rangeGeant.second; it++)
					outFile << it->second.x <<  endl;
				outFile << "\n";
				for (auto it = rangeFiber.first; it != rangeFiber.second; it++) { 		// Iterate through every value with the key
					outFile << it->second.fiberID << " " << it->second.x <<  " " << it->second.energy << endl;
				}
				outFile << "\n";
			}

			if (keyCount == 0) zeros++;				// Averaged Data Does Not Map to a Geant Event

			if (val.flag == true && (keyCount > 1)) flagsTrue++; 	// Flagged Event Corresponds to Ambiguous Event

			if (val.flag == true && (keyCount == 0)) noMap++;		// Flagged Event Corresponds to no Geant Event

			if (keyCount > 1 && val.cannotFix == true) {			// Multiple Hit Event has >3 Local Maxes
				bad++;
			}

		}
	}
	outFile.close();
	matFileX.close();
	matFileY.close();

	cout << "Flags: " << flags << endl;
	cout << "True Flags: " << flagsTrue << endl;
	cout << "No Mapping: " << noMap << endl;
	cout << "More than Three Peaks : " << bad << endl;
	return more + zeros;
}
