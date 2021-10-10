#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <list>

std::random_device dev;
std::mt19937 prng(dev());

using namespace std;
using namespace std::chrono;

int main(int argc, char** argv) {
	string circlesIn = argv[1];
	string circlesOut = argv[2];	
	int photons = atoi(argv[3]);

	ifstream inFile(circlesIn);
	ofstream outFile(circlesOut);
	string txt; 

	uniform_int_distribution<int> D(0,999999);
	list<string> *circleList = new list<string>[1000000];
	int *photonListAll = new int[1000000];
	int *photonList = new int[photons];
	for (int i = 0; i < 1000000; i++) {
		photonListAll[i] = i;
	}
	// Insert Read Data into Array of String Lists
	while (!inFile.eof()) {
		getline(inFile,txt);
		bool breakFlag = false;
		string id = "";
		string circleData = "";
		for (int i = 0; i < txt.length(); i++) {
			if (!breakFlag && txt[i] != ' ') id += txt[i];
			if (breakFlag) circleData += txt[i];
			if (txt[i] == ' ' && breakFlag == false) breakFlag = true;	
		}
		if (!inFile.eof()) circleList[stoi(id)].push_front(circleData);
	}
	
	// Fill list of photons with random ids
	for (int i = 0; i < 200; i++) {
		long ms = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
		prng.seed(ms);
		ofstream outFile(circlesOut + to_string(i) + ".txt");
		if (outFile.is_open()) outFile << "\n";
		for (int j = 0; j < photons; j++) {
			swap(photonList[j], photonListAll[D(prng)]);
		}
		for (int k = 0; k < photons; k++) {
			for (auto it = circleList[photonList[k]].begin(); it != circleList[photonList[k]].end(); it++) {
				if (outFile.is_open()) {
					string data = it->c_str();
					for (int j = 0; j < data.length(); j++) {
						if (data[j] != ' ') outFile << data[j];
						else outFile << "\n";
					}
					outFile << "\n\n";
				}
			}
		}
	}
	outFile.close();
	cout << "done" << endl;
}
