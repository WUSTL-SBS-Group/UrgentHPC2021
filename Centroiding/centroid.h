/*
 * centroid.h
 *
 *  Created on: Aug 10, 2021
 *      Author: jacwh
 */

#ifndef CENTROID_H_
#define CENTROID_H_

#include <string>
#include <iostream>
#include <unordered_map>
#include <map>
#include <fstream>

using namespace std;

// Key Value Struct for Maps
struct mapKey {
	int id, layer;
	string signal, axis;
	// Need to define Equality Operator for Hashing
	bool operator==(const mapKey &o) const{
		        return id == o.id && layer == o.layer && signal == (o.signal) && axis == (o.axis);
	}
};
// Value Struct for Maps
struct fiberVal {
	int fiberID;
	float x, z, energy;
	bool flag, cannotFix;
};
// Hashing Function
template <class T>
inline void hash_combine(std::size_t & s, const T & v)
{
  std::hash<T> h;
  s^= h(v) + 0x9e3779b9 + (s<< 6) + (s>> 2);
}
// Hashing for mapKey type
typedef std::unordered_multimap<mapKey, fiberVal> fiberMap;
namespace std {
template <>
struct hash<mapKey>
{

    std::size_t operator() (const mapKey &node) const
    {
    	std::size_t res = 0;
    	hash_combine(res,node.id);
    	hash_combine(res,node.axis);
    	hash_combine(res,node.layer);
    	hash_combine(res,node.signal);

        return res;
    }
};
}

fiberMap parseFiber(fiberMap map, std::string filename);

fiberMap weightedAverage(fiberMap map, fiberMap averageMap);

fiberMap parseRaw(fiberMap map, std::string filename);

int compareMaps(fiberMap rawMap, fiberMap averageMap, fiberMap fibers);


#endif /* CENTROID_H_ */
