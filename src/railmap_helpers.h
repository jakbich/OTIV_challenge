#ifndef RAIL_MAP_HELPERS_H
#define RAIL_MAP_HELPERS_H

#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <iostream>
#include <fstream>

// Define data structures
struct Node {
    int id;
    float lat, lon, x, y; 
};
struct Edge {
    int node1, node2;
};
struct Solution {
    int target_id, id0, id1;
    float lat, lon, distance;
};

// Function declarations
std::pair<float, float> convertLatLongToXY(float phi, float psi, float phiRef, float psiRef);
std::pair<float, float> convertXYToLatLong(float x, float y, float phi_ref, float psi_ref);
void readRailMap(const std::string filename, std::vector<Node>& nodes, std::vector<Edge>& edges, float &lat_ref, float &lon_ref);
void readTargets(const std::string filename, std::vector<Node>& targets, float& lat_ref, float& lon_ref);
std::pair<float, float> closestPointOnEdge(float x, float y, const Node node1, const Node node2);

#endif // RAIL_MAP_HELPERS_H
