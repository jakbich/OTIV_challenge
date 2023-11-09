#include "railmap_helpers.h"
#include <fstream>
#include <string>
#include <limits>
#include <iostream>
#include <cmath>

int main() {

// Define variables (see data structures in header file)
std::vector<Node> nodes;
std::vector<Node> targets;
std::vector<Edge> edges;
std::vector<Solution> solutions;
float lat_ref, lon_ref;

// Read input files
std::string filename_map = "../inputs/rail_map.txt"; 
std::string filename_targets = "../inputs/targets.txt";
readRailMap(filename_map, nodes, edges, lat_ref, lon_ref);
readTargets(filename_targets, targets, lat_ref, lon_ref);

// Define variables to store closest nodes and closest point during loop
Node node_closest1;
Node node_closest2;
std::pair<float, float> solution_closest_point;

// Iterate through all targets
for (const auto& target : targets) {
    float min_distance = std::numeric_limits<float>::max();

    // Iterate through all edges
    for (const auto& edge : edges) {
        Node node1 = nodes[edge.node1];
        Node node2 = nodes[edge.node2];

        // Find the closest point on the edge to the target
        auto point = closestPointOnEdge(target.x, target.y, node1, node2);
        float distance = sqrt(pow(point.first - target.x, 2) + pow(point.second - target.y, 2));

        // Update the minimum distance and closest nodes if necessary
        if (distance < min_distance) {
            min_distance = distance;
            node_closest1 = node1;
            node_closest2 = node2;
            solution_closest_point = convertXYToLatLong(point.first, point.second, lat_ref, lon_ref);
        }
    }

    // Create solution
    Solution solution = {
        target.id, 
        node_closest1.id, 
        node_closest2.id, 
        solution_closest_point.first, 
        solution_closest_point.second,
        min_distance};

    solutions.push_back(solution);
}

    // Write data, fill with all found solutions
    std::ofstream outFile("../outputs/solutions.txt");
    outFile << "solutions" << " " << solutions.size() <<std::endl;

    for (const auto& solution : solutions) {
        outFile << solution.target_id << " " << solution.id0 << " " << solution.id1 << " " << solution.lat << " " << solution.lon << " " << solution.distance << std::endl;
    }

    outFile.close();
    return 0;
}