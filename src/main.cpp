#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <limits>

// Define a struct for Nodes and Targets
struct Node {
    int id;
    float lat, lon, x, y; // Include Cartesian coordinates as well
};

struct Edge {
    int node1, node2;
};

struct Target {
    int id;
    float lat, lon, x, y;
};


struct Solution {
    int target_id, id0, id1;
    float lat, lon, distance;
};


// Define global variables
std::vector<Node> nodes;
std::vector<Edge> edges;
std::vector<Target> targets;
std::vector<Solution> solutions;

float lat_ref, lon_ref;


// Function to convert latitude and longitude to Cartesian coordinates (x, y)
std::pair<float, float> convertXYToLatLong(float x, float y, float phi_ref, float psi_ref) {
    // Calculate delta psi and delta phi using the reference latitude and longitude in degrees
    float delta_psi = 111412.84 * cos(phi_ref) - 93.5 * cos(3 * phi_ref) + 0.118 * cos(5 * phi_ref);
    float delta_phi = 111132.92 - 559.82 * cos(2 * phi_ref) + 1.175 * cos(4 * phi_ref) - 0.0023 * cos(6 * phi_ref);

    // Convert x, y back to latitude (phi) and longitude (psi) using degrees
    float phi = phi_ref + y / delta_phi;
    float psi = psi_ref + x / delta_psi;

    return std::make_pair(phi, psi);
}




std::pair<float, float> convertLatLongToXY(float phi, float psi, float phiRef, float psiRef) {
    // Calculate deltaPsi and deltaPhi using the reference latitude and longitude in degrees
    float dPsi = 111412.84 * cos(phiRef) - 93.5 * cos(3 * phiRef) + 0.118 * cos(5 * phiRef);
    float dPhi = 111132.92 - 559.82 * cos(2 * phiRef) + 1.175 * cos(4 * phiRef) - 0.0023 * cos(6 * phiRef);

    // Convert latitude and longitude to Cartesian coordinates
    float x = dPsi * (psi - psiRef);
    float y = dPhi * (phi - phiRef);

    return {x, y};
}


void readRailMap(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    int N, E;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    // Read reference coordinates
    file >> line >> lat_ref >> lon_ref;

    // Read nodes
    file >> line >> N;
    nodes.resize(N);
    for (int i = 0; i < N; ++i) {
        file >> nodes[i].id >> nodes[i].lat >> nodes[i].lon;
        auto [x, y] = convertLatLongToXY(nodes[i].lat, nodes[i].lon, lat_ref, lon_ref); // Structured binding
        nodes[i].x = x;
        nodes[i].y = y;
    }

    // Read edges
    file >> line >> E;
    edges.resize(E);
    for (int i = 0; i < E; ++i) {
        file >> edges[i].node1 >> edges[i].node2;
    }
}

void readTargets(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    int numTargets;
    file >> line >> numTargets;

    std::cout << "Number of targets: " << numTargets << std::endl;

    targets.resize(numTargets);

    for (int i = 0; i < numTargets; ++i) {
        file >> targets[i].id >> targets[i].lat >> targets[i].lon;
        auto [x, y] = convertLatLongToXY(targets[i].lat, targets[i].lon, lat_ref, lon_ref);
        targets[i].x = x;
        targets[i].y = y;
        
    }
}


std::pair<float, float> closestPointOnEdge(float x, float y, const Node& node1, const Node& node2) {
    float x1 = node1.x, y1 = node1.y;
    float x2 = node2.x, y2 = node2.y;

    // Calculate the line segment's lengths squared
    float len_sq = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

    // Check if the line segment is actually a point
    if (len_sq == 0.0) return {x1, y1};

    // Calculate projection of the point onto the line (node1, node2)
    float t = ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)) / len_sq;
    t = std::max(0.0f, std::min(1.0f, t));

    // Get the closest point
    float closestX = x1 + t * (x2 - x1);
    float closestY = y1 + t * (y2 - y1);

    return {closestX, closestY};
}







int main() {

    std::string filename_map = "../inputs/rail_map.txt"; 
    std::string filename_targets = "../inputs/targets.txt";


    readRailMap(filename_map);
    readTargets(filename_targets);

    std::cout  << "Read the railmap" << std::endl;

    // Show the first 3 nodes and 3 edges
    std::cout << "First 3 nodes:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << nodes[i].id << " " << nodes[i].lat << " " << nodes[i].lon << " " << nodes[i].x << " " << nodes[i].y << std::endl;
    }
    std::cout << std::endl;

    std::cout << "First 3 edges:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << edges[i].node1 << " " << edges[i].node2 << std::endl;
    }
    
    std::cout << std::endl;

    std::cout << "First 3 targets:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << targets[i].id << " " << targets[i].lat << " " << targets[i].lon << " " << targets[i].x << " " << targets[i].y << std::endl;
    }

    std::cout << std::endl;


Node node_closest1;
Node node_closest2;
std::pair<float, float> solutionClosestPoint;

// Assuming 'targets' is a vector of pairs or a struct with x, y fields
for (const auto& target : targets) {
    float minDistance = std::numeric_limits<float>::max();

    for (const auto& edge : edges) {
        Node node1 = nodes[edge.node1];
        Node node2 = nodes[edge.node2];

        auto point = closestPointOnEdge(target.x, target.y, node1, node2);
        float distance = sqrt(pow(point.first - target.x, 2) + pow(point.second - target.y, 2));


        if (distance < minDistance) {
            minDistance = distance;
            node_closest1 = node1;
            node_closest2 = node2;
            solutionClosestPoint = convertXYToLatLong(point.first, point.second, lat_ref, lon_ref);
        }
    }



    Solution solution = {
        target.id, 
        node_closest1.id, 
        node_closest2.id, 
        solutionClosestPoint.first, 
        solutionClosestPoint.second,
        minDistance};

    solutions.push_back(solution);

}


    // Show the first 3 solutions
    std::cout << "First 3 solutions:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << solutions[i].target_id << " " << solutions[i].id0 << " " << solutions[i].id1 << " " << solutions[i].lat << " " << solutions[i].lon << " " << solutions[i].distance <<::std::endl;
    }

    // Write output to file
    std::ofstream outFile("../outputs/solutions.txt");
    // Write data, fill with all found solutions
    outFile << "solutions" << " " << solutions.size() <<std::endl;

    for (const auto& solution : solutions) {
        outFile << solution.target_id << " " << solution.id0 << " " << solution.id1 << " " << solution.lat << " " << solution.lon << " " << solution.distance << std::endl;
    }


    outFile.close();

    return 0;
}