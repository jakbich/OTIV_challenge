#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

// Define a struct for Nodes and Targets
struct Node {
    int id;
    double lat, lon, x, y; // Include Cartesian coordinates as well
};

struct Edge {
    int node1, node2;
};

// Define global variables
std::vector<Node> nodes;
std::vector<Edge> edges;
double lat_ref, lon_ref;


// Function to convert degrees to radians
double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

// Function to convert latitude and longitude to Cartesian coordinates (x, y)
std::pair<double, double> convertToXY(double phi, double psi, double phiRef, double psiRef) {

    // Convert degrees to radians
    double phiRad = deg2rad(phi);
    double psiRad = deg2rad(psi);
    double phiRefRad = deg2rad(phiRef);
    double psiRefRad = deg2rad(psiRef);

    // Calculate deltaPsi
    double dPsi = 111412.84 * cos(phiRefRad) - 93.5 * cos(3 * phiRefRad) + 0.118 * cos(5 * phiRefRad);

    // Calculate deltaPhi
    double dPhi = 111132.92 - 559.82 * cos(2 * phiRefRad) + 1.175 * cos(4 * phiRefRad) - 0.0023 * cos(6 * phiRefRad);

    // Convert latitude and longitude to Cartesian coordinates
    double x = dPsi * (psiRad - psiRefRad);
    double y = dPhi * (phiRad - phiRefRad);

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

    // Read reference point
    double lat_ref, lon_ref;
    file >> line >> lat_ref >> lon_ref;

    // Read nodes
    file >> line >> N;
    nodes.resize(N);
    for (int i = 0; i < N; ++i) {
        file >> nodes[i].id >> nodes[i].lat >> nodes[i].lon;
        auto [x, y] = convertToXY(nodes[i].lat, nodes[i].lon, lat_ref, lon_ref); // Structured binding
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

// Function to convert lat/lon to Cartesian coordinates
std::pair<double, double> convertToCartesian(double lat, double lon) {
    // Implement conversion logic
}

// Function to find closest point on an edge
std::pair<double, double> closestPointOnEdge(/* Parameters */) {
    // Implement logic to find the closest point on an edge
}

int main() {

    std::string filename = "../inputs/rail_map.txt"; 

    readRailMap(filename);

    std::cout  << "Read the railmap" << std::endl;

    // Show the first 3 nodes and 3 edges
    std::cout << "First 3 nodes:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << nodes[i].id << " " << nodes[i].lat << " " << nodes[i].lon << " " << nodes[i].x << " " << nodes[i].y << std::endl;
    }

    std::cout << "First 3 edges:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << edges[i].node1 << " " << edges[i].node2 << std::endl;
    }
    


    // // Process each target
    // for (const auto& target : targets) {
    //     // Find the closest point on the rail network
    //     for (const auto& edge : edges) {
    //         // Calculate closest point on this edge
    //         auto closestPoint = closestPointOnEdge(/* Parameters */);
    //         // Update if this is the closest so far
    //     }
    // }

    // Write output to file
    std::ofstream outFile("outputs/solutions.txt");
    // Write data
    outFile.close();

    return 0;
}