#include "railmap_helpers.h"

// Function to convert latitude and longitude to Cartesian coordinates (x, y)
std::pair<float, float> convertLatLongToXY(float phi, float psi, float phi_ref, float psi_ref) {

    // Calculate deltaPsi and deltaPhi using the reference latitude and longitude in degrees
    float delta_psi = 111412.84 * cos(phi_ref) - 93.5 * cos(3 * phi_ref) + 0.118 * cos(5 * phi_ref);
    float delta_phi = 111132.92 - 559.82 * cos(2 * phi_ref) + 1.175 * cos(4 * phi_ref) - 0.0023 * cos(6 * phi_ref);

    float x = delta_psi * (psi - psi_ref);
    float y = delta_phi * (phi - phi_ref);

    return {x, y};
}



// Function to convert Cartesian coordinates (x, y) to latitude and longitude
std::pair<float, float> convertXYToLatLong(float x, float y, float phi_ref, float psi_ref) {

    // Calculate delta psi and delta phi using the reference latitude and longitude in degrees
    float delta_psi = 111412.84 * cos(phi_ref) - 93.5 * cos(3 * phi_ref) + 0.118 * cos(5 * phi_ref);
    float delta_phi = 111132.92 - 559.82 * cos(2 * phi_ref) + 1.175 * cos(4 * phi_ref) - 0.0023 * cos(6 * phi_ref);

    // Convert x, y back to latitude (phi) and longitude (psi) using degrees
    float phi = phi_ref + y / delta_phi;
    float psi = psi_ref + x / delta_psi;

    return std::make_pair(phi, psi);
}



// Function to read the rail map
void readRailMap(const std::string filename, std::vector<Node>& nodes, std::vector<Edge>& edges, float &lat_ref, float &lon_ref) {
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
        // Convert latitude and longitude to Cartesian coordinates
        std::pair<double, double> latLongPair = convertLatLongToXY(nodes[i].lat, nodes[i].lon, lat_ref, lon_ref);
        nodes[i].x = latLongPair.first;
        nodes[i].y = latLongPair.second;
    }

    // Read edges
    file >> line >> E;
    edges.resize(E);
    for (int i = 0; i < E; ++i) {
        file >> edges[i].node1 >> edges[i].node2;
    }
}



// Function to read the targets
void readTargets(const std::string filename, std::vector<Node>& targets, float& lat_ref, float& lon_ref) {

    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    int numTargets;
    file >> line >> numTargets;
    targets.resize(numTargets);

    for (int i = 0; i < numTargets; ++i) {
        file >> targets[i].id >> targets[i].lat >> targets[i].lon;
        // Convert latitude and longitude to Cartesian coordinates
        std::pair<double, double> latLongPair = convertLatLongToXY(targets[i].lat, targets[i].lon, lat_ref, lon_ref);
        targets[i].x = latLongPair.first;
        targets[i].y = latLongPair.second;
        
    }
}



// Function to find the closest point on an edge to a target
std::pair<float, float> closestPointOnEdge(float x, float y, const Node node1, const Node node2) {
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


