#include "railmap_helpers.h"

std::pair<float, float> readReferenceCoordinates(const std::string filename) {
/**
 * Reads the reference coordinates from the file.
 * 
 * @param filename The name of the file containing the reference coordinates.
 * @return A pair containing the reference latitude and longitude.
 * 
*/
    std::ifstream file(filename);
    std::string line;
    float lat_ref, lon_ref;

    // Check if file is open
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    // Read reference coordinates
    file >> line >> lat_ref >> lon_ref;
    file.close();

    return {lat_ref, lon_ref};
}



std::pair<float, float> convertLatLongToXY(float phi, float psi, float lat_ref, float lon_ref, float delta_phi, float delta_psi) {
/**
 * Converts latitude and longitude to Cartesian coordinates.
 * 
 * @param phi The latitude in degrees.
 * @param psi The longitude in degrees.
 * @param lat_ref The reference latitude in degrees.
 * @param lon_ref The reference longitude in degrees.
 * @param delta_phi The length of one degree of latitude in meters at reference point.
 * @param delta_psi The length of one degree of longitude in meter at reference point.
 * 
 * @return A pair containing the x and y coordinates.
*/

    float x = delta_psi * (psi - lon_ref);
    float y = delta_phi * (phi - lat_ref);

    return {x, y};
}



std::pair<float, float> convertXYToLatLong(float x, float y, float lat_ref, float lon_ref, float delta_phi, float delta_psi) {
/**
 * Converts Cartesian coordinates to latitude and longitude.
 * 
 * @param x The x-coordinate.
 * @param y The y-coordinate.
 * @param lat_ref The reference latitude in degrees.
 * @param lon_ref The reference longitude in degrees.
 * @param delta_phi The length of one degree of latitude in meters at reference point.
 * @param delta_psi The length of one degree of longitude in meter at reference point.
 * 
 * @return A pair containing the latitude and longitude in degrees.
*/
    float phi = lat_ref + y / delta_phi;
    float psi = lon_ref + x / delta_psi;

    return {phi, psi};
}



void readRailMap(const std::string filename, std::vector<Node>& nodes, std::vector<Edge>& edges, float lat_ref, float lon_ref, float delta_phi, float delta_psi) {
/**
 * Reads the rail map from the file.
 * 
 * @param filename The name of the file containing the rail map.
 * @param nodes A vector to store the nodes.
 * @param edges A vector to store the edges.
 * @param lat_ref The reference latitude in degrees.
 * @param lon_ref The reference longitude in degrees.
 * @param delta_phi The length of one degree of latitude in meters at reference point.
 * @param delta_psi The length of one degree of longitude in meter at reference point. * 
*/

    std::ifstream file(filename);
    std::string line;
    int N, E;

    // Check if file is open
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    // Skip the first line which contains reference coordinates
    std::getline(file, line);

    // Read nodes
    file >> line >> N;
    nodes.resize(N);
    for (int i = 0; i < N; ++i) {
        file >> nodes[i].id >> nodes[i].lat >> nodes[i].lon;
        // Convert latitude and longitude to Cartesian coordinates
        std::pair<double, double> latLongPair = convertLatLongToXY(nodes[i].lat, nodes[i].lon, lat_ref, lon_ref, delta_phi, delta_psi);
        nodes[i].x = latLongPair.first;
        nodes[i].y = latLongPair.second;
    }

    // Read edges
    file >> line >> E;
    edges.resize(E);
    for (int i = 0; i < E; ++i) {
        file >> edges[i].node1 >> edges[i].node2;
    }

    file.close();
}



void readTargets(const std::string filename, std::vector<Node>& targets, float lat_ref, float lon_ref, float delta_phi, float delta_psi) {
/**
 * Reads the targets from the file.
 * 
 * @param filename The name of the file containing the targets.
 * @param targets A vector to store the targets.
 * @param lat_ref The reference latitude in degrees.
 * @param lon_ref The reference longitude in degrees.
 * @param delta_phi The length of one degree of latitude in meters at reference point.
 * @param delta_psi The length of one degree of longitude in meter at reference point. 
*/

    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    // Read targets
    int numTargets;
    file >> line >> numTargets;
    targets.resize(numTargets);

    for (int i = 0; i < numTargets; ++i) {
        file >> targets[i].id >> targets[i].lat >> targets[i].lon;
        // Convert latitude and longitude to Cartesian coordinates
        std::pair<double, double> latLongPair = convertLatLongToXY(targets[i].lat, targets[i].lon, lat_ref, lon_ref, delta_phi, delta_psi);
        targets[i].x = latLongPair.first;
        targets[i].y = latLongPair.second;
        
    }
    file.close();
}




std::pair<float, float> closestPointOnEdge(float x, float y, const Node node1, const Node node2) {
/**
 * Finds the closest point on an edge to a given point.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param node1 The first node of the edge.
 * @param node2 The second node of the edge.
 * 
 * @return A pair containing the x and y coordinates of the closest point.
*/

    float x1 = node1.x, y1 = node1.y;
    float x2 = node2.x, y2 = node2.y;

    // Calculate the line segment's lengths squared
    float len_sq = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

    // Check if the line segment is not actually a point
    if (len_sq == 0.0) return {x1, y1};

    // Calculate projection of the point onto the line (node1, node2)
    float t = ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)) / len_sq;
    t = std::max(0.0f, std::min(1.0f, t));

    // Get the closest point
    float closestX = x1 + t * (x2 - x1);
    float closestY = y1 + t * (y2 - y1);

    return {closestX, closestY};
}


