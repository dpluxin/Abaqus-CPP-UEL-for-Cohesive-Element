// Node.h
#pragma once

// ******** Definition of Node Class ******** //
class Node {
private:
    int id;                             // Unique node identifier
    std::vector<double> coords;         // Coordinates (e.g., x, y, z)
    std::vector<double> displacements;  // Displacement vector
public:
    // constructor
    Node(int nodeId, const std::vector<double>& coordinates)
        : id(nodeId), coords(coordinates), displacements(coordinates.size(), 0.0) {}
    // getter
    int getId() const { return id; }
    const std::vector<double>& getCoordinates() const { return coords; }
    const std::vector<double>& getDisplacements() const { return displacements; }
    // setting functions
    void setDisplacements(const std::vector<double>& disp) {displacements = disp;}
};
