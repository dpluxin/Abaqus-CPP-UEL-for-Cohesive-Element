// Mesh.h
#pragma once

// ******** Definition of Mesh Class ******** //
class Mesh { // use a singleton instance
private:
    Mesh() {} // Private constructor to prevent direct instantiation
    std::mutex mtx;
    std::map<int, std::shared_ptr<Node>> nodes; // Map of nodes by ID
    std::map<int, std::shared_ptr<CohesiveElement>> elements; // Map of elements by ID
public:
    // Delete copy constructor and assignment operator to prevent copies
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    // Singleton instance accessor
    // - static local variable in getInstance, no need to define or manage
    //   the instance explicitly in the global space, instance is automatically
    //   created when getInstance() is called for the first time.
    static Mesh& getInstance() {
        static Mesh instance; // Thread-safe and lazy-initialized
        return instance;
    }

    // definition of node and element map: adding operation is single threaded
    void addNode(int nodeId, const std::vector<double>& coordinates) {
        nodes[nodeId] = std::make_shared<Node>(nodeId, coordinates);
    }
    void addElement(int elementId, const std::vector<int>& nodeIds) {
        std::vector<std::shared_ptr<Node>> elementNodes; // Gather node pointers
        for (int nodeId : nodeIds) {
            if (nodes.find(nodeId) != nodes.end()){
                elementNodes.push_back(nodes[nodeId]);
		    }else {
            std::cerr << "Node ID " << nodeId << " not found!\n";
			}
        }
        elements[elementId] = std::make_shared<CohesiveElement>(elementId, elementNodes);
    }
    // getter of node and element by ID
    std::shared_ptr<Node> getNode(int nodeId) const {
        try {
            return nodes.at(nodeId); // Throws if nodeId is not in the map
        } catch (const std::out_of_range&) {
            return nullptr; // Return nullptr if the nodeId is invalid
        }
    }
    std::shared_ptr<CohesiveElement> getElement(int elementId) const {
        try {
            return elements.at(elementId); // Throws if elementId is not in the map
        } catch (const std::out_of_range&) {
            return nullptr; // Return nullptr if the elementId is invalid
        }
    }

    // implementation function: see the Mesh.inc
    // -- read input file and update mesh components
    void readModelData(char* currentDir);
    // -- output the mesh results to external files
    void outputMesh(char* outDir, int kinc);

    // --------------------------------------------------
    // JUST FOR DEBUGING
    void print();
    void printNodes();
    void printElements();
};

#include "src/Mesh.inc"
