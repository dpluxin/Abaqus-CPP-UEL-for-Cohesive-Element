// Mesh.inc
#pragma once
// include file operation function
#include <unistd.h>    // close, lseek
#include <dirent.h>    // directory
#include <fcntl.h>     // open
#include <sys/mman.h>  // mmap, munmap
#include <sys/stat.h>  // fstat
#include <fstream>     // For std::ofstream
// string related
#include <sstream>     // For std::istringstream
#include <cstring>     // For strcpy ..
#include <ctype.h>     // for isblank

// ******** Template Functions of Mesh Class ******** //
/*
    *** Rules of abaqus naming and numbering:
    1. Parts with UEL should be named as "auelPart1", "auelPart2" ...
    2. For each part, the node and element labels start from 1 in inp,
       but UEL will use global numbering.
    3. If one part has other built-in elements, the element label in UEL
       across different parts can be DIScontinuous, e.g.,
            auelPart1:   1~100 (UEL)    local labels: 1~100
                       101~200 (C3D8)
            auelPart2: 201~300 (UEL)    local labels: 1~100
       This will cause MISMATCH of element labels!!!

    *** schemes to define abaqus UEL model
    1. Currently, "UEL part with UEL only" is allowed.
    2. The UEL parts should be defined in sequence, so that this code
       can extract the number the mesh variable in order.
    3. A Python script is firstly used to update the UEL declaration,
       in this script the following is achieved:
            The local node and element labels will be updated to the
            global ones considered shifting of index.
               > a global counter for nodal and elemental index is used
               > the last UEL integer props is defined for debuging
                 auelPart1: JPROPS[-1] = 1, auelPart2: JPROPS[-1] = 2
    4. This CPP code will be called at the beginning of the analysis,
       and it will read all "auelPart" and create node/element entities.
*/

// read data from inp file and set the mesh variables
void Mesh::readModelData(char* currentDir){

    // define the file name
    char  filePath[CHARLEN];
    strcpy(filePath, currentDir); strcat(filePath, "/");

    // Open the directory for reading
    DIR* dir = opendir(currentDir);
    if (dir == NULL) perror("opendir() error");

    struct dirent* entry;
    // The loop continues until there are no more entries (readdir() returns NULL)
    while ((entry = readdir(dir)) != NULL) {
        // Check if the file ends with ".inp", entry->d_name file name component
        if (strstr(entry->d_name, ".inp") != NULL &&
        strcmp(entry->d_name + strlen(entry->d_name) - 4, ".inp") == 0) {
            // Moves the pointer to the last 4 characters of the string.
            strcat(filePath, entry->d_name);
        }
    }
    // Close the directory
    closedir(dir);

    // read the file based on memory map
    int file_handle = open(filePath, O_RDONLY);
    if (file_handle == -1) std::cerr << "Error opening file." << std::endl;
    // Get the size of the file
    struct stat file_stat;
    if (fstat(file_handle, &file_stat) == -1) std::cerr << "Error getting size." << std::endl;
    size_t file_size = file_stat.st_size;
    // Map the file into memory
    char* file_data = static_cast<char*>(mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, file_handle, 0));
    if (file_data == MAP_FAILED) std::cerr << "Error mapping view of file." << std::endl;

    // process inp file data and set up the mesh map variable
    std::istringstream file_stream(std::string(file_data, file_size));
    std::string line;
    bool parsingNodes = false;
    bool parsingElements = false;

    while (std::getline(file_stream, line)) {

        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n")); // remove leading whitespace
        line.erase(line.find_last_not_of(" \t\r\n") + 1); // remove trailing whitespace

        if (line.empty() || line[0] == '*') { // string is empty or first character is '*'
            // Handle section headers
            if (line == "*Node") {
                parsingNodes = true;
                parsingElements = false;
            } else if (line.find("*ELEMENT") != std::string::npos &&
                       line.find(UELNAME)    != std::string::npos) {
                parsingNodes = false;
                parsingElements = true;
            } else {
                if (line.find("end of auelPart") != std::string::npos) break;
                parsingNodes = false;
                parsingElements = false;
            }
            continue;
        }

        if (parsingNodes) {
            // Parse a node line
            std::istringstream node_stream(line);
            int nodeId;
            double x, y, z;
            char comma;
            if (node_stream >> nodeId >> comma >> x >> comma >> y >> comma >> z) {
                addNode(nodeId, {x, y, z});
            }
        } else if (parsingElements) {
            // Parse an element line
            std::istringstream element_stream(line);
            int elementId, nodeId;
            char comma;
            std::vector<int> nodeIds;
            if (element_stream >> elementId) {
                while (element_stream >> comma >> nodeId) {
                    nodeIds.push_back(nodeId);
                }
                addElement(elementId, nodeIds);
            }
        }

    } // end of while (std::getline(file_stream, line))

    // Unmap and close the file
    munmap(file_data, file_size);
    close(file_handle);

}

// write output files for results visualization
void Mesh::outputMesh(char* outDir, int kinc){

    // output file (path+name)
    char filePath[CHARLEN];

    // convert inc no to string and include in file name
    strcpy(filePath, outDir); strcat(filePath, "CE-");
    // -- include increment in file name
    const int desiredLength = 6;
    std::string incStr = std::to_string(kinc);
    // Add leading zeros
    while (incStr.length() < desiredLength) {
        incStr.insert(0, "0");
    }
    incStr += ".vtk";
    strcat(filePath, incStr.c_str());

    // create or open the file
    std::ofstream vtkFile(filePath, std::ios::out); // create new file or overwrite
    if (!vtkFile.is_open()) std::cerr << "Error: Unable to open file " << filePath << std::endl;

    // VTK Header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Cohesive Element Mesh\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";

    // Write nodes
    vtkFile << "POINTS " << nodes.size() << " double\n";
    for (const auto& [id, node] : nodes) {
        const auto& coords = node->getCoordinates();
        vtkFile << coords[0] << " " << coords[1] << " " << (coords.size() > 2 ? coords[2] : 0.0) << "\n";
    }

    // Write elements
    int totalConnectivity = 0;
    for (const auto& [id, element] : elements) {
        totalConnectivity += element->getNodes().size() + 1; // Element node count + 1 for VTK format
    }

    vtkFile << "\nCELLS " << elements.size() << " " << totalConnectivity << "\n";
    for (const auto& [id, element] : elements) {
        const auto& elementNodes = element->getNodes();
        vtkFile << elementNodes.size(); // Number of nodes in the element
        for (const auto& node : elementNodes) {
            vtkFile << " " << node->getId() - 1; // VTK uses 0-based indexing
        }
        vtkFile << "\n";
    }

    // Write cell types
    vtkFile << "\nCELL_TYPES " << elements.size() << "\n";
    for (const auto& [id, element] : elements) {
        vtkFile << (element->getNodes().size() == 6 ? 13 : 12) << "\n"; // 13 for VTK_WEDGE, 12 for VTK_HEXAHEDRON
    }

    // Write point data (displacements)
    vtkFile << "\nPOINT_DATA " << nodes.size() << "\n";
    vtkFile << "VECTORS displacement double\n";
    for (const auto& [id, node] : nodes) {
        const auto& disp = node->getDisplacements();
        vtkFile << disp[0] << " " << disp[1] << " " << (disp.size() > 2 ? disp[2] : 0.0) << "\n";
    }

    // Write cell data (separation, traction, damage)
    vtkFile << "\nCELL_DATA " << elements.size() << "\n";

    // Separation
    vtkFile << "VECTORS separation double\n";
    for (const auto& [id, element] : elements) {
        const auto& sdvs = element->getSdvs('I');
        std::vector<double> avgSeparation(3, 0.0); // Assuming 3D vectors
        int numIntegrationPoints = sdvs.size();
        for (const auto& sdv : sdvs) {
            for (size_t i = 0; i < avgSeparation.size(); ++i) {
                avgSeparation[i] += sdv->separ[i];
            }
        }
        for (double& val : avgSeparation) {
            val /= numIntegrationPoints;
        }
        vtkFile << avgSeparation[0] << " " << avgSeparation[1] << " " << avgSeparation[2] << "\n";
    }

    // Traction
    vtkFile << "VECTORS traction double\n";
    for (const auto& [id, element] : elements) {
        const auto& sdvs = element->getSdvs('I');
        std::vector<double> avgTraction(3, 0.0); // Assuming 3D vectors
        int numIntegrationPoints = sdvs.size();
        for (const auto& sdv : sdvs) {
            for (size_t i = 0; i < avgTraction.size(); ++i) {
                avgTraction[i] += sdv->tract[i];
            }
        }
        for (double& val : avgTraction) {
            val /= numIntegrationPoints;
        }
        vtkFile << avgTraction[0] << " " << avgTraction[1] << " " << avgTraction[2] << "\n";
    }

    // Damage
    vtkFile << "SCALARS damage double\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (const auto& [id, element] : elements) {
        const auto& sdvs = element->getSdvs('I');
        double avgDamage = 0.0;
        int numIntegrationPoints = sdvs.size();
        for (const auto& sdv : sdvs) {
            avgDamage += sdv->damage;
        }
        avgDamage /= numIntegrationPoints;
        vtkFile << avgDamage << "\n";
    }

    // FailureStatus
    vtkFile << "SCALARS FailureStatus double\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (const auto& [id, element] : elements) {
        const auto& sdvs = element->getSdvs('I');
        SolutionVariable::Status FailureStatus = SolutionVariable::Intact;
        for (const auto& sdv : sdvs) {
            FailureStatus = std::max(FailureStatus, sdv->currentStatus);
        }
        vtkFile << FailureStatus << "\n";
    }

    vtkFile.close();
}

void Mesh::print(){
    // debug printing
    /*
    std::lock_guard<std::mutex> lock(mtx);
    std::cout << "test started" << std::endl;
    for(int i = 0; i <10; i++) std::cout << "test id" << i;
    std::cout << "test finished" << std::endl;
    */
    std::cout << "node size" << nodes.size() << std::endl;
    std::cout << "elem size" << elements.size() << std::endl;

}
void Mesh::printNodes() {
    std::lock_guard<std::mutex> lock(mtx); // Ensure thread safety
    std::cout << "Nodes in the mesh:" << std::endl;
    for (const auto& [id, node] : nodes) {
        std::cout << "Node ID: " << id << ", Coordinates: (";
        const auto& coords = node->getCoordinates();
        for (size_t i = 0; i < coords.size(); ++i) {
            std::cout << coords[i];
            if (i < coords.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }
}
void Mesh::printElements() {
    std::lock_guard<std::mutex> lock(mtx); // Ensure thread safety
    std::cout << "Elements in the mesh:" << std::endl;
    for (const auto& [id, element] : elements) {
        std::cout << "Element ID: " << id << ", Connected Nodes: ";
        const auto& elementNodes = element->getNodes();
        for (size_t i = 0; i < elementNodes.size(); ++i) {
            std::cout << elementNodes[i]->getId();
            if (i < elementNodes.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
}
