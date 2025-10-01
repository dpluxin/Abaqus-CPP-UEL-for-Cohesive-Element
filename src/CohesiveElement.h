// CohesiveElement.h
#pragma once

// ******** Definition of CohesiveElement Class ******** //
class CohesiveElement : public Element {
private:
    int nnode;  // number of nodes of this element: can be 6 (coh3d6) or 8 (coh3d8)
    double localTime; // local time variable to check convergence
    std::vector<std::shared_ptr<Node>> nodes; // List of nodes in the element
    std::vector<std::shared_ptr<SolutionVariable>> conSdvs; // List of converged sdv
    std::vector<std::shared_ptr<SolutionVariable>> itrSdvs; // List of iterating sdv
    //--- number of sdv: coh3d6->3, coh3d8->4
public:
    // constructor, need to call base-class constructor 'Element(elementId)'
    CohesiveElement(int elementId, const std::vector<std::shared_ptr<Node>>& elementNodes)
    : Element(elementId), nnode(elementNodes.size()), localTime(0.0), nodes(elementNodes) {
        for(int i = 0; i < nnode/2; i++){
            conSdvs.push_back(std::make_shared<SolutionVariable>());
            itrSdvs.push_back(std::make_shared<SolutionVariable>());
        }
    }

    // implement Element interface
    // getter
    int getNnode() const override { return nnode; }
    const std::vector<std::shared_ptr<Node>>& getNodes() const override {return nodes;}
    double getLocalTime() const override {return localTime;}
    // setter functions
    void setLocalTime(double inputTime) override {localTime = inputTime;}

    // implementation function: see the CohesiveElement.inc
    void integrate(int kinc, int elemID, int nn, double xc[], double yc[], double zc[],
        double disp[], double paramsR[], double stf[], double qf[], double dt) override;

    // local getter implementation
    const std::vector<std::shared_ptr<SolutionVariable>>& getSdvs(char flag) const {
        if(flag == 'C') return conSdvs; // get the converged sdv
        if(flag == 'I') return itrSdvs; // get the iterating sdv
        throw std::invalid_argument("Invalid flag. Expected 'C' or 'I'.");
    }

    // local integration function for specific element types
    void kq_coh3d8(double xc[], double yc[], double zc[], double disp[],
        double paramsR[], double stf[], double qf[], double dt);
    void kq_coh3d6(double xc[], double yc[], double zc[], double disp[],
        double paramsR[], double stf[], double qf[], double dt);
};

#include "src/CohesiveElement.inc"
