// CohesiveMaterial.h
#pragma once

// ******** Definition of Solution Variable Class ******** //
class SolutionVariable {
public:
    enum Status {           // Enum to represent the damage status
        Intact,             // 0
        DamageInitiation,   // 1
        FinalFailure        // 2
    };
    Status currentStatus = Intact;   // Variable to hold the current status
    double failIndex = 0.0, damage = 0.0, u0 = 0.0, t0 = 0.0, uf = 0.0;
    std::vector<double> tract = std::vector<double>(NSTR, 0.0);
    std::vector<double> separ = std::vector<double>(NSTR, 0.0);
public:
    SolutionVariable() = default;
};

// ******** Definition of CohesiveMaterial Class ******** //
/*  CohesiveMaterial class defines the procedures to calculate the materials tangent
    stiffness matrix and update the SolutionVariable components. The CohesiveMaterial
    instance is NOT required -> static function is defined, and the SolutionVariable
    instance is passed in for updating process (passed-in as a shared pointer).
        Definition of properties and state variables
        --------------------------------------------
        nprops = 10
        k    = para[0];   // penalty stiffness
        t_nc = para[1];   // strengths
        t_sc = para[2];
        Gnc  = para[3];   // fracture toughness
        Gsc  = para[4];
        alpha= para[5];   // mixed mode coefficient
        neta = para[6];   // viscosity parameter
        cof  = para[7];   // compression enhanced mode II fracture
        DMAX = para[8];   // maximum damage, default = 1.0
        ipr  = para[9];   // power law flag: 0-> power, 1->BK
*/
class CohesiveMaterial{
public:
    static void DDSDDE(double u1, double u2, double u3, double dmatrx[NDOF][NDOF],
        double s[], double para[], std::shared_ptr<SolutionVariable> sdv, double dt);

};

#include "src/CohesiveMaterial.inc"
