// uel.h
// declare the classes, variables and functions used in uel.cpp
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Watch out the pointer operations:
        Copy of pointer or object??
++++++++++++++++++++++++++++++++++++++++++++++++++++++  */
#pragma once
#include <string.h>
#include <vector>
#include <memory>
#include <map>
#include <mutex>

// ******** Definition of Parameters ******** d//
#define UELNAME "U"
#define NDOF    3
#define NSTR    3
#define CHARLEN 1024
#define SMALLNUM 0.00000001

// ******** Definition of UEL variables ******** //
class UelVariables{ // some global variable (operate with care)
public:
    static double globalTime;
    static char   workFolder[CHARLEN];
    static char   outFolder[CHARLEN];
};

// ******** Definition of Classes ******** //
#include "src/ToolKit.h"
#include "src/Node.h"
#include "src/CohesiveMaterial.h"
#include "src/Element.h"
#include "src/CohesiveElement.h"
#include "src/Mesh.h"
