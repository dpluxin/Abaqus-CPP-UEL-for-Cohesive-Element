/*--------------------------------------------------
++++++++++++++++++++++++++++++++++++++++++++++++++++
            UEL encapsulated by CPP
    The code is verified in Linux with GCC compiler
         g++ (GCC) 8.x supporting c++ 17
           01/12/2024 by LU XIN @Utokyo
++++++++++++++++++++++++++++++++++++++++++++++++++++
            Important Comments
    1. Rules should be followed to create ABAQUS model
       -> ref. Mesh class (15/01/2025)
--------------------------------------------------*/

// include system head files
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
// SYS API cpu information
//#include <sys/syscall.h>
//#include <sys/types.h>

// Fortran-to-C interface interoperability
#include <aba_for_c.h>
// include user-define head files
#include "uel.h"

/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DECLARE THE GLOBAL VARIABLES OF THE CLASS
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
// initialize global variable in UelVariables class
double UelVariables::globalTime = 0.0;        // Initialize 'globalTime'
char UelVariables::workFolder[CHARLEN] = "";  // Initialize 'workFolder'
char UelVariables::outFolder[CHARLEN] = "";   // Initialize 'outFolder'
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        END DECLARE THE GLOBAL VARIABLES OF THE CLASS
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */


//=============== UEXTERNALDB ===============//
extern "C" void FOR_NAME(uexternaldb,UEXTERNALDB)(int* LOP, int* LRESTART,
    double* TIME, double* DTIME, int* KSTEP, int* KINC){

    // start of the analysis
    if(*LOP == 0){

        // initialize global variable
        //--- get the current working folder
        if (getcwd(UelVariables::workFolder, sizeof(UelVariables::workFolder)) == NULL)
            perror("getcwd() error"); // getcwd
        char* lastSlash = strrchr(UelVariables::workFolder, '/'); // Find the last / in the string
        if (lastSlash != NULL) *lastSlash = '\0';  // Null-terminate at the last /
        //--- set the output folder
        strcpy(UelVariables::outFolder, UelVariables::workFolder);
        strcat(UelVariables::outFolder, "/outputs/");

        //--- read the model data
        // the function is defined in Mesh class
        Mesh::getInstance().readModelData(UelVariables::workFolder);

    }else if(*LOP == 1){ // start of increment
        // set the global time: converged time without time increment
        UelVariables::globalTime = TIME[1];

    }else if(*LOP == 2){ // end of increment
        Mesh::getInstance().outputMesh(UelVariables::outFolder, *KINC);
    }else if(*LOP == 3){ // end of analysis
        //readModelData(UelVariables::workFolder); // caused error when read data at the end of the analysis
    }

}

//=============== UEL ===============//
/*
    ON THE INPUT VARIABLES
        The passed-in variables are declared pointers. The original fortran
    code uses COLUMN-major order, e.g. coords(mcrd,nnode), and when being
    passed in here as a 1D array, the 2D matrix follows the column order.
    COORDS[0] = coords(x,node1), COORDS[1] = coords(y,node1) ...

    ON THE LOCAL VARIABLES
        For local arrays, the CPP built-in array is used for efficiency.
*/
extern "C" void FOR_NAME(uel,UEL)(double* RHS, double* AMATRX, double* SVARS,
    double* ENERGY, int* NDOFEL, int* NRHS, int* NSVARS, double* PROPS,
    int* NPROPS, double* COORDS, int* MCRD, int* NNODE, double* U, double* DU,
    double* V, double* A, int* JTYPE, double* TIME, double* DTIME, int* KSTEP,
    int* KINC,int* JELEM, double* PARAMS, int* NDLOAD,int* JDLTYP, double* ADLMAG,
    double* PREDEF, int* NPREDF, int* LFLAGS,int* MLVARX,double* DDLMAG, int* MDLOAD,
    double* PNEWDT, int* JPROPS, int* NJPROPS, double* PERIOD){

    // Dimension information
    const int nn=*NNODE;
    int nsize=NDOF*nn;
    // Input properties
    const int nparamsI = *NJPROPS;
    const int nparamsR = *NPROPS;
    int paramsI[nparamsI]; // vailid in CPP not C
    double paramsR[nparamsR];
    // coords, displacement, stiffness and force vector
    double xc[nn],yc[nn],zc[nn];
    double disp[nsize];
    double stf[nsize*nsize];
    double qf[nsize];
    double dt = *DTIME;
    // counter
    int i,j;

    // Extract the displacement and coordinates
    for(i=0; i<nsize; i++) disp[i]=U[i];
    for(i=0; i<nn; i++){
        xc[i]=COORDS[0+i*3]; yc[i]=COORDS[1+i*3]; zc[i]=COORDS[2+i*3];
    }
    // Extract element properties
    for(i=0; i<nparamsI; i++) paramsI[i]=JPROPS[i];
    for(i=0; i<nparamsR; i++) paramsR[i]=PROPS[i];

    // Call the integrate function in CohesiveElement Class
    Mesh::getInstance().getElement(*JELEM)->integrate(*KINC, *JELEM, nn,
        xc, yc, zc, disp, paramsR, stf, qf, dt);

    // update the K matrix
    for(i=0; i<nsize*nsize; i++) AMATRX[i]=stf[i];
    // update the RHS
    for(i=0; i<nsize; i++) RHS[i]=-qf[i];

}
