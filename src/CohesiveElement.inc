// CohesiveElement.inc
#pragma once

// ******** Template Functions of CohesiveElement Class ******** //
/*
        All the functions define here have the full access to class
    variables. Not that individual instance (pointer to instance) is
    used to call the functions below (this exists)
*/

void CohesiveElement::integrate(int kinc,int elemID,int nn,double xc[],double yc[],
    double zc[],double disp[],double paramsR[],double stf[],double qf[],double dt){

    /*  Input variables: the passed-in array may lose the dimension information
        --int nsize = NDOF*nn;
        xc[nn], yc[nn], zc[nn], disp[nsize], paramsR[nparamsR]
        stf[nsize*nsize], qf[nsize]
    */

    // local variables
    const int NNODE6 = 6;
    bool last_converged = false; // flag for converged INC, true for 1st new INC

    /* -- analysis control logical variable
        check if last iteration has converged by checking if the global time has advanced;
        global time != local time, last iteration is converged and this is the 1st iteration
        of the new increment since local time will be synchronized -- */
    if(UelVariables::globalTime != localTime ) {
        last_converged = true; // set the converge flag
        setLocalTime(UelVariables::globalTime);
    }

    /* -- update sdvs based on convergence check: deep copy is needed
        update converged sdv with iterating sdv when last iteration is converged
        and revert iterating sdv back to the last converged sdv if otherwise
    */
    if(last_converged){
        // Copy the SolutionVariable objects in itrSdvs to conSdvs
        // - Copy the content of the SolutionVariable
        for (int i = 0; i < itrSdvs.size(); i++) *conSdvs[i] = *itrSdvs[i];
    }else{
        // Copy the SolutionVariable objects in conSdvs to itrSdvs
        // - Copy the content of the SolutionVariable
        for (int i = 0; i < conSdvs.size(); i++) *itrSdvs[i] = *conSdvs[i];
    }

    // update node disp
    for(int i = 0; i < nn; i ++){
        nodes[i]->setDisplacements({disp[i*NDOF], disp[i*NDOF+1], disp[i*NDOF+2]});
    }

    // integrate the cohesive element based on
    if(nn == NNODE6){
        kq_coh3d6(xc,yc,zc,disp,paramsR,stf,qf,dt);
    }else{
        kq_coh3d8(xc,yc,zc,disp,paramsR,stf,qf,dt);
    }

}

// integration of kq_coh3d6 element
void CohesiveElement::kq_coh3d6(double xc[],  double yc[],  double zc[],
    double disp[], double paramsR[], double stf[], double qf[], double dt){
    // ----------- USE GLOBAL VARIABLES -----------
    // NNODE6: number of nodes
    // NDOF: number of DoF
    // NINT: number of IGPs
    // --------------------------------------------
    const int NNODE6 = 6;
    const int nsize = NNODE6*NDOF;
    const int NINT = 3;
    int  i,j,k,num;
    double wg;
    double nmat[NDOF][nsize];
    int stf_size;
    int ip;
    double x1,x2;
    double l1[]={1,0,0};
    double l2[]={0,1,0};
    double e1[NDOF],e2[NDOF],e3[NDOF];
    double xM[NNODE6/2],yM[NNODE6/2],zM[NNODE6/2],g1[NDOF],g2[NDOF],xsum;
    double h[NNODE6/2];
    double area;
    /** cohesive parameter **/
    double u1,u2,u3;
    double p1,p2,p3;
    double dp1d1,dp1d2,dp1d3;
    double dp2d1,dp2d2,dp2d3;
    double dp3d1,dp3d2,dp3d3;
    double dmat[NDOF][NDOF],s[NDOF];

    // size of stiffness matrix
    stf_size=nsize*nsize;
    // initilize the stiffness and nodal force
    for(i=0;i<stf_size;i++) stf[i]=0.0;
    for( i=0;i< nsize;i++) qf[i]=0.0;

    // local base vectors
    for(i=0;i<NNODE6/2;i++){
        xM[i]=0.5*(xc[i]+xc[i+3]);
        yM[i]=0.5*(yc[i]+yc[i+3]);
        zM[i]=0.5*(zc[i]+zc[i+3]);
    }
    // G1
    g1[0]=xM[1]-xM[0]; g1[1]=yM[1]-yM[0]; g1[2]=zM[1]-zM[0];
    // G2
    g2[0]=xM[2]-xM[0]; g2[1]=yM[2]-yM[0]; g2[2]=zM[2]-zM[0];
    // Normal dirction: e1 = (g1xg2)
    cross_product(g1, g2, e1);
    area= vector_magnitude(e1) *0.5;
    normalize(e1, e1);  // e1 = (g1xg2)/|g1xg2|
    // e2 = g1/||g1||
    normalize(g1, e2);
    // e3 = (e1xe2)/||e1xe2||
    cross_product(e1, e2, e3);
    normalize(e3, e3);

    // set CZM model
    // loop for integration
    for(ip=0; ip<NINT; ip++){
        // shape function
        x1=l1[ip]; x2=l2[ip];
        h[0]=x1; h[1]=x2; h[2]=1.0-x1-x2;

        // Set N matrix
        for (int i = 0; i < 3; ++i) { // Iterate over the tangential directions (E1, E2, E3)
            const double* e = (i == 0) ? e1 : (i == 1) ? e2 : e3; // Select the appropriate e vector (e1, e2, e3)
            for (int j = 0; j < 3; ++j) { // Iterate over the components (1, 2, 3)
                for (int k = 0; k < 3; ++k) { // Iterate over h (0, 1, 2)
                    nmat[i][j + 3 * k] = -h[k] * e[j];        // Fill the negative components
                    nmat[i][j + 3 * k + 9] = h[k] * e[j];    // Fill the positive components
                }
            }
        }

        /*** Weight ***/
        wg=area/3.0;
        /*** u1,u2.u3 ***/
        u1=0.0; u2=0.0; u3=0.0;
        for(j=0;j<nsize;j++){
            u1+=nmat[0][j]* disp[j]; u2+=nmat[1][j]* disp[j]; u3+=nmat[2][j]* disp[j];
        }

        // ========================================================
        // call CohesiveMaterial::DDSDDE to material stiffness
        // --- iterating solution variable: passing as a pointer
        //     std::shared_ptr<SolutionVariable> itrSdvs[ip]
        // ========================================================
        CohesiveMaterial::DDSDDE(u1,u2,u3,dmat,s,paramsR,itrSdvs[ip],dt);

        // calculate the stiffness matrix and force vector
        dp1d1=dmat[0][0]; dp1d2=dmat[0][1]; dp1d3=dmat[0][2];
        dp2d1=dmat[1][0]; dp2d2=dmat[1][1]; dp2d3=dmat[1][2];
        dp3d1=dmat[2][0]; dp3d2=dmat[2][1]; dp3d3=dmat[2][2];
        p1   =s[0];       p2   =s[1];       p3   =s[2];

        /*** Stiffness matrix ***/
        for( j=0;j< nsize;j++){
            for( k=0;k<nsize;k++){
                stf[k+nsize*j]+=
                ( dp1d1*nmat[0][k]*nmat[0][j]
                + dp1d2*nmat[0][k]*nmat[1][j]
                + dp1d3*nmat[0][k]*nmat[2][j]
                + dp2d1*nmat[1][k]*nmat[0][j]
                + dp2d2*nmat[1][k]*nmat[1][j]
                + dp2d3*nmat[1][k]*nmat[2][j]
                + dp3d1*nmat[2][k]*nmat[0][j]
                + dp3d2*nmat[2][k]*nmat[1][j]
                + dp3d3*nmat[2][k]*nmat[2][j])*wg;
            }
        }
        /** Add Q force ***/
        for(j=0;j<nsize;j++){
            qf[j]+=(p1*nmat[0][j]+p2*nmat[1][j]+p3*nmat[2][j])*wg;
        }

    }// end of IG iteration

} // end of kq_coh3d6 function


// integration of kq_coh3d8 element
void CohesiveElement::kq_coh3d8(double xc[],  double yc[],  double zc[],
    double disp[], double paramsR[], double stf[], double qf[], double dt){
    // ----------- USE GLOBAL VARIABLES -----------
    // NNODE8: number of nodes
    // NDOF: number of DoF
    // NGAUSS: number of IGPs along each direction
    // --------------------------------------------
    // define the local variables
    const int NNODE8 = 8;
    const int nsize  = NNODE8*NDOF;
    const int NGAUSS = 2;
    int  i,j,k;
    double wg;
    double nmat[NDOF][nsize];
    int stf_size;
    int lx,ly;
    double ri,si;
    double xg[NGAUSS],wgt[NGAUSS];
    double e1[NDOF],e2[NDOF],e3[NDOF];
    double xM[NNODE8/2],yM[NNODE8/2],zM[NNODE8/2],g1[NDOF],g2[NDOF];
    double h[NNODE8/2],hr[NNODE8/2],hs[NNODE8/2],det;
    double rp,sp,rm,sm;
    /** cohesive parameter **/
    double u1,u2,u3;
    double p1,p2,p3;
    double dp1d1,dp1d2,dp1d3;
    double dp2d1,dp2d2,dp2d3;
    double dp3d1,dp3d2,dp3d3;
    double dmat[NDOF][NDOF],s[NDOF];

    // Integration Table: natural coords and weight
    xg[0] = -0.577350269189626; wgt[0] = 1.0;
    xg[1] =  0.577350269189626; wgt[1] = 1.0;
    // size of stiffness matrix
    stf_size = nsize * nsize;
    // initilize the stiffness and nodal force
    for(i=0;i<stf_size;i++) stf[i]=0.0;
    for(i=0;i< nsize;i++) qf[i]=0.0;

    // set CZM model
    /** loop for Gaussian integration **/
    for(lx=0; lx < NGAUSS; lx++){
        ri=xg[lx];
        for(ly=0; ly < NGAUSS; ly++){
            si=xg[ly];
            // --- it is the [lx*2+ly] integration point
            // shape functions
            rp   = 1.0+ri;     sp   = 1.0+si;     rm   = 1.0-ri;     sm   = 1.0-si;
            h[0] = 0.25*rm*sm; h[1] = 0.25*rp*sm; h[2] = 0.25*rp*sp; h[3] = 0.25*rm*sp;
            // derivatives of shape functions
            hr[0]=-0.25*sm;    hr[1]= 0.25*sm;    hr[2]= 0.25*sp;    hr[3]=-0.25*sp;
            hs[0]=-0.25*rm;    hs[1]=-0.25*rp;    hs[2]= 0.25*rp;    hs[3]= 0.25*rm;
            // Mid surface coords
            for(i=0;i<4;i++){
                xM[i]=0.5*(xc[i]+xc[i+4]); yM[i]=0.5*(yc[i]+yc[i+4]); zM[i]=0.5*(zc[i]+zc[i+4]);
            }

            // zero base vectors
            g1[0]=0.0; g1[1]=0.0; g1[2]=0.0;
            g2[0]=0.0; g2[1]=0.0; g2[2]=0.0;
            for(i=0;i<4;i++){
                // G1 1st tangent vector of in-plane direction: dNi/dr * Xi
                g1[0]+=hr[i]*xM[i]; g1[1]+=hr[i]*yM[i]; g1[2]+=hr[i]*zM[i];
                // G1 2nd tangent vector of in-plane direction: dNi/ds * Xi
                g2[0]+=hs[i]*xM[i]; g2[1]+=hs[i]*yM[i]; g2[2]+=hs[i]*zM[i];
            }
            // -- Normal dirction: e1 = (g1xg2)
            cross_product(g1, g2, e1);
            det = vector_magnitude(e1);  // det = |E1|
            normalize(e1, e1);  // e1 = (g1xg2)/|g1xg2|
            // -- e2 = g1/||g1||
            normalize(g1, e2);  // E2 = g1 / |g1|
            // -- e3 = (e1xe2)/||e1xe2||
            cross_product(e1, e2, e3);
            normalize(e3, e3);

            /*** Set N matrix ***/
            for (int dir = 0; dir < 3; ++dir) {
                // Loop over the 3 components (x, y, z)
                for (int i = 0; i < 3; ++i) {
                    // Loop over the 4 terms of h (h[0], h[1], h[2], h[3])
                    for (int j = 0; j < 4; ++j) {
                        // Calculate the matrix indices for the first part of nmat (negative values)
                        nmat[dir][3*j + i] = -h[j] * (dir == 0 ? e1[i] : (dir == 1 ? e2[i] : e3[i]));
                        // Calculate the matrix indices for the second part of nmat (positive values)
                        nmat[dir][3*j + i + 12] = h[j] * (dir == 0 ? e1[i] : (dir == 1 ? e2[i] : e3[i]));
                    }
                }
            }

            /*** Weight ***/
            wg=wgt[lx]*wgt[ly]*det;
            /*** u1,u2,u3 ***/
            u1=0.0; u2=0.0; u3=0.0;
            for(j=0;j<nsize;j++){
                u1+=nmat[0][j]* disp[j]; u2+=nmat[1][j]* disp[j]; u3+=nmat[2][j]* disp[j];
            }

            // ========================================================
            // call CohesiveMaterial::DDSDDE to material stiffness
            // --- iterating solution variable: passing as a pointer
            //     std::shared_ptr<SolutionVariable> itrSdvs[lx*2+ly]
            // ========================================================
            CohesiveMaterial::DDSDDE(u1,u2,u3,dmat,s,paramsR,itrSdvs[lx*2+ly],dt);

            // calculate the stiffness matrix and force vector
            dp1d1=dmat[0][0]; dp1d2=dmat[0][1]; dp1d3=dmat[0][2];
            dp2d1=dmat[1][0]; dp2d2=dmat[1][1]; dp2d3=dmat[1][2];
            dp3d1=dmat[2][0]; dp3d2=dmat[2][1]; dp3d3=dmat[2][2];
            p1   =s[0];       p2   =s[1];       p3   =s[2];
            // Stiffness matrix
            for( j=0;j< nsize;j++){
                for( k=0;k<nsize;k++){
                    stf[k+nsize*j]+=
                    ( dp1d1*nmat[0][k]*nmat[0][j]
                    + dp1d2*nmat[0][k]*nmat[1][j]
                    + dp1d3*nmat[0][k]*nmat[2][j]
                    + dp2d1*nmat[1][k]*nmat[0][j]
                    + dp2d2*nmat[1][k]*nmat[1][j]
                    + dp2d3*nmat[1][k]*nmat[2][j]
                    + dp3d1*nmat[2][k]*nmat[0][j]
                    + dp3d2*nmat[2][k]*nmat[1][j]
                    + dp3d3*nmat[2][k]*nmat[2][j])*wg;
                }
            }
            // Add Q force
            for(j=0;j<nsize;j++){
                qf[j]+=(p1*nmat[0][j]+p2*nmat[1][j]+p3*nmat[2][j])*wg;
            }

        } // end of inner IG iteration
    } // end of outer IG iteration

} // end of kq_coh3d8 function
