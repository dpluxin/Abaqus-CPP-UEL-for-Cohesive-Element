// CohesiveMaterial.inc
#pragma once

// declare local functions
void deemat(double para[], double damage, bool crack_closed, double dmatrx[NSTR][NSTR]);
void failure_check(double para[], double u[], double s[],
    std::shared_ptr<SolutionVariable> sdv, bool& new_dam);
void cohesive_law(double u[], double DMAX, std::shared_ptr<SolutionVariable> sdv);

// ******** Template Functions of CohesiveMaterial Class ******** //

void CohesiveMaterial::DDSDDE(double u1, double u2, double u3, double dmatrx[NDOF][NDOF],
    double s[], double para[], std::shared_ptr<SolutionVariable> sdv, double dt){

    // local variables
    double u[NSTR];
    double conv_damage, NETA, DMAX;
    bool new_dam, is_closed_crack;
    const double DT_CRITICAL = 0.000001;

    u[0] = u1; u[1] = u2; u[2] = u3;
    NETA        = para[6];
    DMAX        = para[8];
    conv_damage = sdv->damage;
    new_dam     = false;
    is_closed_crack   = u[0] < 0.0;

    // ****** 1. iterative tractions ****** //
    // get the D matrix based on converged damage variable
    deemat(para, conv_damage, is_closed_crack, dmatrx);
    // calculate the traction
    s[0]=dmatrx[0][0]*u[0]+dmatrx[0][1]*u[1]+dmatrx[0][2]*u[2];
    s[1]=dmatrx[1][0]*u[0]+dmatrx[1][1]*u[1]+dmatrx[1][2]*u[2];
    s[2]=dmatrx[2][0]*u[0]+dmatrx[2][1]*u[1]+dmatrx[2][2]*u[2];

    // ****** 2. check the failure onset ****** //
    // determine whether the failure check should be called
    if(sdv->currentStatus == SolutionVariable::Intact){
        failure_check(para, u, s, sdv, new_dam);
    }

    // ****** 3. Check elem status ****** //
    // check the damage status after failure check
    if(sdv->currentStatus == SolutionVariable::Intact){
        // no need to update the stiffness and stress matrix
        sdv->separ = {u1, u2, u3};
        sdv->tract = {s[0], s[1], s[2]};
        return;
    }

    // ****** 4. Cohesive Law ****** //
    // material is damaged, check cohesive law to degrade material stiffness.
    // NO need for new_dam material -> damage determined at failure onset
    if(sdv->currentStatus == SolutionVariable::DamageInitiation && !new_dam){
        // use cohesive to update damage variable and consider viscosity
        cohesive_law(u, DMAX, sdv);
        sdv->damage = NETA/(NETA+dt) * conv_damage + dt/(NETA+dt) * sdv->damage;
    }
    if(sdv->currentStatus == SolutionVariable::FinalFailure) sdv->damage = DMAX;

    // ****** 5. update stiffness and stree ****** //
    // - explicit damage marching is applied when dt is too small:
    //   Dmat and s are caculated based on converged damage variable
    //   but the current damage variable is determined by current u
    if(dt < DT_CRITICAL){
        // Explicit: no need to update d matrix and s vector
    } else{
        // dt is larger than critical, implicit updating
        deemat(para, sdv->damage, is_closed_crack, dmatrx);
        s[0]=dmatrx[0][0]*u[0]+dmatrx[0][1]*u[1]+dmatrx[0][2]*u[2];
        s[1]=dmatrx[1][0]*u[0]+dmatrx[1][1]*u[1]+dmatrx[1][2]*u[2];
        s[2]=dmatrx[2][0]*u[0]+dmatrx[2][1]*u[1]+dmatrx[2][2]*u[2];
    }

    // update other solution variable
    sdv->separ = {u1, u2, u3};
    sdv->tract = {s[0], s[1], s[2]};

}

//----- local functions for material model

// ********************************************************* //
void deemat(double para[], double damage, bool crack_closed, double dmatrx[NSTR][NSTR]){
    // local variables
    double Dnn, Dtt, Dll, DMAX;
    const double RESIDUAL_MODULUS = 0.00001;
    Dnn = para[0]; Dtt = para[0]; Dll = para[0]; DMAX = para[8];


    // calculate degraded stiffness
    if(!crack_closed) Dnn = Dnn * (1.0 - damage);
    Dtt = Dtt * (1.0 - damage); Dll = Dll * (1.0 - damage);
    // consider the maximum damage and small residual stiffness
    if ((damage >= (DMAX-SMALLNUM)) && !crack_closed) Dnn = RESIDUAL_MODULUS;
    if ( damage >= (DMAX-SMALLNUM)) { Dtt=RESIDUAL_MODULUS; Dll=RESIDUAL_MODULUS;}

    // update the D matrxi
    dmatrx[0][0] = Dnn; dmatrx[0][1] = 0.0; dmatrx[0][2] = 0.0;
    dmatrx[1][0] = 0.0; dmatrx[1][1] = Dtt; dmatrx[1][2] = 0.0;
    dmatrx[2][0] = 0.0; dmatrx[2][1] = 0.0; dmatrx[2][2] = Dll;
}

// ********************************************************* //
void failure_check(double para[], double u[], double s[],
    std::shared_ptr<SolutionVariable> sdv, bool& new_dam){
    // local variable
    double t_n, u_n, t_eff, u_eff;
    double t_nc, t_sc, tau_se;
    double Gnc, Gsc, alpha, cof, DMAX, mix_flag;
    double Gn, Gt, Gl, Gmc, lambda;
    // extract from input properties
    t_nc = para[1]; t_sc = para[2];
    Gnc  = para[3]; Gsc  = para[4];
    alpha= para[5]; cof  = para[7];
    DMAX = para[8]; mix_flag = para[9];
    // consider the compressive status
    t_n = std::max(s[0],0.0);  u_n = std::max(u[0],0.0);

    // enhance the shear strength and thoughness under compression
    if(cof > 0.0 && s[0] < 0.0){
        t_sc = t_sc + cof*(-s[0]);
        Gsc  = Gsc*(t_sc/para[2])*(t_sc/para[2]);
    }

    // calculate the failure index for tensile failure
    sdv->failIndex = std::sqrt( (t_n/t_nc)*(t_n/t_nc) + (s[1]/t_sc)*(s[1]/t_sc) +
                               (s[2]/t_sc)*(s[2]/t_sc) );
    // check to see if the failure onset criterion is reached;
    if(sdv->failIndex >= 1.0){
        // update sdv and new damage flag
        sdv->currentStatus = SolutionVariable::DamageInitiation;
        new_dam = true;
        // calculate the effective separation and traction at failure onset
        // - energy components
        Gn = 0.5 * t_n * u_n; Gt = 0.5 * s[1] * u[1]; Gl = 0.5 * s[2] * u[2];
        // - the mode mixity
        lambda = (Gt + Gl) / (Gn + Gt + Gl);
        // - mix_mode thoughness
        Gmc = (mix_flag > 0.0) ? (Gnc + (Gsc-Gnc)*pow(lambda,alpha))
            : (pow(pow((1.0-lambda)/Gnc,alpha)+pow(lambda/Gsc,alpha),-1.0/alpha));
        // - effective separation
        u_eff = std::sqrt( u_n*u_n + u[1]*u[1] + u[2]*u[2] );
        // effective traction is calculated such that
        // 0.5 * u_eff * t_eff = total strain energy density
        t_eff = 2.0 * (Gn + Gt + Gl) / u_eff;

        // adjust effective separation and traction at failure onset,
        // - adjusted with overshoot
        u_eff = u_eff / sdv->failIndex; t_eff = t_eff / sdv->failIndex;
        // - set up the failure onset separation
        sdv->u0 = u_eff; sdv->t0 = t_eff; sdv->uf = 2.0 * Gmc / t_eff;

        // update damage WITHOUT cohesive law
        // - brittle failure
        if(sdv->uf <= sdv->u0){
            sdv->damage = DMAX; sdv->currentStatus = SolutionVariable::FinalFailure;
        } else{
            sdv->damage = 1.0 - 1.0 / sdv->failIndex;
        }

    } // end check failure onset

}

// ********************************************************* //
void cohesive_law(double u[], double DMAX, std::shared_ptr<SolutionVariable> sdv){
    // local variables
    double u_n, u_eff, dm_tmp;

    // check correctness of cohesive parameter
    if(sdv->u0 <= SMALLNUM || sdv->u0 <= SMALLNUM){
        std::cerr << "zero cohesive parameters, CohesiveMaterial.inc" << std::endl;
        std::exit(1);
    }

    // - brittle failure
    if(sdv->uf <= sdv->u0 + SMALLNUM){
        sdv->damage = DMAX; sdv->currentStatus = SolutionVariable::FinalFailure;
    }

    // - damage progress
    u_n   = std::max(u[0],0.0);
    u_eff = std::sqrt(u_n*u_n + u[1]*u[1] + u[2]*u[2]);
    if(u_eff > SMALLNUM){
        // linear softening law to calculate the damage factor
        dm_tmp = (sdv->uf / u_eff) * (u_eff-sdv->u0) / (sdv->uf-sdv->u0);
        // update damage variable if it becomes larger
        if(dm_tmp > sdv->damage) sdv->damage = dm_tmp;
        // check final failure
        if(sdv->damage > DMAX - SMALLNUM) {
            sdv->damage = DMAX; sdv->currentStatus = SolutionVariable::FinalFailure;
        }
    } // - end damage progress

}
