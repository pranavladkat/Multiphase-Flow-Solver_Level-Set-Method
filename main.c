static char help[] = "Solves 2D multiphase flow problem using level set formulation.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <engine.h>


typedef struct{
    Vec U, V, P, rho, mu;
}Solution;

typedef struct{

    Solution soln;
    Solution solstar;
    Solution solnp1;
    Solution solnm1;
    Vec RHSu, RHSv;

    PetscReal h, dt, dt_old, Re, Bond, time;
    PetscReal smallval;

    Vec RHS;
    Mat Jac;
    PetscInt counter, M, N;

    Vec phi_0, phi, phi_new, H;
    Vec normalX, normalY;
    PetscReal Volume0, VolumeNew, Length, dphi;

    double rho_ratio, mu_ratio;
    double gravity, sigma, R;

}UserContext;



/* User defined functions*/
extern PetscErrorCode CreateStructures(DM, UserContext *);
extern PetscErrorCode DestroyStructures(UserContext *);
extern PetscErrorCode PredictVelocity(DM, UserContext *);
extern PetscErrorCode WriteVec(Vec, char const *);
extern PetscErrorCode ComputeExactSol(DM, UserContext *);
extern PetscErrorCode ComputeRHS(DM, UserContext *);
extern PetscErrorCode ComputeMatrix(DM, UserContext *);
extern PetscErrorCode WriteMat(Mat,char const *);
extern PetscErrorCode Projection_Step(KSP, UserContext *);
extern PetscErrorCode UpdateVelocity(DM, UserContext *);
extern PetscErrorCode DefineLevelSet(DM, UserContext *);
extern PetscErrorCode Reinitialize(DM, UserContext *);
extern PetscReal sign(PetscReal);
extern PetscErrorCode DefineVariables(DM, UserContext *);
extern PetscErrorCode AdvectInterface(DM, UserContext *);
extern PetscErrorCode WriteOutput(UserContext *);
extern PetscErrorCode ComputeTimeStep(UserContext *);
extern PetscErrorCode ConserveMass(DM, UserContext*);
extern PetscErrorCode Screenshot(UserContext *, Engine *);



#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args){

    KSP ksp;
    DM da;
    UserContext user;
    PetscErrorCode ierr;

    /* start matlab engine */
    //Engine *ep = engOpen(NULL);

    /* Define parameters */
    user.R = 1.0;      /* Radius of bubble */
    user.M = 60; user.N = user.M;
    user.h = 7*user.R/(user.M-1);
    user.rho_ratio = 1.0/40;
    user.mu_ratio = 1.0/40;
    user.Re = 100.0;
    user.Bond = 200.0;
    user.smallval = user.h*user.h;
    int itmax = 2000;
    user.counter = 0;
    user.time = 0.0;

    ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                        user.M,user.M,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da); CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
    ierr = CreateStructures(da,&user); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"Reynolds no.= %g Bond no.= %g\n",user.Re,user.Bond);

    ierr = DefineLevelSet(da,&user); CHKERRQ(ierr);
    ierr = Reinitialize(da,&user); CHKERRQ(ierr);

    for(user.counter = 0; user.counter < itmax; user.counter++){

        PetscPrintf(PETSC_COMM_WORLD,"\n\nIteration: %d\n",user.counter+1);

        ierr = ComputeTimeStep(&user); CHKERRQ(ierr);
        ierr = AdvectInterface(da,&user); CHKERRQ(ierr);
        ierr = Reinitialize(da,&user); CHKERRQ(ierr);
        ierr = DefineVariables(da,&user); CHKERRQ(ierr);
        ierr = ConserveMass(da,&user); CHKERRQ(ierr);
        ierr = PredictVelocity(da,&user); CHKERRQ(ierr);
        ierr = ComputeRHS(da,&user); CHKERRQ(ierr);
        ierr = ComputeMatrix(da,&user); CHKERRQ(ierr);
        ierr = Projection_Step(ksp,&user); CHKERRQ(ierr);
        ierr = UpdateVelocity(da,&user); CHKERRQ(ierr);
        //ierr = Screenshot(&user,ep); CHKERRQ(ierr);

        if(user.time >= 4.5){
            break;
        }

    }

    ierr = WriteOutput(&user); CHKERRQ(ierr);

    ierr = DestroyStructures(&user); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    //ierr = DMDestroy(&da); CHKERRQ(ierr);

    ierr = PetscFinalize();



    return(0);
}




/* Create vectors */
#undef __FUNCT__
#define __FUNCT__ "CreateStructures"
PetscErrorCode CreateStructures(DM da, UserContext *user){

    PetscErrorCode ierr;

    ierr = DMCreateGlobalVector(da,&user->soln.U); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->soln.V); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->soln.P); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->soln.rho); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->soln.mu); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solstar.U); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solstar.V); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solstar.P); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnp1.U); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnp1.V); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnp1.P); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnp1.rho); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnp1.mu); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->RHS); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnm1.U); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->solnm1.V); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->phi); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->phi_0); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->phi_new); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->H); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->normalX); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->normalY); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->RHSu); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&user->RHSv); CHKERRQ(ierr);

    ierr = DMCreateMatrix(da,&user->Jac); CHKERRQ(ierr);
    return(0);
}



/* Destroy vectors */
#undef __FUNCT__
#define __FUNCT__ "DestroyStructures"
PetscErrorCode DestroyStructures(UserContext *user){

    PetscErrorCode ierr;

    ierr = VecDestroy(&user->soln.U); CHKERRQ(ierr);
    ierr = VecDestroy(&user->soln.V); CHKERRQ(ierr);
    ierr = VecDestroy(&user->soln.P); CHKERRQ(ierr);
    ierr = VecDestroy(&user->soln.rho); CHKERRQ(ierr);
    ierr = VecDestroy(&user->soln.mu); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solstar.U); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solstar.V); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solstar.P); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnp1.U); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnp1.V); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnp1.P); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnp1.rho); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnp1.mu); CHKERRQ(ierr);
    ierr = VecDestroy(&user->RHS); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnm1.U); CHKERRQ(ierr);
    ierr = VecDestroy(&user->solnm1.V); CHKERRQ(ierr);
    ierr = VecDestroy(&user->phi); CHKERRQ(ierr);
    ierr = VecDestroy(&user->phi_0); CHKERRQ(ierr);
    ierr = VecDestroy(&user->phi_new); CHKERRQ(ierr);
    ierr = VecDestroy(&user->H); CHKERRQ(ierr);
    ierr = VecDestroy(&user->normalX); CHKERRQ(ierr);
    ierr = VecDestroy(&user->normalY); CHKERRQ(ierr);
    ierr = VecDestroy(&user->RHSu); CHKERRQ(ierr);
    ierr = VecDestroy(&user->RHSv); CHKERRQ(ierr);

    ierr = MatDestroy(&user->Jac); CHKERRQ(ierr);
    return(0);
}




#undef __FUNCT__
#define __FUNCT__ "WriteVec"
PetscErrorCode WriteVec(Vec vec, char const *name){
    PetscViewer viewer;
    PetscErrorCode ierr;

    /* Create "Output" directory */
    struct stat st = {0};
    if(stat("Output",&st) == -1)
        mkdir("Output",0777);

    char filename[64] = "Output/Vec_"; char pfix[12] = ".m";
    strcat(filename,name); strcat(filename,pfix);
    ierr = PetscObjectSetName((PetscObject)vec,name); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
    ierr = VecView(vec,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "PredictVelocity"
extern PetscErrorCode PredictVelocity(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscReal **_Un, **_Vn, **_Pn, **_Ustr, **_Vstr, **_Unm1, **_Vnm1, **_Pstr, **_mu, **_rho;
    PetscReal **_N1, **_N2, **_phi, **_H, **_RHSu, **_RHSv;
    PetscReal h,dt, dt_old;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal ux, uy, vx, vy;
    PetscReal convect, viscous, curvature, delta, surftension;
    double pi = 3.1415926535, eps = 1.5*user->h;
    PetscReal GV1,GV2;
    PetscReal phi_x, phi_y, phi_xx, phi_yy, phi_xy;
    PetscReal mu1, mu2,mu3,mu4;

    h = user->h; dt = user->dt; dt_old = user->dt_old;

    ierr = DMDAVecGetArray(da,user->soln.U,&_Un); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.V,&_Vn); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.P,&_Pn); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.U,&_Ustr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.V,&_Vstr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solnm1.U,&_Unm1); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solnm1.V,&_Vnm1); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.P,&_Pstr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.mu,&_mu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->normalX,&_N1); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->normalY,&_N2); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->H,&_H); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->RHSu,&_RHSu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->RHSv,&_RHSv); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++){
            if(i > 0 && i < user->M-1 && j > 0 && j < user->M-1){


                /* convective term */
                if(user->counter == 0){
                    /* central difference */
                    ux = 0.5*(_Un[j][i+1] - _Un[j][i-1])/h;
                    uy = 0.5*(_Un[j+1][i] - _Un[j-1][i])/h;

                    convect = _Un[j][i]*ux + _Vn[j][i]*uy;

                }else if(user->counter > 0){
                    /* Adam Bashford */

                        convect = 0.5*(((dt + 2*dt_old)/dt_old)*(_Un[j][i]*(0.5*(_Un[j][i+1] - _Un[j][i-1])/h) + _Vn[j][i]*(0.5*(_Un[j+1][i] - _Un[j-1][i])/h))
                                -(dt/dt_old)*(_Unm1[j][i]*(0.5*(_Unm1[j][i+1] - _Unm1[j][i-1])/h) + _Vnm1[j][i]*(0.5*(_Unm1[j+1][i] - _Unm1[j-1][i])/h)));

                }

                mu1 = 0.5*(_mu[j][i+1] + _mu[j][i]);
                mu2 = 0.5*(_mu[j][i-1] + _mu[j][i]);
                mu3 = 0.5*(_mu[j+1][i] + _mu[j][i]);
                mu4 = 0.5*(_mu[j-1][i] + _mu[j][i]);

                viscous = (2*mu1*((_Un[j][i+1] - _Un[j][i])/h) - 2*mu2*((_Un[j][i] - _Un[j][i-1])/h))/h
                        + (mu3*((_Un[j+1][i] - _Un[j][i])/h) - mu4*((_Un[j][i] - _Un[j-1][i])/h))/h
                        + (mu3*((_Vn[j+1][i+1] - _Vn[j+1][i-1] + _Vn[j][i+1] - _Vn[j][i-1])/(4*h))
                        - mu4*((_Vn[j][i+1] - _Vn[j][i-1] + _Vn[j-1][i+1] - _Vn[j-1][i-1])/(4*h)))/h;


                viscous = viscous/(user->Re*_rho[j][i]);


                /* surface tension term */
                phi_x = (_phi[j][i+1] - _phi[j][i-1])/(2*h);
                phi_y = (_phi[j+1][i] - _phi[j-1][i])/(2*h);
                phi_xx = (_phi[j][i+1] + _phi[j][i-1] - 2*_phi[j][i])/(h*h);
                phi_yy = (_phi[j+1][i] + _phi[j-1][i] - 2*_phi[j][i])/(h*h);
                phi_xy = (_phi[j+1][i+1] - _phi[j-1][i+1] - _phi[j+1][i-1] + _phi[j-1][i-1])/(4*h*h);

                GV1 = (_phi[j][i+1] - _phi[j][i-1])/(2*h);
                GV2 = (_phi[j+1][i] - _phi[j-1][i])/(2*h);

                curvature = (phi_xx*(phi_y*phi_y) + phi_yy*(phi_x*phi_x) - 2*phi_xy*phi_x*phi_y) / pow((phi_x*phi_x + phi_y*phi_y + user->smallval),1.5);

                if(curvature < -1.0/h){
                    curvature = -1.0/h;
                }else if(curvature > 1.0/h){
                    curvature = 1.0/h;
                }


                if( fabs(_phi[j][i]) <= 1.5*h){
                    delta = 0.5*(1 + cos(pi*_phi[j][i]/eps))/eps;
                }else if(fabs(_phi[j][i]) > 1.5*h){
                    delta = 0.0;
                }

                surftension = curvature*delta*GV1/(user->Bond*_rho[j][i]);


                /* compute U star */

                if(user->counter == 0){
                    /* Backward Euler */
                    _Ustr[j][i] = _Un[j][i] + dt*(-convect + viscous - surftension);
                    _RHSu[j][i] = -convect + viscous - surftension;
                }
                else if(user->counter > 0){
                    /* Adam Bashford */
                    _Ustr[j][i] = _Un[j][i] + dt*(0.5*(((dt + 2*dt_old)/dt_old)*(-convect + viscous - surftension) - (dt/dt_old)*_RHSu[j][i]));
                    _RHSu[j][i] = -convect + viscous - surftension;
                }


                /* convective term */
                if(user->counter == 0){
                    vx = 0.5*(_Vn[j][i+1] - _Vn[j][i-1])/h;
                    vy = 0.5*(_Vn[j+1][i] - _Vn[j-1][i])/h;

                    convect = _Un[j][i]*vx + _Vn[j][i]*vy;
                }
                else if(user->counter > 0){
                        convect = 0.5*(((dt + 2*dt_old)/dt_old)*(_Un[j][i]*(0.5*(_Vn[j][i+1] - _Vn[j][i-1])/h) + _Vn[j][i]*(0.5*(_Vn[j+1][i] - _Vn[j-1][i])/h))
                                -(dt/dt_old)*(_Unm1[j][i]*(0.5*(_Vnm1[j][i+1] - _Vnm1[j][i-1])/h) + _Vnm1[j][i]*(0.5*(_Vnm1[j+1][i] - _Vnm1[j-1][i])/h)));
                }


                viscous = (mu1*((_Vn[j][i+1] - _Vn[j][i])/h) - mu2*((_Vn[j][i] - _Vn[j][i-1])/h))/h
                        + (mu1*((_Un[j+1][i+1] - _Un[j-1][i+1] + _Un[j+1][i] - _Un[j-1][i])/(4*h))
                        - mu2*((_Un[j+1][i] - _Un[j-1][i] + _Un[j+1][i-1] - _Un[j-1][i-1])/(4*h)))/h
                        + (2*mu3*((_Vn[j+1][i] - _Vn[j][i])/h) - 2*mu4*((_Vn[j][i] - _Vn[j-1][i])/h))/h;


                viscous = viscous/(user->Re*_rho[j][i]);

                /* surface tension term */
                surftension = curvature*delta*GV2/(user->Bond*_rho[j][i]);

                /* compute V star */
                if(user->counter == 0){
                    _Vstr[j][i] = _Vn[j][i] + dt*(-convect + viscous - surftension -1);
                    _RHSv[j][i] = -convect + viscous - surftension -1;
                }
                else if(user->counter > 0){
                    /* Adam Bashford */
                    _Vstr[j][i] = _Vn[j][i] + dt*(  0.5*(((dt + 2*dt_old)/dt_old)*(-convect + viscous - surftension -1) - (dt/dt_old)*_RHSv[j][i]));
                    _RHSv[j][i] = -convect + viscous - surftension -1;
                }

            }else if(i == 0 || i == user->M-1 || j == 0 || j == user->M-1){

                /* 4 corners */
                if(i == 0 && j == 0){
                    _Ustr[j][i] = _Un[j][i];
                    _Vstr[j][i] = _Vn[j][i];
                }else if(i == 0 && j == user->M-1){
                    _Ustr[j][i] = _Un[j][i];
                    _Vstr[j][i] = _Vn[j][i];
                }else if(i == user->M-1 && j == 0){
                    _Ustr[j][i] = _Un[j][i];
                    _Vstr[j][i] = _Vn[j][i];
                }else if(i == user->M-1 && j == user->M-1){
                    _Ustr[j][i] = _Un[j][i];
                    _Vstr[j][i] = _Vn[j][i];
                }

                /* normal B.C. */
                if(i == 0 || i == user->M-1){
                    _Ustr[j][i] = _Un[j][i];
                }

                if(j == 0 || j == user->M-1){
                    _Vstr[j][i] = _Vn[j][i];
                }

                /* tangent b.c */
                if(j == 0 || j == user->M-1){
                    if(i != 0 && i != user->M-1){
                        _Ustr[j][i] = _Un[j][i] + dt*(_Pstr[j][i+1] - _Pstr[j][i-1])/(2*h);
                    }
                }

                if(i == 0 && i == user->M-1){
                    if(j != 0 && j != user->M-1){
                        _Vstr[j][i] = _Vn[j][i] + dt*(_Pstr[j+1][i] -  _Pstr[j-1][i])/(2*h);
                    }
                }

            }

        }
    }

    ierr = DMDAVecRestoreArray(da,user->soln.U,&_Un); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.V,&_Vn); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.P,&_Pn); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.U,&_Ustr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.V,&_Vstr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solnm1.U,&_Unm1); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solnm1.V,&_Vnm1); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.P,&_Pstr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.mu,&_mu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->normalX,&_N1); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->normalY,&_N2); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->H,&_H); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->RHSu,&_RHSu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->RHSv,&_RHSv); CHKERRQ(ierr);

    return(0);
}





#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscReal **_Ustr, **_Vstr, **_b;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal dt = user->dt, h = user->h, ux, vy;

    ierr = DMDAVecGetArray(da,user->solstar.U,&_Ustr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.V,&_Vstr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->RHS,&_b); CHKERRQ(ierr);
    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    for(i = xs; i < xs+xm; i++){
        for(j = ys; j < ys+ym; j++){

            if(i == 0){
                ux = 0.5*(-3*_Ustr[j][i] + 4*_Ustr[j][i+1] - _Ustr[j][i+2])/h;
            }else if(i == user->M-1){
                ux = 0.5*(3*_Ustr[j][i] - 4*_Ustr[j][i-1] + _Ustr[j][i-2])/h;
            }else if(i > 0 && i < user->M-1 ){
                ux = (0.5*(_Ustr[j][i+1] - _Ustr[j][i-1]))/h;
                //ux = (0.5*(_Ustr[j][i] + _Ustr[j][i+1]) - 0.5*(_Ustr[j][i] + _Ustr[j][i-1]))/h;
            }

            if(j == 0){
                vy = 0.5*(-3*_Vstr[j][i] + 4*_Vstr[j+1][i] - _Vstr[j+2][i])/h;
            }else if (j == user->M-1){
                vy = 0.5*(3*_Vstr[j][i] - 4*_Vstr[j-1][i] + _Vstr[j-2][i])/h;
            }else if(j > 0 && j < user->M-1){
                vy = (0.5*(_Vstr[j+1][i] - _Vstr[j-1][i]))/h;
                //vy = (0.5*(_Ustr[j][i] + _Ustr[j+1][i]) - 0.5*(_Ustr[j][i] + _Ustr[j-1][i]))/h;
            }

            _b[j][i] = (ux + vy)/dt;

        }
    }

    ierr = VecAssemblyBegin(user->RHS); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(user->RHS); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(da,user->solstar.U,&_Ustr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.V,&_Vstr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->RHS,&_b); CHKERRQ(ierr);

    MatNullSpace nullspace;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullspace,user->RHS); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace); CHKERRQ(ierr);


    return(0);
}



#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal v[5];
    MatStencil row, col[5];
    PetscReal h = user->h, h2 = h*h;
    PetscReal a1, a2, a3, a4, a5;
    PetscReal **_rho;

    ierr = DMDAVecGetArray(da,user->soln.rho,&_rho);
    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    if(user->counter > 0){
        ierr = MatZeroEntries(user->Jac);
    }


    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {

            row.i = i; row.j = j;

            /* compute coefficients */
            if(i == 0){
                a2 = 0;
                a4 = (2.0/(_rho[j][i] + _rho[j][i+1]))/h2;
            }else if(i == user->M-1){
                a2 = (2.0/(_rho[j][i] + _rho[j][i-1]))/h2;
                a4 = 0;
            }else if(i != 0 && i != user->M-1){
                a2 = (2.0/(_rho[j][i] + _rho[j][i-1]))/h2;
                a4 = (2.0/(_rho[j][i] + _rho[j][i+1]))/h2;
            }

            if(j == 0){
                a1 = 0;
                a5 = (2.0/(_rho[j][i] + _rho[j+1][i]))/h2;
            }else if(j == user->N-1){
                a1 = (2.0/(_rho[j][i] + _rho[j-1][i]))/h2;
                a5 = 0;
            }else if(j != 0 && j != user->N-1){
                a1 = (2.0/(_rho[j][i] + _rho[j-1][i]))/h2;
                a5 = (2.0/(_rho[j][i] + _rho[j+1][i]))/h2;
            }

            a3 = -(a1 + a2 + a4 + a5);

            /* assemble matrix */

            if(i == 0 && j == 0){
                v[0] = a3;          col[0].i = i;   col[0].j = j;
                v[1] = a4;          col[1].i = i+1; col[1].j = j;
                v[2] = a5;          col[2].i = i;   col[2].j = j+1;
                ierr = MatSetValuesStencil(user->Jac,1,&row,3,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i == 0 && j != 0 && j != user->N-1){
                v[0] = a1;          col[0].i = i;   col[0].j = j-1;
                v[1] = a3;          col[1].i = i;   col[1].j = j;
                v[2] = a4;          col[2].i = i+1; col[2].j = j;
                v[3] = a5;          col[3].i = i;   col[3].j = j+1;
                ierr = MatSetValuesStencil(user->Jac,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i == 0 && j == user->N-1){
                v[0] = a1;          col[0].i = i;   col[0].j = j-1;
                v[1] = a3;          col[1].i = i;   col[1].j = j;
                v[2] = a4;          col[2].i = i+1; col[2].j = j;
                ierr = MatSetValuesStencil(user->Jac,1,&row,3,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i != 0 && i != user->M-1 && j == 0){
                v[0] = a2;          col[0].i = i-1; col[0].j = j;
                v[1] = a3;          col[1].i = i;   col[1].j = j;
                v[2] = a4;          col[2].i = i+1; col[2].j = j;
                v[3] = a5;          col[3].i = i;   col[3].j = j+1;
                ierr = MatSetValuesStencil(user->Jac,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i == user->M-1 && j == 0){
                v[0] = a2;          col[0].i = i-1; col[0].j = j;
                v[1] = a3;          col[1].i = i;   col[1].j = j;
                v[2] = a5;          col[2].i = i;   col[2].j = j+1;
                ierr = MatSetValuesStencil(user->Jac,1,&row,3,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i == user->M-1 && j != 0 && j != user->N-1){
                v[0] = a1;          col[0].i = i;   col[0].j = j-1;
                v[1] = a2;          col[1].i = i-1; col[1].j = j;
                v[2] = a3;          col[2].i = i;   col[2].j = j;
                v[3] = a5;          col[3].i = i;   col[3].j = j+1;
                ierr = MatSetValuesStencil(user->Jac,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i == user->M-1 && j == user->N-1){
                v[0] = a1;          col[0].i = i;   col[0].j = j-1;
                v[1] = a2;          col[1].i = i-1; col[1].j = j;
                v[2] = a3;          col[2].i = i;   col[2].j = j;
                ierr = MatSetValuesStencil(user->Jac,1,&row,3,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i != 0 && i != user->M-1 && j == user->N-1){
                v[0] = a1;          col[0].i = i;   col[0].j = j-1;
                v[1] = a2;          col[1].i = i-1; col[1].j = j;
                v[2] = a3;          col[2].i = i;   col[2].j = j;
                v[3] = a4;          col[3].i = i+1; col[3].j = j;
                ierr = MatSetValuesStencil(user->Jac,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
            else if(i > 0 && i < user->M-1 && j > 0 && j < user->N-1){
                v[0] = a1;          col[0].i = i;   col[0].j = j-1;
                v[1] = a2;          col[1].i = i-1; col[1].j = j;
                v[2] = a3;          col[2].i = i;   col[2].j = j;
                v[3] = a4;          col[3].i = i+1; col[3].j = j;
                v[4] = a5;          col[4].i = i;   col[4].j = j+1;
                ierr = MatSetValuesStencil(user->Jac,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

        }
    }

    ierr = MatAssemblyBegin(user->Jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(user->Jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = MatSetOption(user->Jac,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);

    MatNullSpace nullspace;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); CHKERRQ(ierr);
    ierr = MatSetNullSpace(user->Jac,nullspace); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace); CHKERRQ(ierr);

    return(0);
}




#undef __FUNCT__
#define __FUNCT__ "WriteMat"
PetscErrorCode WriteMat(Mat mat,char const *name){

    PetscViewer viewer;
    PetscErrorCode ierr;

    /* Create "Output" directory */
    struct stat st = {0};
    if(stat("Output",&st) == -1)
        mkdir("Output",0777);

    char filename[64] = "Output/Mat_"; char pfix[12] = ".m";
    strcat(filename,name); strcat(filename,pfix);
    ierr = PetscObjectSetName((PetscObject)mat,name); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
    ierr = MatView(mat,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "Projection_Step"
PetscErrorCode Projection_Step(KSP ksp, UserContext *user){

    PetscErrorCode ierr;
    MatNullSpace nullsp;
    PetscInt itn;

    ierr = KSPSetOperators(ksp,user->Jac,user->Jac); CHKERRQ(ierr);

    if (user->counter > 0){
        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    }


    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp); CHKERRQ(ierr);
    ierr = KSPSetNullSpace(ksp, nullsp); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.e-9,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
    ierr = KSPSolve(ksp,user->RHS,user->solstar.P); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullsp); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&itn); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"KSP: %d\n",itn); CHKERRQ(ierr);

    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullsp); CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullsp,user->solstar.P); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullsp); CHKERRQ(ierr);

    return(0);

}




#undef __FUNCT__
#define __FUNCT__ "UpdateVelocity"
extern PetscErrorCode UpdateVelocity(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscReal **_Unp1, **_Vnp1, **_Ustr, **_Vstr, **_Pnp1, **_Pstr, **_Pn, **_Un, **_Vn, **_rho;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal dt = user->dt, h = user->h;

    ierr = DMDAVecGetArray(da,user->solnp1.U,&_Unp1); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solnp1.V,&_Vnp1); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solnp1.P,&_Pnp1); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.U,&_Ustr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.V,&_Vstr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->solstar.P,&_Pstr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.P,&_Pn); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.U,&_Un); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.V,&_Vn); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            if(i > 0 && i < user->M-1 && j > 0 && j < user->M-1){

                _Unp1[j][i] = _Ustr[j][i] - dt*(1.0/_rho[j][i])*((_Pstr[j][i+1] - _Pstr[j][i-1])/(2*h));

                _Vnp1[j][i] = _Vstr[j][i] - dt*(1.0/_rho[j][i])*((_Pstr[j+1][i] - _Pstr[j-1][i])/(2*h));

                _Pnp1[j][i] = _Pstr[j][i];
            }
            else if(i ==0 || j == 0 || i == user->M-1 || j == user->M-1){

                _Unp1[j][i] = _Un[j][i];
                _Vnp1[j][i] = _Vn[j][i];
                _Pnp1[j][i] = _Pstr[j][i];
            }

        }
    }

    ierr = DMDAVecRestoreArray(da,user->solnp1.U,&_Unp1); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solnp1.V,&_Vnp1); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solnp1.P,&_Pnp1); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.U,&_Ustr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.V,&_Vstr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->solstar.P,&_Pstr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.P,&_Pn); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.U,&_Un); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.V,&_Vn); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);

    ierr = VecCopy(user->soln.U,user->solnm1.U); CHKERRQ(ierr);
    ierr = VecCopy(user->soln.V,user->solnm1.V); CHKERRQ(ierr);

    ierr = VecCopy(user->solnp1.U,user->soln.U); CHKERRQ(ierr);
    ierr = VecCopy(user->solnp1.V,user->soln.V); CHKERRQ(ierr);
    ierr = VecCopy(user->solnp1.P,user->soln.P); CHKERRQ(ierr);

    return(0);
}




#undef __FUNCT__
#define __FUNCT__ "DefineLevelSet"
PetscErrorCode DefineLevelSet(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscInt i,j,xs,ys,xm,ym;
    DM cda;
    Vec xy;
    DMDACoor2d **_xy;
    PetscReal **_phi, R = user->R;
    PetscReal xmin, xmax, ymin,ymax, xbubble, ybubble;

    xmin = 0.0; xmax = 7*R; ymin = 0.0; ymax = 7*R;
    xbubble = xmin + fabs(xmax-xmin)/2; ybubble = ymin + fabs(xmax-xmin)/4;

    ierr = DMDASetUniformCoordinates(da,xmin,xmax,ymin,ymax,0.0,1.0); CHKERRQ(ierr);

    /* Get coordinates */
    ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(da,&xy); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda,xy,&_xy); CHKERRQ(ierr);

    /* Get Vec arrays */
    ierr = DMDAVecGetArray(da,user->phi_0,&_phi); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    /* Compute exact solution */
    for(i = xs; i < xs+xm; i++){
        for(j = ys; j < ys+ym; j++){

            _phi[j][i] = pow((_xy[j][i].x - xbubble),2) + pow((_xy[j][i].y - ybubble),2) - R*R;

        }
    }

    /* Restore Vec arrays */
    ierr = DMDAVecRestoreArray(da,user->phi_0,&_phi); CHKERRQ(ierr);

    ierr = VecDestroy(&xy); CHKERRQ(ierr);
    ierr = DMDestroy(&cda); CHKERRQ(ierr);

    return(0);
}






#undef __FUNCT__
#define __FUNCT__ "Reinitialize"
PetscErrorCode Reinitialize(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscReal **_phi_0, **_phi, **_phi_new, **_Un, **_Vn;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal s, a, b, c, d, a_plus, a_minus, b_plus, b_minus,c_plus,
            c_minus, d_plus, d_minus, s_plus, s_minus, D, phi_x, phi_y;

    PetscReal phix, phiy, eps = user->h*user->h;

    PetscReal h = user->h, norm;

    int t, itmax = 5*user->M;   		// max iterations;
    PetscReal dt = 0.5*h;      			// time step size(dTau)

    if(user->counter == 0){
        itmax = 5*user->M;   		// max iterations
    }else if(user->counter > 0){
        itmax = 50;   		// max iterations
    }


    ierr = VecCopy(user->phi_0,user->phi); CHKERRQ(ierr);


    ierr = DMDAVecGetArray(da,user->phi_0,&_phi_0); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->phi_new,&_phi_new); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.U,&_Un); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.V,&_Vn); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    for(t = 0; t < itmax; t++){

        for (j=ys; j<ys+ym; j++) {
            for (i=xs; i<xs+xm; i++) {

                /* Compute a, b, c, d */
                if (i == 0){ // B.C.on a using Ghost Node
                    a = (_phi[j][i + 1] - _phi[j][i]) / h;
                }else{
                    a = (_phi[j][i] - _phi[j][i - 1]) / h;
                }

                if (i == user->M - 1){  // B.C.on b using Ghost Node
                    b = (_phi[j][i] - _phi[j][i - 1]) / h;
                }else{
                    b = (_phi[j][i + 1] - _phi[j][i]) / h;
                }

                if (j == 0){ //B.C.on c using Ghost Node
                    c = (_phi[j + 1][i] - _phi[j][i]) / h;
                }else{
                    c = (_phi[j][i] - _phi[j - 1][i]) / h;
                }
                if (j == user->N - 1){	//B.C.on d using Ghost Node
                    d = (_phi[j][i] - _phi[j - 1][i]) / h;
                }else{
                    d = (_phi[j + 1][i] - _phi[j][i]) / h;
                }

                // Compute phi_x
                if (i == 0){
                    phi_x = (_phi[j][i + 1] - _phi[j][i]) / h;
                }
                else if (i == user->N - 1){
                    phi_x = (_phi[j][i] - _phi[j][i - 1]) / h;
                }
                else{ phi_x = (_phi[j][i + 1] - _phi[j][i - 1]) / (2 * h); }

                // Compute phi_y
                if (j == 0){
                    phi_y = (_phi[j + 1][i] - _phi[j][i]) / h;
                }
                else if (j == user->N - 1){
                    phi_y = (_phi[j][i] - _phi[j - 1][i]) / h;
                }
                else{ phi_y = (_phi[j + 1][i] - _phi[j - 1][i]) / (2 * h); }

                /* mollified sign function */
                s = _phi[j][i] / (sqrt(pow(_phi[j][i], 2) + (pow(phi_x, 2) + pow(phi_y, 2))*pow(h, 2)));

                /* Compute positive and negative contribution */
                a_plus = fmax(a, 0.0);
                a_minus = fmin(a, 0.0);
                b_plus = fmax(b, 0.0);
                b_minus = fmin(b, 0.0);
                c_plus = fmax(c, 0.0);
                c_minus = fmin(c, 0.0);
                d_plus = fmax(d, 0.0);
                d_minus = fmin(d, 0.0);
                s_plus = fmax(s, 0.0);
                s_minus = fmin(s, 0.0);


                int flag = 0;
                /* if interface exists -- Use Subcell Fix by Russo and Smereka */
                if(i > 0 && i < user->M-1 && j > 0 && j < user->N-1){
                    if (_phi[j][i]*_phi[j][i+1] < 0 || _phi[j][i]*_phi[j][i-1] < 0 || _phi[j][i]*_phi[j+1][i] < 0 || _phi[j][i]*_phi[j-1][i] < 0){

                        phix = ((_phi_0[j][i] + _phi_0[j][i+1])/2 - (_phi_0[j][i] + _phi_0[j][i-1])/2)/h;
                        phiy = ((_phi_0[j][i] + _phi_0[j+1][i])/2 - (_phi_0[j][i] + _phi_0[j-1][i])/2)/h;

                        D = _phi_0[j][i] / sqrt(pow(phix,2) + pow(phiy,2) + eps);

                        _phi_new[j][i] = _phi[j][i] - (dt / h)*(sign(_phi_0[j][i])*fabs(_phi[j][i]) - D);
                        flag = 1;

                    }
                }

                if(flag == 0){
                    _phi_new[j][i] = _phi[j][i]
                            - dt*(s_plus*(pow((fmax(pow(a_plus, 2), pow(b_minus, 2))) + (fmax(pow(c_plus, 2), pow(d_minus, 2))), 0.5) - 1)
                                  + s_minus*(pow((fmax(pow(a_minus, 2), pow(b_plus, 2))) + (fmax(pow(c_minus, 2), pow(d_plus, 2))), 0.5) - 1));
                }

            }
        }

        ierr = VecAXPY(user->phi,-1.0,user->phi_new); CHKERRQ(ierr);
        ierr = VecNorm(user->phi,NORM_2,&norm);
        if(norm <= 0.000001){
            ierr = VecCopy(user->phi_new,user->phi); CHKERRQ(ierr);
            break;
        }
        ierr = VecCopy(user->phi_new,user->phi); CHKERRQ(ierr);


    }

    PetscPrintf(PETSC_COMM_WORLD,"SDF: %d and norm = %e\n",t,norm);

    ierr = DMDAVecRestoreArray(da,user->phi_0,&_phi_0); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->phi_new,&_phi_new); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.U,&_Un); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.V,&_Vn); CHKERRQ(ierr);

    return(0);
}





#undef __FUNCT__
#define __FUNCT__ "sign"
PetscReal sign(PetscReal sgn){
    if (sgn < 0){ return (-1.0); }
    else if (sgn > 0){ return (1.0); }
    else return (0.0);
}




#undef __FUNCT__
#define __FUNCT__ "DefineVariables"
PetscErrorCode DefineVariables(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal **_H, **_phi, **_Mu, **_rho;
    PetscReal eps = 1.5*user->h;
    double pi = 3.1415926535, rho_ratio = user->rho_ratio, mu_ratio = user->mu_ratio;

    ierr = DMDAVecGetArray(da,user->soln.mu,&_Mu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->H,&_H); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {

            /* compute heaviside */
            if(_phi[j][i] < -eps)
                _H[j][i] = 0;
            else if(fabs(_phi[j][i]) <= eps)
                _H[j][i] = 0.5*(1.0 + _phi[j][i]/eps + sin(pi*_phi[j][i]/eps)/pi);
            else if(_phi[j][i] > eps)
                _H[j][i] = 1;

            /* Define Viscosity */
            _Mu[j][i] = _H[j][i] + (mu_ratio)*(1.0 - _H[j][i]);

            /* Define Density */
            _rho[j][i] = _H[j][i] + (rho_ratio)*(1.0 - _H[j][i]);

        }
    }

    ierr = DMDAVecRestoreArray(da,user->soln.mu,&_Mu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.rho,&_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->H,&_H); CHKERRQ(ierr);

    return(0);
}



#undef __FUNCT__
#define __FUNCT__ "AdvectInterface"
PetscErrorCode AdvectInterface(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscReal **_phi_new, **_phi, **_U, **_V;
    PetscReal phix, phiy;
    PetscReal h = user->h, dt = user->dt;
    PetscInt i,j,xs,ys,xm,ym;

    ierr = DMDAVecGetArray(da,user->phi_new,&_phi_new); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.U,&_U); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user->soln.V,&_V); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {

            if(i == 0){
                phix = (_phi[j][i+1] - _phi[j][i]) / h;
            }else if(i == user->M-1){
                phix = (_phi[j][i] - _phi[j][i - 1]) / h;
            }else if(i !=0 && i != user->M-1){
                if(_U[j][i] > 0){
                    phix = (_phi[j][i] - _phi[j][i-1])/h;
                }else if(_U[j][i] < 0){
                    phix = (_phi[j][i+1] - _phi[j][i])/h;
                }else if(_U[j][i] == 0){
                    phix = 0.0;
                }
            }


            if(j == 0){
                phiy = (_phi[j+1][i] - _phi[j][i]) / h;
            }else if(j == user->N-1){
                phiy = (_phi[j][i] - _phi[j-1][i]) / h;
            }else if(j != 0 && j != user->N-1){
                if(_V[j][i] > 0){
                    phiy = (_phi[j][i] - _phi[j-1][i])/h;
                }else if(_V[j][i] < 0){
                    phiy = (_phi[j+1][i] - _phi[j][i])/h;
                }else if(_V[j][i] == 0){
                    phiy = 0.0;
                }
            }

                _phi_new[j][i] = _phi[j][i] - dt*(_U[j][i]*phix + _V[j][i]*phiy);
        }
    }

    ierr = DMDAVecRestoreArray(da,user->phi_new,&_phi_new); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->phi,&_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.U,&_U); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user->soln.V,&_V); CHKERRQ(ierr);

    ierr = VecCopy(user->phi_new,user->phi_0); CHKERRQ(ierr);

    return(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeTimeStep"
extern PetscErrorCode ComputeTimeStep(UserContext *user){

    PetscErrorCode ierr;
    PetscReal dts, dtv, dtv1, dtv2, dtc, dtc1, dtc2, dtf, dtf1, dtf2;
    PetscReal rhoWater = 1.0, rhoAir = 0.001226;
    PetscReal muWater = 0.01137, muAir = 0.000178;
    double pi = 3.1415926535;
    PetscReal umax, umin, vmax, vmin,u,v, RHSumax, RHSumin, RHSvmax,RHSvmin,f1,f2;

    dts = sqrt(((rhoWater + rhoAir)*user->Bond)/(8.0*pi))*pow(user->h,1.5);

    dtv1 = (3.0/14.0*user->Re)*rhoWater*(user->h*user->h)/muWater;

    dtv2 = (3.0/14.0*user->Re)*rhoAir*(user->h*user->h)/muAir;

    ierr = VecMax(user->soln.U,NULL,&umax); CHKERRQ(ierr);
    ierr = VecMin(user->soln.U,NULL,&umin); CHKERRQ(ierr);
    ierr = VecMax(user->soln.V,NULL,&vmax); CHKERRQ(ierr);
    ierr = VecMin(user->soln.V,NULL,&vmin); CHKERRQ(ierr);
    ierr = VecMax(user->RHSu,NULL,&RHSumax); CHKERRQ(ierr);
    ierr = VecMin(user->RHSu,NULL,&RHSumin); CHKERRQ(ierr);
    ierr = VecMax(user->RHSv,NULL,&RHSvmax); CHKERRQ(ierr);
    ierr = VecMin(user->RHSv,NULL,&RHSvmin); CHKERRQ(ierr);

    u = fmax(fabs(umax),fabs(umin));
    v = fmax(fabs(vmax),fabs(vmin));

    dtc1 = user->h/u;
    dtc2 = user->h/v;

    dtv = fmin(dtv1,dtv2);
    dtc = fmin(dtc1,dtc2);

    f1 = fmax(fabs(RHSumax),fabs(RHSumin));
    f2 = fmax(fabs(RHSvmax),fabs(RHSvmin));

    dtf1 = 2*user->h/f1;
    dtf2 = 2*user->h/f2;

    dtf = fmin(dtf1,dtf2);

    user->dt_old = user->dt;

    user->dt = 0.5*fmin(fmin(fmin(dts,dtv),dtc),dtf);

    user->time = user->time + user->dt;

    PetscPrintf(PETSC_COMM_WORLD,"Physical Time: %g\n",user->time);
    PetscPrintf(PETSC_COMM_WORLD,"dt: %g\n",user->dt);
    PetscPrintf(PETSC_COMM_WORLD,"umax: %g, vmax: %g\n",u,v);

    return(0);
}




#undef __FUNCT__
#define __FUNCT__ "WriteOutput"
extern PetscErrorCode WriteOutput(UserContext *user){

    PetscErrorCode ierr;

    /* normal to interface */
    //WriteVec(user->normalX,"normalX"); CHKERRQ(ierr);
    //WriteVec(user->normalY,"normalY"); CHKERRQ(ierr);

    /* write intermediate velocity field */
    //ierr = WriteVec(user->solstar.U,"Ustar"); CHKERRQ(ierr);
    //ierr = WriteVec(user->solstar.V,"Vstar"); CHKERRQ(ierr);

    /* RHS of projection step */
    //WriteVec(user->RHS,"b"); CHKERRQ(ierr);

    /* write jacobian matrix */
    //ierr = WriteMat(user->Jac,"Jac"); CHKERRQ(ierr);

    /* write the KSP solution */
    //ierr = WriteVec(user->solstar.P,"KSP"); CHKERRQ(ierr);

    /* velocity at n+1 */
    ierr = WriteVec(user->solnp1.U,"Unp1"); CHKERRQ(ierr);
    ierr = WriteVec(user->solnp1.V,"Vnp1"); CHKERRQ(ierr);
    //ierr = WriteVec(user->solnp1.P,"Pnp1"); CHKERRQ(ierr);

    /* define level set */
    //ierr = WriteVec(user->phi_0,"levelset"); CHKERRQ(ierr);

    /* reinitialize */
    //ierr = WriteVec(user->phi,"phi1"); CHKERRQ(ierr);

    /* define variables */
    //WriteVec(user->H,"Heaviside");
    //WriteVec(user->soln.mu,"Mu");
    //WriteVec(user->soln.rho,"Rho");

    /* advect level set */
    ierr = WriteVec(user->phi_new,"phinp1"); CHKERRQ(ierr);

    return(0);
}





#undef __FUNCT__
#define __FUNCT__ "ConserveMass"
PetscErrorCode ConserveMass(DM da, UserContext *user){

    PetscErrorCode ierr;
    PetscReal h = user->h;
    PetscInt i,j,xs,ys,xm,ym;
    PetscReal **_H, **_phi, delta;
    double pi = 3.1415926535, eps = 2*user->h;

    ierr = DMDAVecGetArray(da,user->H,&_H);
    ierr = DMDAVecGetArray(da,user->phi,&_phi);

    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); CHKERRQ(ierr);

    user->VolumeNew = 0.0; user->Length = 0.0;

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {

            /* compute volume */
            if(user->counter == 0){
                user->Volume0 += (1.0 - _H[j][i])*h*h;
            }

            else if(user->counter > 0){
                user->VolumeNew += (1.0 - _H[j][i])*h*h;

                /* compute length */
                if( fabs(_phi[j][i]) <= 1.5*h){
                    delta = 0.5*(1 + cos(pi*_phi[j][i]/eps))/eps;
                }else if(fabs(_phi[j][i]) > 1.5*h){
                    delta = 0.0;
                }
                user->Length += delta*h*h;

            }
        }
    }

    if(user->counter > 0){
        /* compute dphi */
        user->dphi = (user->VolumeNew - user->Volume0) / user->Length;

        for (j=ys; j<ys+ym; j++) {
            for (i=xs; i<xs+xm; i++) {
                /* update phi */
                _phi[j][i] += user->dphi;
            }
        }

    }


    if(user->counter == 0){
        PetscPrintf(PETSC_COMM_WORLD,"Initial volume: %g\n",user->Volume0);
    }else if(user->counter > 0){
        PetscPrintf(PETSC_COMM_WORLD,"New volume: %g and Length: %g\n",user->VolumeNew, user->Length);
        PetscPrintf(PETSC_COMM_WORLD,"volume change: %g and dphi: %g\n",user->VolumeNew - user->Volume0,user->dphi);
    }

    ierr = DMDAVecRestoreArray(da,user->H,&_H);
    ierr = DMDAVecRestoreArray(da,user->phi,&_phi);


    if(user->counter > 0){
        FILE *volume, *length;
        volume =  fopen("volume.dat","a+");
        length =  fopen("length.dat","a+");
        fprintf(volume,"%g\t",user->VolumeNew);
        fprintf(length,"%g\t",user->Length);
        fclose(volume);
        fclose(length);
    }

    return(0);
}





#undef __FUNCT__
#define __FUNCT__ "Screenshot"
PetscErrorCode Screenshot(UserContext *user,Engine *ep){

    PetscErrorCode ierr;

    /* velocity at n+1 */
    ierr = WriteVec(user->solnp1.U,"Unp1"); CHKERRQ(ierr);
    ierr = WriteVec(user->solnp1.V,"Vnp1"); CHKERRQ(ierr);
    ierr = WriteVec(user->solnp1.P,"Pnp1"); CHKERRQ(ierr);
    /* level set */
    ierr = WriteVec(user->phi_new,"phinp1"); CHKERRQ(ierr);

    char counter[16]; sprintf(counter,"counter = %4d;",user->counter);
    char time[16]; sprintf(time,"time = %.3f;",user->time);

    engEvalString(ep,counter);
    engEvalString(ep,time);

    if(user->counter == 0){
        engEvalString(ep,"cd Output;");
    }
    engEvalString(ep,"screenshot;");


    return(0);
}


