/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the tracer flux.

  Compute the tracer flux along one row of computational 
  zones for the HD and MHD  modules according to Spitzer (1962):
  \f[
     \vec{F}_c = \frac{q}{|\vec{F}_{\rm class}| + q}\vec{F}_{\rm class}
  \f] 
  where \f$ \vec{F}_{\rm class} = \rho \nu_{C} \nabla C\f$ is the
  tracer diffusion flux,
       
  \authors A. Mignone (mignone@to.infn.it)\n
           T. Matsakos
           A. Dutta
  \date    Mar 27, 2023
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "local_pluto.h"

/* ********************************************************************* */
void RHS_TRACER_Flux (double ****TracerField, const Sweep *sweep,  
              double **tracer_flux, int beg, int end, Grid *grid)
/*! 
 * Compute the thermal conduction flux, sweep->par_flx.
 *
 * \param [in]     TracerField   4D array containing the dimensionless 
 *                               3D tracer fields
 * \param [in] sweep   pointer to a Sweep structure
 * \param [out]    tracer_flux  the flux due to the tracer source
 *                         the time step.
 * \param [in]     beg     initial index of computation
 * \param [in]     end     final   index of computation
 * \param [in]     grid    pointer to an array of Grid structures
 *
 * \return This function has no return value.                       
 *********************************************************************** */
{
  int  i, trc, nv;
  double Flux;
  double vi[NVAR];
  double **vc = sweep->vn;
  static double ***gradTRC;
  
  double del_u = 2*g_inputParam[U_FLOW]; // CGS
  double chi   = g_inputParam[LENGTH]*del_u/g_inputParam[REYNOLDS]; 
  
  double nu_dye = (chi/(UNIT_LENGTH*UNIT_VELOCITY));

/* -----------------------------------------------------------
   1. Allocate memory, compute tracer gradient in the
      required direction and set 2 out of the 3 coordinates.
   ----------------------------------------------------------- */

  if (gradTRC == NULL) {
    gradTRC = ARRAY_3D(NTRACER, NMAX_POINT, 3, double);
  }

/* ----------------------------------------------- 
   2. Compute Tracer Difussion Flux (trcflx).
   ----------------------------------------------- */
  
  for (trc = 0; trc < NTRACER; trc++){  
    GetTracerGradient (TracerField[trc], gradTRC[trc], beg, end, grid);
    for (i = beg; i <= end; i++){

    /* -- 3a. Compute interface values -- */

      NVAR_LOOP(nv) vi[nv]  = (vc[i][nv]*grid->dx[g_dir][i] + vc[i+1][nv]*grid->dx[g_dir][i+1])/(grid->dx[g_dir][i]+grid->dx[g_dir][i+1]);
    
    /* -- 3b. Compute the Tracer flux -- */
       
      Flux        = vi[RHO]*nu_dye*gradTRC[trc][i][g_dir];  
      tracer_flux[i][trc] = Flux;
    }
  }
}


/* ********************************************************************* */
void GetTracerGradient (double ***Field, double **gradField, 
                  int beg, int end, Grid *grid)
/*!
 *   Compute the gradient of a 3D scalar quantity C in the direction
 *   given by g_dir.
 *   Return a 1D array (dField/dx, dField/dy, dField/dz) along that direction 
 *   computed at cell interfaces, e.g.
 *
 *   if g_dir == IDIR  --> compute 
 *  
 *    [ dField/dl1, dField/dl2, dField/dl3 ] at interface (i+1/2,j,k)
 *
 *   if g_dir == JDIR  --> compute 
 *  
 *    [ dField/dl1, dField/dl2, dField/dl3 ] at interface (i,j+1/2,k)
 * 
 *   if g_dir == KDIR  --> compute 
 *  
 *    [ dField/dl1, dField/dl2, dField/dl3 ] at interface (i,j,k+1/2)
 *   
 *
 *   Here dl1, dl2 and dl3 are the line element in the 
 *   thre directions:
 *
 *   Cartesian:   {dl1, dl2, dl3} = {dx,       dy, dz}
 *   Cylindrical: {dl1, dl2, dl3} = {dr,       dz, - }
 *   Polar:       {dl1, dl2, dl3} = {dr,   r.dphi, dz} 
 *   Spherical:   {dl1, dl2, dl3} = {dr, r.dtheta, r.sin(theta).dphi}
 *   
 * LAST MODIFIED
 *
 *   14 Apr 2011 by T. Matsakos, A. Mignone
 *   27 Mar 2023 by A. Dutta
 *
 *
 *********************************************************************** */
{
  int  i,j,k;
  double *r, *rp;
  double *inv_dx,  *inv_dy,  *inv_dz;
  double *inv_dxi, *inv_dyi, *inv_dzi;
  double dl1, dl2, dl3, theta, r_1, s_1;
  double dx1, dx2, dx3;
  
  inv_dx  = grid->inv_dx[IDIR]; inv_dxi = grid->inv_dxi[IDIR];
  inv_dy  = grid->inv_dx[JDIR]; inv_dyi = grid->inv_dxi[JDIR];
  inv_dz  = grid->inv_dx[KDIR]; inv_dzi = grid->inv_dxi[KDIR];

  r  = grid->x[IDIR];
  rp = grid->xr[IDIR];

  i = g_i;
  j = g_j;
  k = g_k;

  if (g_dir == IDIR) {

    #if GEOMETRY == SPHERICAL
    theta = grid->x[JDIR][j];
    s_1   = 1.0/sin(theta);
    #endif
    DIM_EXPAND(                 ,
      dl2 = dx2 = inv_dy[j];  ,
      dl3 = dx3 = inv_dz[k];
    )
    for (i = beg; i <= end; i++){
      dl1 = inv_dxi[i];  
      #if GEOMETRY == POLAR && DIMENSIONS >= 2
      dl2 = dx2/rp[i]; 
      #elif GEOMETRY == SPHERICAL
      DIM_EXPAND(                      ,
               dl2 = dx2/rp[i];      ,
               dl3 = dx3*s_1/rp[i];)
      #endif
      DIM_EXPAND( 
        gradField[i][0] = (Field[k][j][i+1] - Field[k][j][i])*dl1;         ,
        gradField[i][1] = 0.25*(  Field[k][j+1][i] + Field[k][j+1][i+1]
                            - Field[k][j-1][i] - Field[k][j-1][i+1])*dl2;  ,
        gradField[i][2] = 0.25*(  Field[k+1][j][i] + Field[k+1][j][i+1] 
                            - Field[k-1][j][i] - Field[k-1][j][i+1])*dl3;
      )
    }

  }else if (g_dir == JDIR) {

    r_1  = 1.0/r[i];
    DIM_EXPAND(
      dl1 = dx1 = inv_dx[i];   ,
                               ,
      dl3 = dx3 = inv_dz[k];
    )
    for (j = beg; j <= end; j++){
      dl2 = inv_dyi[j]; 
      #if GEOMETRY == POLAR
      dl2 *= r_1;
      #elif GEOMETRY == SPHERICAL
      DIM_EXPAND(               ,
               dl2  *= r_1;   ,
               theta = grid->xr[JDIR][j];
               dl3   = dx3*r_1/sin(theta);)
      #endif
      DIM_EXPAND( 
               gradField[j][0] = 0.25*(  Field[k][j][i+1] + Field[k][j+1][i+1]
                                   - Field[k][j][i-1] - Field[k][j+1][i-1])*dl1;   ,
               gradField[j][1] = (Field[k][j+1][i] - Field[k][j][i])*dl2;              ,
               gradField[j][2] = 0.25*(  Field[k+1][j][i] + Field[k+1][j+1][i]
                                   - Field[k-1][j][i] - Field[k-1][j+1][i])*dl3;
      )
    }
  
  }else if (g_dir == KDIR){

    dl1 = inv_dx[i];            
    dl2 = inv_dy[j]; 
    #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
    r_1  = 1.0/r[i];
    dl2 *= r_1; 
    #endif
    #if GEOMETRY == SPHERICAL
    theta = grid->x[JDIR][j];
    s_1   = 1.0/sin(theta);
    #endif

    for (k = beg; k <= end; k++){
      dl3 = inv_dzi[k]; 
      #if GEOMETRY == SPHERICAL
       dl3 *= r_1*s_1;
      #endif
      gradField[k][0] = 0.25*(  Field[k][j][i+1] + Field[k+1][j][i+1]
                          - Field[k][j][i-1] - Field[k+1][j][i-1])*dl1;
      gradField[k][1] = 0.25*(  Field[k][j+1][i] + Field[k+1][j+1][i]
                          - Field[k][j-1][i] - Field[k+1][j-1][i])*dl2;
      gradField[k][2] = (Field[k+1][j][i] - Field[k][j][i])*dl3;
    }
  }
}
