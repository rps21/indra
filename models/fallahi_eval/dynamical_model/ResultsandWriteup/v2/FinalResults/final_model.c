/*   
**   FINALMODEL_cvode.c
**	 
**   Cvode-Mex implementation of BioNetGen model 'FINALMODEL'.
**
**   Code Adapted from templates provided by Mathworks and Sundials.
**   QUESTIONS about the code generator?  Email justinshogg@gmail.com
**
**   Requires the CVODE libraries:  sundials_cvode and sundials_nvecserial.
**   https://computation.llnl.gov/casc/sundials/main.html
**
**-----------------------------------------------------------------------------
**
**   COMPILE in MATLAB:
**   mex -L<path_to_cvode_libraries> -I<path_to_cvode_includes>  ...
**          -lsundials_nvecserial -lsundials_cvode -lm FINALMODEL_cvode.c
**
**   note1: if cvode is in your library path, you can omit path specifications.
**
**   note2: if linker complains about lib stdc++, try removing "-lstdc++"
**     from the mex configuration file "gccopts.sh".  This should be in the
**     matlab bin folder.
** 
**-----------------------------------------------------------------------------
**
**   EXECUTE in MATLAB:
**   [error_status, species_out, observables_out]
**        = FINALMODEL_cvode( timepoints, species_init, parameters )
**
**   timepoints      : column vector of time points returned by integrator.
**   parameters      : row vector of 120 parameters.
**   species_init    : row vector of 81 initial species populations.
**
**   error_status    : 0 if the integrator exits without error, non-zero otherwise.
**   species_out     : species population trajectories
**                        (columns correspond to states, rows correspond to time).
**   observables_out : observable trajectories
**                        (columns correspond to observables, rows correspond to time).
*/

/* Library headers */
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>             /* prototypes for CVODE  */
#include <nvector/nvector_serial.h>  /* serial N_Vector       */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_spgmr.h>       /* prototype for CVSpgmr */

/* Problem Dimensions */
#define __N_PARAMETERS__   120
#define __N_EXPRESSIONS__  120
#define __N_OBSERVABLES__  15
#define __N_RATELAWS__     180
#define __N_SPECIES__      81

/* core function declarations */
void  mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
int   check_flag  ( void *flagvalue, char *funcname, int opt );
void  calc_expressions ( N_Vector expressions, double * parameters );
void  calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions );
void  calc_ratelaws    ( N_Vector ratelaws,  N_Vector species, N_Vector expressions, N_Vector observables );
int   calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data );

/* user-defined function declarations */


/* user-defined function definitions  */


/* Calculate expressions */
void
calc_expressions ( N_Vector expressions, double * parameters )
{
    NV_Ith_S(expressions,0) = parameters[0];
    NV_Ith_S(expressions,1) = parameters[1];
    NV_Ith_S(expressions,2) = parameters[2];
    NV_Ith_S(expressions,3) = parameters[3];
    NV_Ith_S(expressions,4) = parameters[4];
    NV_Ith_S(expressions,5) = parameters[5];
    NV_Ith_S(expressions,6) = parameters[6];
    NV_Ith_S(expressions,7) = parameters[7];
    NV_Ith_S(expressions,8) = parameters[8];
    NV_Ith_S(expressions,9) = parameters[9];
    NV_Ith_S(expressions,10) = parameters[10];
    NV_Ith_S(expressions,11) = parameters[11];
    NV_Ith_S(expressions,12) = parameters[12];
    NV_Ith_S(expressions,13) = parameters[13];
    NV_Ith_S(expressions,14) = parameters[14];
    NV_Ith_S(expressions,15) = parameters[15];
    NV_Ith_S(expressions,16) = parameters[16];
    NV_Ith_S(expressions,17) = parameters[17];
    NV_Ith_S(expressions,18) = parameters[18];
    NV_Ith_S(expressions,19) = parameters[19];
    NV_Ith_S(expressions,20) = parameters[20];
    NV_Ith_S(expressions,21) = parameters[21];
    NV_Ith_S(expressions,22) = parameters[22];
    NV_Ith_S(expressions,23) = parameters[23];
    NV_Ith_S(expressions,24) = parameters[24];
    NV_Ith_S(expressions,25) = parameters[25];
    NV_Ith_S(expressions,26) = parameters[26];
    NV_Ith_S(expressions,27) = parameters[27];
    NV_Ith_S(expressions,28) = parameters[28];
    NV_Ith_S(expressions,29) = parameters[29];
    NV_Ith_S(expressions,30) = parameters[30];
    NV_Ith_S(expressions,31) = parameters[31];
    NV_Ith_S(expressions,32) = parameters[32];
    NV_Ith_S(expressions,33) = parameters[33];
    NV_Ith_S(expressions,34) = parameters[34];
    NV_Ith_S(expressions,35) = parameters[35];
    NV_Ith_S(expressions,36) = parameters[36];
    NV_Ith_S(expressions,37) = parameters[37];
    NV_Ith_S(expressions,38) = parameters[38];
    NV_Ith_S(expressions,39) = parameters[39];
    NV_Ith_S(expressions,40) = parameters[40];
    NV_Ith_S(expressions,41) = parameters[41];
    NV_Ith_S(expressions,42) = parameters[42];
    NV_Ith_S(expressions,43) = parameters[43];
    NV_Ith_S(expressions,44) = parameters[44];
    NV_Ith_S(expressions,45) = parameters[45];
    NV_Ith_S(expressions,46) = parameters[46];
    NV_Ith_S(expressions,47) = parameters[47];
    NV_Ith_S(expressions,48) = parameters[48];
    NV_Ith_S(expressions,49) = parameters[49];
    NV_Ith_S(expressions,50) = parameters[50];
    NV_Ith_S(expressions,51) = parameters[51];
    NV_Ith_S(expressions,52) = parameters[52];
    NV_Ith_S(expressions,53) = parameters[53];
    NV_Ith_S(expressions,54) = parameters[54];
    NV_Ith_S(expressions,55) = parameters[55];
    NV_Ith_S(expressions,56) = parameters[56];
    NV_Ith_S(expressions,57) = parameters[57];
    NV_Ith_S(expressions,58) = parameters[58];
    NV_Ith_S(expressions,59) = parameters[59];
    NV_Ith_S(expressions,60) = parameters[60];
    NV_Ith_S(expressions,61) = parameters[61];
    NV_Ith_S(expressions,62) = parameters[62];
    NV_Ith_S(expressions,63) = parameters[63];
    NV_Ith_S(expressions,64) = parameters[64];
    NV_Ith_S(expressions,65) = parameters[65];
    NV_Ith_S(expressions,66) = parameters[66];
    NV_Ith_S(expressions,67) = parameters[67];
    NV_Ith_S(expressions,68) = parameters[68];
    NV_Ith_S(expressions,69) = parameters[69];
    NV_Ith_S(expressions,70) = parameters[70];
    NV_Ith_S(expressions,71) = parameters[71];
    NV_Ith_S(expressions,72) = parameters[72];
    NV_Ith_S(expressions,73) = parameters[73];
    NV_Ith_S(expressions,74) = parameters[74];
    NV_Ith_S(expressions,75) = parameters[75];
    NV_Ith_S(expressions,76) = parameters[76];
    NV_Ith_S(expressions,77) = parameters[77];
    NV_Ith_S(expressions,78) = parameters[78];
    NV_Ith_S(expressions,79) = parameters[79];
    NV_Ith_S(expressions,80) = parameters[80];
    NV_Ith_S(expressions,81) = parameters[81];
    NV_Ith_S(expressions,82) = parameters[82];
    NV_Ith_S(expressions,83) = parameters[83];
    NV_Ith_S(expressions,84) = parameters[84];
    NV_Ith_S(expressions,85) = parameters[85];
    NV_Ith_S(expressions,86) = parameters[86];
    NV_Ith_S(expressions,87) = parameters[87];
    NV_Ith_S(expressions,88) = parameters[88];
    NV_Ith_S(expressions,89) = parameters[89];
    NV_Ith_S(expressions,90) = parameters[90];
    NV_Ith_S(expressions,91) = parameters[91];
    NV_Ith_S(expressions,92) = parameters[92];
    NV_Ith_S(expressions,93) = parameters[93];
    NV_Ith_S(expressions,94) = parameters[94];
    NV_Ith_S(expressions,95) = parameters[95];
    NV_Ith_S(expressions,96) = parameters[96];
    NV_Ith_S(expressions,97) = parameters[97];
    NV_Ith_S(expressions,98) = parameters[98];
    NV_Ith_S(expressions,99) = parameters[99];
    NV_Ith_S(expressions,100) = parameters[100];
    NV_Ith_S(expressions,101) = parameters[101];
    NV_Ith_S(expressions,102) = parameters[102];
    NV_Ith_S(expressions,103) = parameters[103];
    NV_Ith_S(expressions,104) = parameters[104];
    NV_Ith_S(expressions,105) = parameters[105];
    NV_Ith_S(expressions,106) = parameters[106];
    NV_Ith_S(expressions,107) = parameters[107];
    NV_Ith_S(expressions,108) = parameters[108];
    NV_Ith_S(expressions,109) = parameters[109];
    NV_Ith_S(expressions,110) = parameters[110];
    NV_Ith_S(expressions,111) = parameters[111];
    NV_Ith_S(expressions,112) = parameters[112];
    NV_Ith_S(expressions,113) = parameters[113];
    NV_Ith_S(expressions,114) = parameters[114];
    NV_Ith_S(expressions,115) = parameters[115];
    NV_Ith_S(expressions,116) = parameters[116];
    NV_Ith_S(expressions,117) = parameters[117];
    NV_Ith_S(expressions,118) = parameters[118];
    NV_Ith_S(expressions,119) = parameters[119];
   
}

/* Calculate observables */
void
calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions )
{
    NV_Ith_S(observables,0) = NV_Ith_S(species,29) +NV_Ith_S(species,34) +NV_Ith_S(species,38) +NV_Ith_S(species,39) +NV_Ith_S(species,43) +NV_Ith_S(species,45) +NV_Ith_S(species,46) +NV_Ith_S(species,69) +NV_Ith_S(species,70);
    NV_Ith_S(observables,1) = NV_Ith_S(species,67) +NV_Ith_S(species,69) +NV_Ith_S(species,71) +NV_Ith_S(species,78) +NV_Ith_S(species,79) +NV_Ith_S(species,80);
    NV_Ith_S(observables,2) = NV_Ith_S(species,6) +NV_Ith_S(species,20) +NV_Ith_S(species,21) +NV_Ith_S(species,22);
    NV_Ith_S(observables,3) = NV_Ith_S(species,53) +NV_Ith_S(species,64);
    NV_Ith_S(observables,4) = NV_Ith_S(species,13) +NV_Ith_S(species,34) +NV_Ith_S(species,35);
    NV_Ith_S(observables,5) = NV_Ith_S(species,12);
    NV_Ith_S(observables,6) = NV_Ith_S(species,61);
    NV_Ith_S(observables,7) = NV_Ith_S(species,54) +NV_Ith_S(species,65);
    NV_Ith_S(observables,8) = NV_Ith_S(species,28) +NV_Ith_S(species,36) +NV_Ith_S(species,37) +NV_Ith_S(species,42);
    NV_Ith_S(observables,9) = NV_Ith_S(species,31);
    NV_Ith_S(observables,10) = NV_Ith_S(species,57);
    NV_Ith_S(observables,11) = NV_Ith_S(species,18) +NV_Ith_S(species,24) +NV_Ith_S(species,25) +NV_Ith_S(species,26) +NV_Ith_S(species,32) +NV_Ith_S(species,46) +NV_Ith_S(species,48) +NV_Ith_S(species,50) +NV_Ith_S(species,56) +NV_Ith_S(species,58) +NV_Ith_S(species,63) +NV_Ith_S(species,68) +NV_Ith_S(species,70) +NV_Ith_S(species,72) +NV_Ith_S(species,74) +NV_Ith_S(species,77);
    NV_Ith_S(observables,12) = NV_Ith_S(species,25) +NV_Ith_S(species,46) +NV_Ith_S(species,48) +NV_Ith_S(species,56) +NV_Ith_S(species,63) +NV_Ith_S(species,68) +NV_Ith_S(species,70) +NV_Ith_S(species,72) +NV_Ith_S(species,74) +NV_Ith_S(species,77);
    NV_Ith_S(observables,13) = NV_Ith_S(species,26);
    NV_Ith_S(observables,14) = NV_Ith_S(species,24) +NV_Ith_S(species,32) +NV_Ith_S(species,50) +NV_Ith_S(species,58);

}

/* Calculate ratelaws */
void
calc_ratelaws ( N_Vector ratelaws, N_Vector species, N_Vector expressions, N_Vector observables )
{  
    NV_Ith_S(ratelaws,0) = NV_Ith_S(expressions,85)*NV_Ith_S(species,0)*NV_Ith_S(species,1);
    NV_Ith_S(ratelaws,1) = NV_Ith_S(expressions,43)*NV_Ith_S(species,6)*NV_Ith_S(species,7);
    NV_Ith_S(ratelaws,2) = NV_Ith_S(expressions,46)*NV_Ith_S(species,6)*NV_Ith_S(species,8);
    NV_Ith_S(ratelaws,3) = NV_Ith_S(expressions,49)*NV_Ith_S(species,6)*NV_Ith_S(species,9);
    NV_Ith_S(ratelaws,4) = NV_Ith_S(expressions,109)*NV_Ith_S(species,15)*NV_Ith_S(species,16);
    NV_Ith_S(ratelaws,5) = NV_Ith_S(expressions,23)*NV_Ith_S(species,18)*NV_Ith_S(species,16);
    NV_Ith_S(ratelaws,6) = NV_Ith_S(expressions,25)*NV_Ith_S(species,12);
    NV_Ith_S(ratelaws,7) = NV_Ith_S(expressions,26)*NV_Ith_S(species,13);
    NV_Ith_S(ratelaws,8) = NV_Ith_S(expressions,19)*NV_Ith_S(species,18)*NV_Ith_S(species,5);
    NV_Ith_S(ratelaws,9) = NV_Ith_S(expressions,21)*NV_Ith_S(species,18)*NV_Ith_S(species,6);
    NV_Ith_S(ratelaws,10) = NV_Ith_S(expressions,27)*NV_Ith_S(species,19);
    NV_Ith_S(ratelaws,11) = NV_Ith_S(expressions,86)*NV_Ith_S(species,19);
    NV_Ith_S(ratelaws,12) = NV_Ith_S(expressions,45)*NV_Ith_S(species,20);
    NV_Ith_S(ratelaws,13) = NV_Ith_S(expressions,44)*NV_Ith_S(species,20);
    NV_Ith_S(ratelaws,14) = NV_Ith_S(expressions,48)*NV_Ith_S(species,21);
    NV_Ith_S(ratelaws,15) = NV_Ith_S(expressions,47)*NV_Ith_S(species,21);
    NV_Ith_S(ratelaws,16) = NV_Ith_S(expressions,51)*NV_Ith_S(species,22);
    NV_Ith_S(ratelaws,17) = NV_Ith_S(expressions,50)*NV_Ith_S(species,22);
    NV_Ith_S(ratelaws,18) = 2.0*NV_Ith_S(expressions,110)*NV_Ith_S(species,23);
    NV_Ith_S(ratelaws,19) = NV_Ith_S(expressions,111)*NV_Ith_S(species,23);
    NV_Ith_S(ratelaws,20) = NV_Ith_S(expressions,23)*NV_Ith_S(species,18)*NV_Ith_S(species,23);
    NV_Ith_S(ratelaws,21) = NV_Ith_S(expressions,24)*NV_Ith_S(species,24);
    NV_Ith_S(ratelaws,22) = NV_Ith_S(expressions,20)*NV_Ith_S(species,25);
    NV_Ith_S(ratelaws,23) = NV_Ith_S(expressions,22)*NV_Ith_S(species,26);
    NV_Ith_S(ratelaws,24) = NV_Ith_S(expressions,28)*NV_Ith_S(species,27)*NV_Ith_S(species,2);
    NV_Ith_S(ratelaws,25) = NV_Ith_S(expressions,32)*NV_Ith_S(species,13)*NV_Ith_S(species,29);
    NV_Ith_S(ratelaws,26) = NV_Ith_S(expressions,35)*NV_Ith_S(species,13)*NV_Ith_S(species,30);
    NV_Ith_S(ratelaws,27) = NV_Ith_S(expressions,52)*NV_Ith_S(species,28)*NV_Ith_S(species,8);
    NV_Ith_S(ratelaws,28) = NV_Ith_S(expressions,55)*NV_Ith_S(species,28)*NV_Ith_S(species,9);
    NV_Ith_S(ratelaws,29) = NV_Ith_S(expressions,61)*NV_Ith_S(species,29)*NV_Ith_S(species,10);
    NV_Ith_S(ratelaws,30) = NV_Ith_S(expressions,70)*NV_Ith_S(species,29)*NV_Ith_S(species,11);
    NV_Ith_S(ratelaws,31) = NV_Ith_S(expressions,76)*NV_Ith_S(species,30)*NV_Ith_S(species,10);
    NV_Ith_S(ratelaws,32) = NV_Ith_S(expressions,86)*NV_Ith_S(species,27);
    NV_Ith_S(ratelaws,33) = NV_Ith_S(expressions,97)*NV_Ith_S(species,14)*NV_Ith_S(species,28);
    NV_Ith_S(ratelaws,34) = NV_Ith_S(expressions,100)*NV_Ith_S(species,14)*NV_Ith_S(species,29);
    NV_Ith_S(ratelaws,35) = NV_Ith_S(expressions,103)*NV_Ith_S(species,14)*NV_Ith_S(species,30);
    NV_Ith_S(ratelaws,36) = NV_Ith_S(expressions,58)*NV_Ith_S(species,29)*NV_Ith_S(species,5);
    NV_Ith_S(ratelaws,37) = NV_Ith_S(expressions,58)*NV_Ith_S(species,29)*NV_Ith_S(species,25);
    NV_Ith_S(ratelaws,38) = NV_Ith_S(expressions,73)*NV_Ith_S(species,30)*NV_Ith_S(species,5);
    NV_Ith_S(ratelaws,39) = NV_Ith_S(expressions,73)*NV_Ith_S(species,30)*NV_Ith_S(species,25);
    NV_Ith_S(ratelaws,40) = NV_Ith_S(expressions,110)*NV_Ith_S(species,32);
    NV_Ith_S(ratelaws,41) = NV_Ith_S(expressions,112)*NV_Ith_S(species,31)*NV_Ith_S(species,17);
    NV_Ith_S(ratelaws,42) = NV_Ith_S(expressions,23)*NV_Ith_S(species,18)*NV_Ith_S(species,31);
    NV_Ith_S(ratelaws,43) = NV_Ith_S(expressions,24)*NV_Ith_S(species,32);
    NV_Ith_S(ratelaws,44) = NV_Ith_S(expressions,28)*NV_Ith_S(species,41)*NV_Ith_S(species,2);
    NV_Ith_S(ratelaws,45) = NV_Ith_S(expressions,29)*NV_Ith_S(species,33);
    NV_Ith_S(ratelaws,46) = NV_Ith_S(expressions,30)*NV_Ith_S(species,33)*NV_Ith_S(species,3);
    NV_Ith_S(ratelaws,47) = NV_Ith_S(expressions,34)*NV_Ith_S(species,34);
    NV_Ith_S(ratelaws,48) = NV_Ith_S(expressions,33)*NV_Ith_S(species,34);
    NV_Ith_S(ratelaws,49) = NV_Ith_S(expressions,37)*NV_Ith_S(species,35);
    NV_Ith_S(ratelaws,50) = NV_Ith_S(expressions,36)*NV_Ith_S(species,35);
    NV_Ith_S(ratelaws,51) = NV_Ith_S(expressions,54)*NV_Ith_S(species,36);
    NV_Ith_S(ratelaws,52) = NV_Ith_S(expressions,53)*NV_Ith_S(species,36);
    NV_Ith_S(ratelaws,53) = NV_Ith_S(expressions,57)*NV_Ith_S(species,37);
    NV_Ith_S(ratelaws,54) = NV_Ith_S(expressions,56)*NV_Ith_S(species,37);
    NV_Ith_S(ratelaws,55) = NV_Ith_S(expressions,63)*NV_Ith_S(species,38);
    NV_Ith_S(ratelaws,56) = NV_Ith_S(expressions,62)*NV_Ith_S(species,38);
    NV_Ith_S(ratelaws,57) = NV_Ith_S(expressions,72)*NV_Ith_S(species,39);
    NV_Ith_S(ratelaws,58) = NV_Ith_S(expressions,71)*NV_Ith_S(species,39);
    NV_Ith_S(ratelaws,59) = NV_Ith_S(expressions,78)*NV_Ith_S(species,40);
    NV_Ith_S(ratelaws,60) = NV_Ith_S(expressions,77)*NV_Ith_S(species,40);
    NV_Ith_S(ratelaws,61) = NV_Ith_S(expressions,85)*NV_Ith_S(species,0)*NV_Ith_S(species,41);
    NV_Ith_S(ratelaws,62) = NV_Ith_S(expressions,86)*NV_Ith_S(species,33);
    NV_Ith_S(ratelaws,63) = NV_Ith_S(expressions,99)*NV_Ith_S(species,42);
    NV_Ith_S(ratelaws,64) = NV_Ith_S(expressions,99)*NV_Ith_S(species,42);
    NV_Ith_S(ratelaws,65) = NV_Ith_S(expressions,102)*NV_Ith_S(species,43);
    NV_Ith_S(ratelaws,66) = NV_Ith_S(expressions,101)*NV_Ith_S(species,43);
    NV_Ith_S(ratelaws,67) = NV_Ith_S(expressions,105)*NV_Ith_S(species,44);
    NV_Ith_S(ratelaws,68) = NV_Ith_S(expressions,104)*NV_Ith_S(species,44);
    NV_Ith_S(ratelaws,69) = NV_Ith_S(expressions,47)*NV_Ith_S(species,45);
    NV_Ith_S(ratelaws,70) = NV_Ith_S(expressions,47)*NV_Ith_S(species,46);
    NV_Ith_S(ratelaws,71) = NV_Ith_S(expressions,50)*NV_Ith_S(species,47);
    NV_Ith_S(ratelaws,72) = NV_Ith_S(expressions,50)*NV_Ith_S(species,48);
    NV_Ith_S(ratelaws,73) = NV_Ith_S(expressions,60)*NV_Ith_S(species,45);
    NV_Ith_S(ratelaws,74) = NV_Ith_S(expressions,60)*NV_Ith_S(species,46);
    NV_Ith_S(ratelaws,75) = NV_Ith_S(expressions,59)*NV_Ith_S(species,45);
    NV_Ith_S(ratelaws,76) = NV_Ith_S(expressions,59)*NV_Ith_S(species,46);
    NV_Ith_S(ratelaws,77) = NV_Ith_S(expressions,75)*NV_Ith_S(species,47);
    NV_Ith_S(ratelaws,78) = NV_Ith_S(expressions,75)*NV_Ith_S(species,48);
    NV_Ith_S(ratelaws,79) = NV_Ith_S(expressions,74)*NV_Ith_S(species,47);
    NV_Ith_S(ratelaws,80) = NV_Ith_S(expressions,74)*NV_Ith_S(species,48);
    NV_Ith_S(ratelaws,81) = 2.0*NV_Ith_S(expressions,113)*NV_Ith_S(species,49);
    NV_Ith_S(ratelaws,82) = NV_Ith_S(expressions,114)*NV_Ith_S(species,49);
    NV_Ith_S(ratelaws,83) = NV_Ith_S(expressions,23)*NV_Ith_S(species,18)*NV_Ith_S(species,49);
    NV_Ith_S(ratelaws,84) = NV_Ith_S(expressions,24)*NV_Ith_S(species,50);
    NV_Ith_S(ratelaws,85) = NV_Ith_S(expressions,26)*NV_Ith_S(species,34);
    NV_Ith_S(ratelaws,86) = NV_Ith_S(expressions,26)*NV_Ith_S(species,35);
    NV_Ith_S(ratelaws,87) = NV_Ith_S(expressions,29)*NV_Ith_S(species,51);
    NV_Ith_S(ratelaws,88) = NV_Ith_S(expressions,29)*NV_Ith_S(species,52);
    NV_Ith_S(ratelaws,89) = NV_Ith_S(expressions,30)*NV_Ith_S(species,51)*NV_Ith_S(species,3);
    NV_Ith_S(ratelaws,90) = NV_Ith_S(expressions,31)*NV_Ith_S(species,52);
    NV_Ith_S(ratelaws,91) = NV_Ith_S(expressions,38)*NV_Ith_S(species,52)*NV_Ith_S(species,4);
    NV_Ith_S(ratelaws,92) = NV_Ith_S(expressions,40)*NV_Ith_S(species,53);
    NV_Ith_S(ratelaws,93) = NV_Ith_S(expressions,41)*NV_Ith_S(species,54);
    NV_Ith_S(ratelaws,94) = NV_Ith_S(expressions,85)*NV_Ith_S(species,0)*NV_Ith_S(species,51);
    NV_Ith_S(ratelaws,95) = NV_Ith_S(expressions,86)*NV_Ith_S(species,52);
    NV_Ith_S(ratelaws,96) = NV_Ith_S(expressions,88)*NV_Ith_S(species,14)*NV_Ith_S(species,55);
    NV_Ith_S(ratelaws,97) = NV_Ith_S(expressions,88)*NV_Ith_S(species,14)*NV_Ith_S(species,56);
    NV_Ith_S(ratelaws,98) = NV_Ith_S(expressions,91)*NV_Ith_S(species,14)*NV_Ith_S(species,53);
    NV_Ith_S(ratelaws,99) = NV_Ith_S(expressions,106)*NV_Ith_S(species,14)*NV_Ith_S(species,54);
    NV_Ith_S(ratelaws,100) = NV_Ith_S(expressions,113)*NV_Ith_S(species,58);
    NV_Ith_S(ratelaws,101) = NV_Ith_S(expressions,115)*NV_Ith_S(species,57)*NV_Ith_S(species,10);
    NV_Ith_S(ratelaws,102) = NV_Ith_S(expressions,118)*NV_Ith_S(species,57);
    NV_Ith_S(ratelaws,103) = NV_Ith_S(expressions,24)*NV_Ith_S(species,58);
    NV_Ith_S(ratelaws,104) = NV_Ith_S(expressions,119)*NV_Ith_S(species,53);
    NV_Ith_S(ratelaws,105) = NV_Ith_S(expressions,19)*NV_Ith_S(species,18)*NV_Ith_S(species,55);
    NV_Ith_S(ratelaws,106) = NV_Ith_S(expressions,20)*NV_Ith_S(species,56);
    NV_Ith_S(ratelaws,107) = NV_Ith_S(expressions,28)*NV_Ith_S(species,27)*NV_Ith_S(species,59);
    NV_Ith_S(ratelaws,108) = NV_Ith_S(expressions,28)*NV_Ith_S(species,41)*NV_Ith_S(species,59);
    NV_Ith_S(ratelaws,109) = NV_Ith_S(expressions,29)*NV_Ith_S(species,60);
    NV_Ith_S(ratelaws,110) = NV_Ith_S(expressions,31)*NV_Ith_S(species,59);
    NV_Ith_S(ratelaws,111) = NV_Ith_S(expressions,31)*NV_Ith_S(species,60);
    NV_Ith_S(ratelaws,112) = NV_Ith_S(expressions,38)*NV_Ith_S(species,59)*NV_Ith_S(species,4);
    NV_Ith_S(ratelaws,113) = NV_Ith_S(expressions,38)*NV_Ith_S(species,60)*NV_Ith_S(species,4);
    NV_Ith_S(ratelaws,114) = NV_Ith_S(expressions,40)*NV_Ith_S(species,64);
    NV_Ith_S(ratelaws,115) = NV_Ith_S(expressions,41)*NV_Ith_S(species,65);
    NV_Ith_S(ratelaws,116) = NV_Ith_S(expressions,42)*NV_Ith_S(species,12)*NV_Ith_S(species,61);
    NV_Ith_S(ratelaws,117) = NV_Ith_S(expressions,85)*NV_Ith_S(species,0)*NV_Ith_S(species,60);
    NV_Ith_S(ratelaws,118) = 2.0*NV_Ith_S(expressions,89)*NV_Ith_S(species,62);
    NV_Ith_S(ratelaws,119) = 2.0*NV_Ith_S(expressions,89)*NV_Ith_S(species,63);
    NV_Ith_S(ratelaws,120) = NV_Ith_S(expressions,90)*NV_Ith_S(species,62);
    NV_Ith_S(ratelaws,121) = NV_Ith_S(expressions,90)*NV_Ith_S(species,63);
    NV_Ith_S(ratelaws,122) = NV_Ith_S(expressions,93)*NV_Ith_S(species,64);
    NV_Ith_S(ratelaws,123) = NV_Ith_S(expressions,92)*NV_Ith_S(species,64);
    NV_Ith_S(ratelaws,124) = NV_Ith_S(expressions,108)*NV_Ith_S(species,65);
    NV_Ith_S(ratelaws,125) = NV_Ith_S(expressions,107)*NV_Ith_S(species,65);
    NV_Ith_S(ratelaws,126) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,5);
    NV_Ith_S(ratelaws,127) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,25);
    NV_Ith_S(ratelaws,128) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,45);
    NV_Ith_S(ratelaws,129) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,46);
    NV_Ith_S(ratelaws,130) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,47);
    NV_Ith_S(ratelaws,131) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,48);
    NV_Ith_S(ratelaws,132) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,55);
    NV_Ith_S(ratelaws,133) = NV_Ith_S(expressions,39)*NV_Ith_S(species,61)*NV_Ith_S(species,56);
    NV_Ith_S(ratelaws,134) = NV_Ith_S(expressions,116)*NV_Ith_S(species,66);
    NV_Ith_S(ratelaws,135) = NV_Ith_S(expressions,117)*NV_Ith_S(species,66);
    NV_Ith_S(ratelaws,136) = NV_Ith_S(expressions,115)*NV_Ith_S(species,66);
    NV_Ith_S(ratelaws,137) = NV_Ith_S(expressions,119)*NV_Ith_S(species,64);
    NV_Ith_S(ratelaws,138) = NV_Ith_S(expressions,19)*NV_Ith_S(species,18)*NV_Ith_S(species,62);
    NV_Ith_S(ratelaws,139) = NV_Ith_S(expressions,20)*NV_Ith_S(species,63);
    NV_Ith_S(ratelaws,140) = NV_Ith_S(expressions,88)*NV_Ith_S(species,14)*NV_Ith_S(species,73);
    NV_Ith_S(ratelaws,141) = NV_Ith_S(expressions,88)*NV_Ith_S(species,14)*NV_Ith_S(species,74);
    NV_Ith_S(ratelaws,142) = NV_Ith_S(expressions,92)*NV_Ith_S(species,75);
    NV_Ith_S(ratelaws,143) = NV_Ith_S(expressions,43)*NV_Ith_S(species,67)*NV_Ith_S(species,7);
    NV_Ith_S(ratelaws,144) = NV_Ith_S(expressions,46)*NV_Ith_S(species,67)*NV_Ith_S(species,8);
    NV_Ith_S(ratelaws,145) = NV_Ith_S(expressions,47)*NV_Ith_S(species,69);
    NV_Ith_S(ratelaws,146) = NV_Ith_S(expressions,47)*NV_Ith_S(species,70);
    NV_Ith_S(ratelaws,147) = NV_Ith_S(expressions,49)*NV_Ith_S(species,67)*NV_Ith_S(species,9);
    NV_Ith_S(ratelaws,148) = NV_Ith_S(expressions,50)*NV_Ith_S(species,71);
    NV_Ith_S(ratelaws,149) = NV_Ith_S(expressions,50)*NV_Ith_S(species,72);
    NV_Ith_S(ratelaws,150) = NV_Ith_S(expressions,58)*NV_Ith_S(species,29)*NV_Ith_S(species,67);
    NV_Ith_S(ratelaws,151) = NV_Ith_S(expressions,58)*NV_Ith_S(species,29)*NV_Ith_S(species,68);
    NV_Ith_S(ratelaws,152) = NV_Ith_S(expressions,60)*NV_Ith_S(species,69);
    NV_Ith_S(ratelaws,153) = NV_Ith_S(expressions,60)*NV_Ith_S(species,70);
    NV_Ith_S(ratelaws,154) = NV_Ith_S(expressions,59)*NV_Ith_S(species,69);
    NV_Ith_S(ratelaws,155) = NV_Ith_S(expressions,59)*NV_Ith_S(species,70);
    NV_Ith_S(ratelaws,156) = NV_Ith_S(expressions,73)*NV_Ith_S(species,30)*NV_Ith_S(species,67);
    NV_Ith_S(ratelaws,157) = NV_Ith_S(expressions,73)*NV_Ith_S(species,30)*NV_Ith_S(species,68);
    NV_Ith_S(ratelaws,158) = NV_Ith_S(expressions,75)*NV_Ith_S(species,71);
    NV_Ith_S(ratelaws,159) = NV_Ith_S(expressions,75)*NV_Ith_S(species,72);
    NV_Ith_S(ratelaws,160) = NV_Ith_S(expressions,74)*NV_Ith_S(species,71);
    NV_Ith_S(ratelaws,161) = NV_Ith_S(expressions,74)*NV_Ith_S(species,72);
    NV_Ith_S(ratelaws,162) = NV_Ith_S(expressions,19)*NV_Ith_S(species,18)*NV_Ith_S(species,67);
    NV_Ith_S(ratelaws,163) = NV_Ith_S(expressions,19)*NV_Ith_S(species,18)*NV_Ith_S(species,73);
    NV_Ith_S(ratelaws,164) = NV_Ith_S(expressions,20)*NV_Ith_S(species,68);
    NV_Ith_S(ratelaws,165) = NV_Ith_S(expressions,20)*NV_Ith_S(species,74);
    NV_Ith_S(ratelaws,166) = 2.0*NV_Ith_S(expressions,89)*NV_Ith_S(species,76);
    NV_Ith_S(ratelaws,167) = 2.0*NV_Ith_S(expressions,89)*NV_Ith_S(species,77);
    NV_Ith_S(ratelaws,168) = NV_Ith_S(expressions,90)*NV_Ith_S(species,76);
    NV_Ith_S(ratelaws,169) = NV_Ith_S(expressions,90)*NV_Ith_S(species,77);
    NV_Ith_S(ratelaws,170) = NV_Ith_S(expressions,45)*NV_Ith_S(species,78);
    NV_Ith_S(ratelaws,171) = NV_Ith_S(expressions,44)*NV_Ith_S(species,78);
    NV_Ith_S(ratelaws,172) = NV_Ith_S(expressions,48)*NV_Ith_S(species,79);
    NV_Ith_S(ratelaws,173) = NV_Ith_S(expressions,47)*NV_Ith_S(species,79);
    NV_Ith_S(ratelaws,174) = NV_Ith_S(expressions,51)*NV_Ith_S(species,80);
    NV_Ith_S(ratelaws,175) = NV_Ith_S(expressions,50)*NV_Ith_S(species,80);
    NV_Ith_S(ratelaws,176) = NV_Ith_S(expressions,59)*NV_Ith_S(species,79);
    NV_Ith_S(ratelaws,177) = NV_Ith_S(expressions,74)*NV_Ith_S(species,80);
    NV_Ith_S(ratelaws,178) = NV_Ith_S(expressions,19)*NV_Ith_S(species,18)*NV_Ith_S(species,76);
    NV_Ith_S(ratelaws,179) = NV_Ith_S(expressions,20)*NV_Ith_S(species,77);

}


/* Calculate species derivatives */
int
calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data )
{
    int         return_val;
    N_Vector *  temp_data;
    
    N_Vector    expressions;
    N_Vector    observables;
    N_Vector    ratelaws;

    /* cast temp_data */
    temp_data = (N_Vector*)f_data;
     
    /* sget ratelaws Vector */
    expressions = temp_data[0];
    observables = temp_data[1];
    ratelaws    = temp_data[2];
       
    /* calculate observables */
    calc_observables( observables, species, expressions );
    
    /* calculate ratelaws */
    calc_ratelaws( ratelaws, species, expressions, observables );
                        
    /* calculate derivatives */
    NV_Ith_S(Dspecies,0) = -NV_Ith_S(ratelaws,0) +NV_Ith_S(ratelaws,11) +NV_Ith_S(ratelaws,32) -NV_Ith_S(ratelaws,61) +NV_Ith_S(ratelaws,62) -NV_Ith_S(ratelaws,94) +NV_Ith_S(ratelaws,95) -NV_Ith_S(ratelaws,117);
    NV_Ith_S(Dspecies,1) = -NV_Ith_S(ratelaws,0) +NV_Ith_S(ratelaws,11);
    NV_Ith_S(Dspecies,2) = -NV_Ith_S(ratelaws,24) -NV_Ith_S(ratelaws,44) +NV_Ith_S(ratelaws,45) +NV_Ith_S(ratelaws,87) +NV_Ith_S(ratelaws,110);
    NV_Ith_S(Dspecies,3) = -NV_Ith_S(ratelaws,46) -NV_Ith_S(ratelaws,89) +NV_Ith_S(ratelaws,90) +NV_Ith_S(ratelaws,110) +NV_Ith_S(ratelaws,111);
    NV_Ith_S(Dspecies,4) = -NV_Ith_S(ratelaws,91) -NV_Ith_S(ratelaws,112) -NV_Ith_S(ratelaws,113) +NV_Ith_S(ratelaws,116);
    NV_Ith_S(Dspecies,5) = -NV_Ith_S(ratelaws,8) +NV_Ith_S(ratelaws,22) -NV_Ith_S(ratelaws,36) -NV_Ith_S(ratelaws,38) +NV_Ith_S(ratelaws,69) +NV_Ith_S(ratelaws,71) +NV_Ith_S(ratelaws,75) +NV_Ith_S(ratelaws,79) +NV_Ith_S(ratelaws,120) -NV_Ith_S(ratelaws,126);
    NV_Ith_S(Dspecies,6) = -NV_Ith_S(ratelaws,1) -NV_Ith_S(ratelaws,2) -NV_Ith_S(ratelaws,3) -NV_Ith_S(ratelaws,9) +NV_Ith_S(ratelaws,12) +NV_Ith_S(ratelaws,13) +NV_Ith_S(ratelaws,14) +NV_Ith_S(ratelaws,15) +NV_Ith_S(ratelaws,16) +NV_Ith_S(ratelaws,17) +NV_Ith_S(ratelaws,23);
    NV_Ith_S(Dspecies,7) = -NV_Ith_S(ratelaws,1) +NV_Ith_S(ratelaws,13) +NV_Ith_S(ratelaws,63) -NV_Ith_S(ratelaws,143) +NV_Ith_S(ratelaws,171);
    NV_Ith_S(Dspecies,8) = -NV_Ith_S(ratelaws,2) +NV_Ith_S(ratelaws,15) -NV_Ith_S(ratelaws,27) +NV_Ith_S(ratelaws,47) +NV_Ith_S(ratelaws,52) +NV_Ith_S(ratelaws,65) -NV_Ith_S(ratelaws,144) +NV_Ith_S(ratelaws,173) +NV_Ith_S(ratelaws,176);
    NV_Ith_S(Dspecies,9) = -NV_Ith_S(ratelaws,3) +NV_Ith_S(ratelaws,17) -NV_Ith_S(ratelaws,28) +NV_Ith_S(ratelaws,49) +NV_Ith_S(ratelaws,54) +NV_Ith_S(ratelaws,67) -NV_Ith_S(ratelaws,147) +NV_Ith_S(ratelaws,175) +NV_Ith_S(ratelaws,177);
    NV_Ith_S(Dspecies,10) = -NV_Ith_S(ratelaws,29) -NV_Ith_S(ratelaws,31) +NV_Ith_S(ratelaws,56) +NV_Ith_S(ratelaws,60) -NV_Ith_S(ratelaws,101) +NV_Ith_S(ratelaws,104) +NV_Ith_S(ratelaws,122) +NV_Ith_S(ratelaws,134) +NV_Ith_S(ratelaws,136) +NV_Ith_S(ratelaws,142);
    NV_Ith_S(Dspecies,11) = -NV_Ith_S(ratelaws,30) +NV_Ith_S(ratelaws,58) +NV_Ith_S(ratelaws,124);
    NV_Ith_S(Dspecies,12) = -NV_Ith_S(ratelaws,6) +NV_Ith_S(ratelaws,93) +NV_Ith_S(ratelaws,115);
    NV_Ith_S(Dspecies,13) = -NV_Ith_S(ratelaws,7) -NV_Ith_S(ratelaws,25) -NV_Ith_S(ratelaws,26) +NV_Ith_S(ratelaws,47) +NV_Ith_S(ratelaws,48) +NV_Ith_S(ratelaws,49) +NV_Ith_S(ratelaws,50) +NV_Ith_S(ratelaws,92) +NV_Ith_S(ratelaws,114);
    NV_Ith_S(Dspecies,14) = -NV_Ith_S(ratelaws,33) -NV_Ith_S(ratelaws,34) -NV_Ith_S(ratelaws,35) +NV_Ith_S(ratelaws,63) +NV_Ith_S(ratelaws,64) +NV_Ith_S(ratelaws,65) +NV_Ith_S(ratelaws,66) +NV_Ith_S(ratelaws,67) +NV_Ith_S(ratelaws,68) -NV_Ith_S(ratelaws,96) -NV_Ith_S(ratelaws,97) -NV_Ith_S(ratelaws,98) -NV_Ith_S(ratelaws,99) +NV_Ith_S(ratelaws,118) +NV_Ith_S(ratelaws,119) +NV_Ith_S(ratelaws,120) +NV_Ith_S(ratelaws,121) +NV_Ith_S(ratelaws,122) +NV_Ith_S(ratelaws,123) +NV_Ith_S(ratelaws,124) +NV_Ith_S(ratelaws,125) -NV_Ith_S(ratelaws,140) -NV_Ith_S(ratelaws,141) +NV_Ith_S(ratelaws,142) +NV_Ith_S(ratelaws,166) +NV_Ith_S(ratelaws,167) +NV_Ith_S(ratelaws,168) +NV_Ith_S(ratelaws,169);
    NV_Ith_S(Dspecies,15) = -NV_Ith_S(ratelaws,4) +NV_Ith_S(ratelaws,18) +NV_Ith_S(ratelaws,19) +NV_Ith_S(ratelaws,40);
    NV_Ith_S(Dspecies,16) = -NV_Ith_S(ratelaws,4) -NV_Ith_S(ratelaws,5) +NV_Ith_S(ratelaws,18) +NV_Ith_S(ratelaws,21);
    NV_Ith_S(Dspecies,17) = -NV_Ith_S(ratelaws,41) +NV_Ith_S(ratelaws,81) +NV_Ith_S(ratelaws,100) +NV_Ith_S(ratelaws,102);
    NV_Ith_S(Dspecies,18) = -NV_Ith_S(ratelaws,5) -NV_Ith_S(ratelaws,8) -NV_Ith_S(ratelaws,9) -NV_Ith_S(ratelaws,20) +NV_Ith_S(ratelaws,21) +NV_Ith_S(ratelaws,22) +NV_Ith_S(ratelaws,23) -NV_Ith_S(ratelaws,42) +NV_Ith_S(ratelaws,43) -NV_Ith_S(ratelaws,83) +NV_Ith_S(ratelaws,84) +NV_Ith_S(ratelaws,103) -NV_Ith_S(ratelaws,105) +NV_Ith_S(ratelaws,106) -NV_Ith_S(ratelaws,138) +NV_Ith_S(ratelaws,139) -NV_Ith_S(ratelaws,162) -NV_Ith_S(ratelaws,163) +NV_Ith_S(ratelaws,164) +NV_Ith_S(ratelaws,165) -NV_Ith_S(ratelaws,178) +NV_Ith_S(ratelaws,179);
    NV_Ith_S(Dspecies,19) = NV_Ith_S(ratelaws,0) -NV_Ith_S(ratelaws,10) -NV_Ith_S(ratelaws,11);
    NV_Ith_S(Dspecies,20) = NV_Ith_S(ratelaws,1) -NV_Ith_S(ratelaws,12) -NV_Ith_S(ratelaws,13);
    NV_Ith_S(Dspecies,21) = NV_Ith_S(ratelaws,2) -NV_Ith_S(ratelaws,14) -NV_Ith_S(ratelaws,15);
    NV_Ith_S(Dspecies,22) = NV_Ith_S(ratelaws,3) -NV_Ith_S(ratelaws,16) -NV_Ith_S(ratelaws,17);
    NV_Ith_S(Dspecies,23) = NV_Ith_S(ratelaws,4) -NV_Ith_S(ratelaws,18) -NV_Ith_S(ratelaws,19) -NV_Ith_S(ratelaws,20) +NV_Ith_S(ratelaws,43);
    NV_Ith_S(Dspecies,24) = NV_Ith_S(ratelaws,5) -NV_Ith_S(ratelaws,21) +NV_Ith_S(ratelaws,40);
    NV_Ith_S(Dspecies,25) = NV_Ith_S(ratelaws,8) -NV_Ith_S(ratelaws,22) -NV_Ith_S(ratelaws,37) -NV_Ith_S(ratelaws,39) +NV_Ith_S(ratelaws,70) +NV_Ith_S(ratelaws,72) +NV_Ith_S(ratelaws,76) +NV_Ith_S(ratelaws,80) +NV_Ith_S(ratelaws,121) -NV_Ith_S(ratelaws,127);
    NV_Ith_S(Dspecies,26) = NV_Ith_S(ratelaws,9) -NV_Ith_S(ratelaws,23);
    NV_Ith_S(Dspecies,27) = NV_Ith_S(ratelaws,10) -NV_Ith_S(ratelaws,24) -NV_Ith_S(ratelaws,32) +NV_Ith_S(ratelaws,45) +NV_Ith_S(ratelaws,61) +NV_Ith_S(ratelaws,88) -NV_Ith_S(ratelaws,107);
    NV_Ith_S(Dspecies,28) = NV_Ith_S(ratelaws,12) -NV_Ith_S(ratelaws,27) -NV_Ith_S(ratelaws,28) -NV_Ith_S(ratelaws,33) +NV_Ith_S(ratelaws,51) +NV_Ith_S(ratelaws,52) +NV_Ith_S(ratelaws,53) +NV_Ith_S(ratelaws,54) +NV_Ith_S(ratelaws,64) +NV_Ith_S(ratelaws,170);
    NV_Ith_S(Dspecies,29) = NV_Ith_S(ratelaws,14) -NV_Ith_S(ratelaws,25) -NV_Ith_S(ratelaws,29) -NV_Ith_S(ratelaws,30) -NV_Ith_S(ratelaws,34) -NV_Ith_S(ratelaws,36) -NV_Ith_S(ratelaws,37) +NV_Ith_S(ratelaws,48) +NV_Ith_S(ratelaws,51) +NV_Ith_S(ratelaws,55) +NV_Ith_S(ratelaws,56) +NV_Ith_S(ratelaws,57) +NV_Ith_S(ratelaws,58) +NV_Ith_S(ratelaws,66) +NV_Ith_S(ratelaws,69) +NV_Ith_S(ratelaws,70) +NV_Ith_S(ratelaws,73) +NV_Ith_S(ratelaws,74) +NV_Ith_S(ratelaws,75) +NV_Ith_S(ratelaws,76) +NV_Ith_S(ratelaws,145) +NV_Ith_S(ratelaws,146) -NV_Ith_S(ratelaws,150) -NV_Ith_S(ratelaws,151) +NV_Ith_S(ratelaws,152) +NV_Ith_S(ratelaws,153) +NV_Ith_S(ratelaws,154) +NV_Ith_S(ratelaws,155) +NV_Ith_S(ratelaws,172);
    NV_Ith_S(Dspecies,30) = NV_Ith_S(ratelaws,16) -NV_Ith_S(ratelaws,26) -NV_Ith_S(ratelaws,31) -NV_Ith_S(ratelaws,35) -NV_Ith_S(ratelaws,38) -NV_Ith_S(ratelaws,39) +NV_Ith_S(ratelaws,50) +NV_Ith_S(ratelaws,53) +NV_Ith_S(ratelaws,59) +NV_Ith_S(ratelaws,60) +NV_Ith_S(ratelaws,68) +NV_Ith_S(ratelaws,71) +NV_Ith_S(ratelaws,72) +NV_Ith_S(ratelaws,77) +NV_Ith_S(ratelaws,78) +NV_Ith_S(ratelaws,79) +NV_Ith_S(ratelaws,80) +NV_Ith_S(ratelaws,148) +NV_Ith_S(ratelaws,149) -NV_Ith_S(ratelaws,156) -NV_Ith_S(ratelaws,157) +NV_Ith_S(ratelaws,158) +NV_Ith_S(ratelaws,159) +NV_Ith_S(ratelaws,160) +NV_Ith_S(ratelaws,161) +NV_Ith_S(ratelaws,174);
    NV_Ith_S(Dspecies,31) = NV_Ith_S(ratelaws,19) -NV_Ith_S(ratelaws,41) -NV_Ith_S(ratelaws,42) +NV_Ith_S(ratelaws,81) +NV_Ith_S(ratelaws,82) +NV_Ith_S(ratelaws,84);
    NV_Ith_S(Dspecies,32) = NV_Ith_S(ratelaws,20) -NV_Ith_S(ratelaws,40) -NV_Ith_S(ratelaws,43);
    NV_Ith_S(Dspecies,33) = NV_Ith_S(ratelaws,24) -NV_Ith_S(ratelaws,45) -NV_Ith_S(ratelaws,46) -NV_Ith_S(ratelaws,62) +NV_Ith_S(ratelaws,90) +NV_Ith_S(ratelaws,94);
    NV_Ith_S(Dspecies,34) = NV_Ith_S(ratelaws,25) -NV_Ith_S(ratelaws,47) -NV_Ith_S(ratelaws,48) -NV_Ith_S(ratelaws,85);
    NV_Ith_S(Dspecies,35) = NV_Ith_S(ratelaws,26) -NV_Ith_S(ratelaws,49) -NV_Ith_S(ratelaws,50) -NV_Ith_S(ratelaws,86);
    NV_Ith_S(Dspecies,36) = NV_Ith_S(ratelaws,27) -NV_Ith_S(ratelaws,51) -NV_Ith_S(ratelaws,52);
    NV_Ith_S(Dspecies,37) = NV_Ith_S(ratelaws,28) -NV_Ith_S(ratelaws,53) -NV_Ith_S(ratelaws,54);
    NV_Ith_S(Dspecies,38) = NV_Ith_S(ratelaws,29) -NV_Ith_S(ratelaws,55) -NV_Ith_S(ratelaws,56);
    NV_Ith_S(Dspecies,39) = NV_Ith_S(ratelaws,30) -NV_Ith_S(ratelaws,57) -NV_Ith_S(ratelaws,58);
    NV_Ith_S(Dspecies,40) = NV_Ith_S(ratelaws,31) -NV_Ith_S(ratelaws,59) -NV_Ith_S(ratelaws,60);
    NV_Ith_S(Dspecies,41) = NV_Ith_S(ratelaws,32) -NV_Ith_S(ratelaws,44) -NV_Ith_S(ratelaws,61) +NV_Ith_S(ratelaws,87) -NV_Ith_S(ratelaws,108) +NV_Ith_S(ratelaws,109);
    NV_Ith_S(Dspecies,42) = NV_Ith_S(ratelaws,33) -NV_Ith_S(ratelaws,63) -NV_Ith_S(ratelaws,64);
    NV_Ith_S(Dspecies,43) = NV_Ith_S(ratelaws,34) -NV_Ith_S(ratelaws,65) -NV_Ith_S(ratelaws,66);
    NV_Ith_S(Dspecies,44) = NV_Ith_S(ratelaws,35) -NV_Ith_S(ratelaws,67) -NV_Ith_S(ratelaws,68);
    NV_Ith_S(Dspecies,45) = NV_Ith_S(ratelaws,36) -NV_Ith_S(ratelaws,69) -NV_Ith_S(ratelaws,73) -NV_Ith_S(ratelaws,75) -NV_Ith_S(ratelaws,128);
    NV_Ith_S(Dspecies,46) = NV_Ith_S(ratelaws,37) -NV_Ith_S(ratelaws,70) -NV_Ith_S(ratelaws,74) -NV_Ith_S(ratelaws,76) -NV_Ith_S(ratelaws,129);
    NV_Ith_S(Dspecies,47) = NV_Ith_S(ratelaws,38) -NV_Ith_S(ratelaws,71) -NV_Ith_S(ratelaws,77) -NV_Ith_S(ratelaws,79) -NV_Ith_S(ratelaws,130);
    NV_Ith_S(Dspecies,48) = NV_Ith_S(ratelaws,39) -NV_Ith_S(ratelaws,72) -NV_Ith_S(ratelaws,78) -NV_Ith_S(ratelaws,80) -NV_Ith_S(ratelaws,131);
    NV_Ith_S(Dspecies,49) = NV_Ith_S(ratelaws,41) -NV_Ith_S(ratelaws,81) -NV_Ith_S(ratelaws,82) -NV_Ith_S(ratelaws,83) +NV_Ith_S(ratelaws,103);
    NV_Ith_S(Dspecies,50) = NV_Ith_S(ratelaws,42) -NV_Ith_S(ratelaws,84) +NV_Ith_S(ratelaws,100);
    NV_Ith_S(Dspecies,51) = NV_Ith_S(ratelaws,44) +NV_Ith_S(ratelaws,62) -NV_Ith_S(ratelaws,87) -NV_Ith_S(ratelaws,89) -NV_Ith_S(ratelaws,94) +NV_Ith_S(ratelaws,111);
    NV_Ith_S(Dspecies,52) = NV_Ith_S(ratelaws,46) -NV_Ith_S(ratelaws,88) -NV_Ith_S(ratelaws,90) -NV_Ith_S(ratelaws,95) +NV_Ith_S(ratelaws,107) +NV_Ith_S(ratelaws,117);
    NV_Ith_S(Dspecies,53) = NV_Ith_S(ratelaws,55) +NV_Ith_S(ratelaws,59) -NV_Ith_S(ratelaws,98) -NV_Ith_S(ratelaws,104) +NV_Ith_S(ratelaws,123) +NV_Ith_S(ratelaws,135);
    NV_Ith_S(Dspecies,54) = NV_Ith_S(ratelaws,57) -NV_Ith_S(ratelaws,99) +NV_Ith_S(ratelaws,125);
    NV_Ith_S(Dspecies,55) = NV_Ith_S(ratelaws,73) +NV_Ith_S(ratelaws,77) -NV_Ith_S(ratelaws,96) -NV_Ith_S(ratelaws,105) +NV_Ith_S(ratelaws,106) +NV_Ith_S(ratelaws,118) -NV_Ith_S(ratelaws,132);
    NV_Ith_S(Dspecies,56) = NV_Ith_S(ratelaws,74) +NV_Ith_S(ratelaws,78) -NV_Ith_S(ratelaws,97) +NV_Ith_S(ratelaws,105) -NV_Ith_S(ratelaws,106) +NV_Ith_S(ratelaws,119) -NV_Ith_S(ratelaws,133);
    NV_Ith_S(Dspecies,57) = NV_Ith_S(ratelaws,82) -NV_Ith_S(ratelaws,101) -NV_Ith_S(ratelaws,102) +NV_Ith_S(ratelaws,134) +NV_Ith_S(ratelaws,135) +NV_Ith_S(ratelaws,136);
    NV_Ith_S(Dspecies,58) = NV_Ith_S(ratelaws,83) -NV_Ith_S(ratelaws,100) -NV_Ith_S(ratelaws,103);
    NV_Ith_S(Dspecies,59) = NV_Ith_S(ratelaws,88) -NV_Ith_S(ratelaws,107) -NV_Ith_S(ratelaws,108) +NV_Ith_S(ratelaws,109) -NV_Ith_S(ratelaws,110);
    NV_Ith_S(Dspecies,60) = NV_Ith_S(ratelaws,89) +NV_Ith_S(ratelaws,95) +NV_Ith_S(ratelaws,108) -NV_Ith_S(ratelaws,109) -NV_Ith_S(ratelaws,111) -NV_Ith_S(ratelaws,117);
    NV_Ith_S(Dspecies,61) = NV_Ith_S(ratelaws,91) +NV_Ith_S(ratelaws,112) +NV_Ith_S(ratelaws,113) -NV_Ith_S(ratelaws,116);
    NV_Ith_S(Dspecies,62) = NV_Ith_S(ratelaws,96) -NV_Ith_S(ratelaws,118) -NV_Ith_S(ratelaws,120) -NV_Ith_S(ratelaws,138) +NV_Ith_S(ratelaws,139);
    NV_Ith_S(Dspecies,63) = NV_Ith_S(ratelaws,97) -NV_Ith_S(ratelaws,119) -NV_Ith_S(ratelaws,121) +NV_Ith_S(ratelaws,138) -NV_Ith_S(ratelaws,139);
    NV_Ith_S(Dspecies,64) = NV_Ith_S(ratelaws,98) -NV_Ith_S(ratelaws,122) -NV_Ith_S(ratelaws,123) -NV_Ith_S(ratelaws,137);
    NV_Ith_S(Dspecies,65) = NV_Ith_S(ratelaws,99) -NV_Ith_S(ratelaws,124) -NV_Ith_S(ratelaws,125);
    NV_Ith_S(Dspecies,66) = NV_Ith_S(ratelaws,101) -NV_Ith_S(ratelaws,134) -NV_Ith_S(ratelaws,135) -NV_Ith_S(ratelaws,136);
    NV_Ith_S(Dspecies,67) = NV_Ith_S(ratelaws,126) -NV_Ith_S(ratelaws,143) -NV_Ith_S(ratelaws,144) +NV_Ith_S(ratelaws,145) -NV_Ith_S(ratelaws,147) +NV_Ith_S(ratelaws,148) -NV_Ith_S(ratelaws,150) +NV_Ith_S(ratelaws,154) -NV_Ith_S(ratelaws,156) +NV_Ith_S(ratelaws,160) -NV_Ith_S(ratelaws,162) +NV_Ith_S(ratelaws,164) +NV_Ith_S(ratelaws,168) +NV_Ith_S(ratelaws,170) +NV_Ith_S(ratelaws,171) +NV_Ith_S(ratelaws,172) +NV_Ith_S(ratelaws,173) +NV_Ith_S(ratelaws,174) +NV_Ith_S(ratelaws,175) +NV_Ith_S(ratelaws,176) +NV_Ith_S(ratelaws,177);
    NV_Ith_S(Dspecies,68) = NV_Ith_S(ratelaws,127) +NV_Ith_S(ratelaws,146) +NV_Ith_S(ratelaws,149) -NV_Ith_S(ratelaws,151) +NV_Ith_S(ratelaws,155) -NV_Ith_S(ratelaws,157) +NV_Ith_S(ratelaws,161) +NV_Ith_S(ratelaws,162) -NV_Ith_S(ratelaws,164) +NV_Ith_S(ratelaws,169);
    NV_Ith_S(Dspecies,69) = NV_Ith_S(ratelaws,128) -NV_Ith_S(ratelaws,145) +NV_Ith_S(ratelaws,150) -NV_Ith_S(ratelaws,152) -NV_Ith_S(ratelaws,154);
    NV_Ith_S(Dspecies,70) = NV_Ith_S(ratelaws,129) -NV_Ith_S(ratelaws,146) +NV_Ith_S(ratelaws,151) -NV_Ith_S(ratelaws,153) -NV_Ith_S(ratelaws,155);
    NV_Ith_S(Dspecies,71) = NV_Ith_S(ratelaws,130) -NV_Ith_S(ratelaws,148) +NV_Ith_S(ratelaws,156) -NV_Ith_S(ratelaws,158) -NV_Ith_S(ratelaws,160);
    NV_Ith_S(Dspecies,72) = NV_Ith_S(ratelaws,131) -NV_Ith_S(ratelaws,149) +NV_Ith_S(ratelaws,157) -NV_Ith_S(ratelaws,159) -NV_Ith_S(ratelaws,161);
    NV_Ith_S(Dspecies,73) = NV_Ith_S(ratelaws,132) -NV_Ith_S(ratelaws,140) +NV_Ith_S(ratelaws,152) +NV_Ith_S(ratelaws,158) -NV_Ith_S(ratelaws,163) +NV_Ith_S(ratelaws,165) +NV_Ith_S(ratelaws,166);
    NV_Ith_S(Dspecies,74) = NV_Ith_S(ratelaws,133) -NV_Ith_S(ratelaws,141) +NV_Ith_S(ratelaws,153) +NV_Ith_S(ratelaws,159) +NV_Ith_S(ratelaws,163) -NV_Ith_S(ratelaws,165) +NV_Ith_S(ratelaws,167);
    NV_Ith_S(Dspecies,75) = NV_Ith_S(ratelaws,137) -NV_Ith_S(ratelaws,142);
    NV_Ith_S(Dspecies,76) = NV_Ith_S(ratelaws,140) -NV_Ith_S(ratelaws,166) -NV_Ith_S(ratelaws,168) -NV_Ith_S(ratelaws,178) +NV_Ith_S(ratelaws,179);
    NV_Ith_S(Dspecies,77) = NV_Ith_S(ratelaws,141) -NV_Ith_S(ratelaws,167) -NV_Ith_S(ratelaws,169) +NV_Ith_S(ratelaws,178) -NV_Ith_S(ratelaws,179);
    NV_Ith_S(Dspecies,78) = NV_Ith_S(ratelaws,143) -NV_Ith_S(ratelaws,170) -NV_Ith_S(ratelaws,171);
    NV_Ith_S(Dspecies,79) = NV_Ith_S(ratelaws,144) -NV_Ith_S(ratelaws,172) -NV_Ith_S(ratelaws,173) -NV_Ith_S(ratelaws,176);
    NV_Ith_S(Dspecies,80) = NV_Ith_S(ratelaws,147) -NV_Ith_S(ratelaws,174) -NV_Ith_S(ratelaws,175) -NV_Ith_S(ratelaws,177);


    return(0);
}


/*
**   ========
**   main MEX
**   ========
*/
void mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[] )
{
    /* variables */
    double *  return_status;
    double *  species_out;
    double *  observables_out;    
    double *  parameters;
    double *  species_init;
    double *  timepoints; 
    size_t    n_timepoints;
    size_t    i;
    size_t    j;

    /* intermediate data vectors */
    N_Vector  expressions;
    N_Vector  observables;
    N_Vector  ratelaws;

    /* array to hold pointers to data vectors */
    N_Vector  temp_data[3];
    
    /* CVODE specific variables */
    realtype  reltol;
    realtype  abstol;
    realtype  time;
    N_Vector  species;
    void *    cvode_mem;
    int       flag;

    /* check number of input/output arguments */
    if (nlhs != 3)
    {  mexErrMsgTxt("syntax: [err_flag, species_out, obsv_out] = network_mex( timepoints, species_init, params )");  }
    if (nrhs != 3)
    {  mexErrMsgTxt("syntax: [err_flag, species_out, obsv_out] = network_mex( timepoints, species_init, params )");  }


    /* make sure timepoints has correct dimensions */
    if ( (mxGetM(prhs[0]) < 2)  ||  (mxGetN(prhs[0]) != 1) )
    {  mexErrMsgTxt("TIMEPOINTS must be a column vector with 2 or more elements.");  }

    /* make sure species_init has correct dimensions */
    if ( (mxGetM(prhs[1]) != 1)  ||  (mxGetN(prhs[1]) != __N_SPECIES__) )
    {  mexErrMsgTxt("SPECIES_INIT must be a row vector with 81 elements.");  } 

    /* make sure params has correct dimensions */
    if ( (mxGetM(prhs[2]) != 1)  ||  (mxGetN(prhs[2]) != __N_PARAMETERS__) )
    {  mexErrMsgTxt("PARAMS must be a column vector with 120 elements.");  }

    /* get pointers to input arrays */
    timepoints   = mxGetPr(prhs[0]);
    species_init = mxGetPr(prhs[1]);
    parameters   = mxGetPr(prhs[2]);

    /* get number of timepoints */
    n_timepoints = mxGetM(prhs[0]);

    /* Create an mxArray for output trajectories */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL );
    plhs[1] = mxCreateDoubleMatrix(n_timepoints, __N_SPECIES__, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_timepoints, __N_OBSERVABLES__, mxREAL);

    /* get pointers to output arrays */
    return_status   = mxGetPr(plhs[0]);
    species_out     = mxGetPr(plhs[1]);
    observables_out = mxGetPr(plhs[2]);    
   
    /* initialize intermediate data vectors */
    expressions  = NULL;
    expressions = N_VNew_Serial(__N_EXPRESSIONS__);
    if (check_flag((void *)expressions, "N_VNew_Serial", 0))
    {
        return_status[0] = 1;
        return;
    }

    observables = NULL;
    observables = N_VNew_Serial(__N_OBSERVABLES__);
    if (check_flag((void *)observables, "N_VNew_Serial", 0))
    {
        N_VDestroy_Serial(expressions);
        return_status[0] = 1;
        return;
    }

    ratelaws    = NULL; 
    ratelaws = N_VNew_Serial(__N_RATELAWS__);
    if (check_flag((void *)ratelaws, "N_VNew_Serial", 0))
    {   
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);        
        return_status[0] = 1;
        return;
    }
    
    /* set up pointers to intermediate data vectors */
    temp_data[0] = expressions;
    temp_data[1] = observables;
    temp_data[2] = ratelaws;

    /* calculate expressions (expressions are constant, so only do this once!) */
    calc_expressions( expressions, parameters );

        
    /* SOLVE model equations! */
    species   = NULL;
    cvode_mem = NULL;

    /* Set the scalar relative tolerance */
    reltol = 1e-08;
    abstol = 1e-06;

    /* Create serial vector for Species */
    species = N_VNew_Serial(__N_SPECIES__);
    if (check_flag((void *)species, "N_VNew_Serial", 0))
    {  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);
        return_status[0] = 1;
        return;
    }
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   NV_Ith_S(species,i) = species_init[i];   }
    
    /* write initial species populations into species_out */
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   species_out[i*n_timepoints] = species_init[i];   }
    
    /* write initial observables populations into species_out */ 
    calc_observables( observables, species, expressions );  
    for ( i = 0; i < __N_OBSERVABLES__; i++ )
    {   observables_out[i*n_timepoints] = NV_Ith_S(observables,i);   }

    /*   Call CVodeCreate to create the solver memory:    
     *   CV_ADAMS or CV_BDF is the linear multistep method
     *   CV_FUNCTIONAL or CV_NEWTON is the nonlinear solver iteration
     *   A pointer to the integrator problem memory is returned and stored in cvode_mem.
     */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  



    /*   Call CVodeInit to initialize the integrator memory:     
     *   cvode_mem is the pointer to the integrator memory returned by CVodeCreate
     *   rhs_func  is the user's right hand side function in y'=f(t,y)
     *   T0        is the initial time
     *   y         is the initial dependent variable vector
     */
    flag = CVodeInit(cvode_mem, calc_species_deriv, timepoints[0], species);
    if (check_flag(&flag, "CVodeInit", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  
   
    /* Set scalar relative and absolute tolerances */
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                     
   
    /* pass params to rhs_func */
    flag = CVodeSetUserData(cvode_mem, &temp_data);
    if (check_flag(&flag, "CVodeSetFdata", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  
    
    /* select linear solver */
    flag = CVDense(cvode_mem, __N_SPECIES__);
    if (check_flag(&flag, "CVDense", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  
    
    flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  

    flag = CVodeSetMaxErrTestFails(cvode_mem, 7);
    if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  

    flag = CVodeSetMaxConvFails(cvode_mem, 10);
    if (check_flag(&flag, "CVodeSetMaxConvFails", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  

    flag = CVodeSetMaxStep(cvode_mem, 0.0);
    if (check_flag(&flag, "CVodeSetMaxStep", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        return_status[0] = 1;          
        return;                        
    }                                  

    /* integrate to each timepoint */
    for ( i=1;  i < n_timepoints;  i++ )
    {
        flag = CVode(cvode_mem, timepoints[i], species, &time, CV_NORMAL);
        if (check_flag(&flag, "CVode", 1))
        {
            N_VDestroy_Serial(expressions);
            N_VDestroy_Serial(observables);           
            N_VDestroy_Serial(ratelaws);
            N_VDestroy_Serial(species);
            CVodeFree(&cvode_mem);
            return_status[0] = 1; 
            return;
        }

        /* copy species output from nvector to matlab array */
        for ( j = 0; j < __N_SPECIES__; j++ )
        {   species_out[j*n_timepoints + i] = NV_Ith_S(species,j);   }
        
        /* copy observables output from nvector to matlab array */
        calc_observables( observables, species, expressions );         
        for ( j = 0; j < __N_OBSERVABLES__; j++ )
        {   observables_out[j*n_timepoints + i] = NV_Ith_S(observables,j);   }      
    }
 
    /* Free vectors */
    N_VDestroy_Serial(expressions);
    N_VDestroy_Serial(observables);  
    N_VDestroy_Serial(ratelaws);        
    N_VDestroy_Serial(species);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return;
}


/*  Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */
int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        mexPrintf( "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname );    
        return(1);
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *) flagvalue;
        if (*errflag < 0)
        {
            mexPrintf( "\nSUNDIALS_ERROR: %s() failed with flag = %d\n", funcname, *errflag );
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        mexPrintf( "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n", funcname );
        return(1);
    }

    return(0);
}
