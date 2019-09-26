# EnKF_Analysis
EnKF analysis routines in Fortran 90.  Stochastich and SQRT formulations with subspace inversion. 

Analysis routines:

analysis.F90:   Routine updated 02/02/2007 including options for perturbed
                observations or symmetric square root schemes

analysis2_EnOI.F90: Ensemble OI version of the perturbed observations analysis

Note that if you using the square root scheme with the random rotation of anomalies AND the local analysis,
you have to make sure you use the same rotation for each grid point. You may then compute the random matrix
outside the analysis routine, and pass the same matrix in each local call.
