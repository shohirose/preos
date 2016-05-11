# eos
MatLab codes for the Peng-Robinson equation of state.

The main functions are as follows.

(1) Calculate z-factor and fugacity cofficients for pure component systems.
- fugcoef_purecomp.m : return the z-factor which minimizes the gibbs free energy if multiple roots are found.
- fugcoef_purecomp_liquid.m : return the minimum z-factor if multiple roots are found.
- fugcoef_purecomp_vapor.m: return the maximum z-factor if multiple roots are found.

(2) Calculate z-factor and fugacity cofficients for multi-component systems.
- fugcoef_multicomp.m : return the z-factor which minimizes the gibbs free energy if multiple roots are found.
- fugcoef_multicomp_liquid.m : return the minimum z-factor if multiple roots are found.
- fugcoef_multicomp_vapor.m: return the maximum z-factor if multiple roots are found.

(3) Calculate Bubble point pressure.
- pressbub_purecomp_ss.m : For pure component systems. Successive substitution method is used.
- pressbub_purecomp_newton.m : For pure component systems. Newton method is used.
- pressbub_multicomp.m : For multi-component systems. Successive substitution method is used.

(4) Calculate Dew point pressure.
- pressdew_multicomp.m : For multi-component systems. Successive substitution method is used.
