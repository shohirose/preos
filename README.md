# eos
MatLab codes for the equation of state.

The main functions are as follows.

(1) Calculate z-factor and fugacity cofficients for pure component systems.
- fugcoef_purecomp.m : return the z-factor which minimizes the gibbs free energy if multiple roots are found and corresponding fugacity coefficient.
- fugcoef_purecomp_liquid.m : return the minimum z-factor if multiple roots are found and corresponding fugacity coefficient.
- fugcoef_purecomp_vapor.m: return the maximum z-factor if multiple roots are found and corresponding fugacity coefficient.

(2) Calculate z-factor and fugacity cofficients for multi-component systems.
- fugcoef_multicomp.m : return the z-factor which minimizes the gibbs free energy if multiple roots are found and corresponding fugacity coefficient.
- fugcoef_multicomp_liquid.m : return the minimum z-factor if multiple roots are found and corresponding fugacity coefficient.
- fugcoef_multicomp_vapor.m: return the maximum z-factor if multiple roots are found and corresponding fugacity coefficient.
