%% Phase mole fraction calculation
% The method used here is based on
% Okuno et al., 2010, A New Algorithm for Rachford-Rice For Multiphase
% Compositional Simulation, SPE Journal, 14 (2): 313-325

% Overall composition, 7 components.
comp_overall = [
    0.204322077
    0.070970999
    0.267194323
    0.296291965
    0.067046081
    0.062489248
    0.031685307 ];
% Equilibrium constant, 3 phases.
K = [
    1.234669887  1.527133414
    0.897277011  0.024564880
    2.295257081  1.463482405
    1.589548999  1.160905462
    0.233493486  0.241662899
    0.020381086  0.148152826
    1.407156410  14.31280108 ];

% tolerance and the maximu iteration.
tol = 1e-6;
maxiter = 20;

[phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);

disp("phasefrac is:"), disp(phasefrac);
disp("comp is:"), disp(comp);