%% CALCULATE THE FUGACITY COEFFICIENT AND Z-FACTOR OF PURE COMPONENT
% -------------------------------------------------------------------------
% The Definition of Variables.
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% fugcoef : fugacity coefficient
% zfactor : compressibility factor
% -------------------------------------------------------------------------
% In this function, an appropriate z-factor is automatically chosen
% according to gibbs free energy if multiple roots are found.
function [fugcoef, zfactor] = fugacitycoef_purecomp(press, temp, pressc, tempc, acentric)

[A, B] = calcab_purecomp(press, temp, pressc, tempc, acentric);

zfactor = calczfactor(A, B);

if (size(zfactor,1) > 1)
    zfactor = choosezfactor(zfactor, A, B);
end

fugcoef = calcfugcoef_purecomp(zfactor, A, B);

end

%% CHOOSE AN APPROPRIATE Z-FACTOR
% Calculate dimensionless excess gibbs free energy, and return the z
% factor which minimizes the gibbs free energy.
function minzfactor = choosezfactor(zfactor, A, B)

gibbsenergy = [];
for i = 1:size(zfactor,1)
    fugcoef = calcfugcoef_purecomp(zfactor(i), A, B);
    g = log(fugcoef);
    gibbsenergy = cat(1,gibbsenergy,g);
end

[~, index] = sort(gibbsenergy);
minzfactor = zfactor(index(1));

end