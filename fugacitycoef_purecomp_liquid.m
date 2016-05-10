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
% In this function, the minimum z-factor is chosen if multiple roots are
% found.
function [fugcoef, zfactor] = fugacitycoef_purecomp_liquid(press, temp, pressc, tempc, acentric)

[A, B] = calcab_purecomp(press, temp, pressc, tempc, acentric);

zfactor = calczfactor(A, B);

if (size(zfactor,1) > 1)
    zfactor = min(zfactor);
end

fugcoef = calcfugcoef_purecomp(zfactor, A, B);

end