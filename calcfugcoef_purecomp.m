%% CALCULATE THE FUGACITY COEFFICIENT OF PURE COMPONENT
function fugcoef = calcfugcoef_purecomp(zfactor, A, B)

c1 = 1 + sqrt(2);
c2 = 1 - sqrt(2);

if (zfactor < B)
    error('Z factor must be larger than b in PR-EOS.\n');
end

fugcoef = (zfactor - 1) - log(zfactor - B) ...
    - A/(2*sqrt(2)*B)*log((zfactor + c1*B)/(zfactor + c2*B));

fugcoef = exp(fugcoef);

end