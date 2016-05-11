%% CALCULATE BUBBLE POINT PRESSURE BY USING SUCCESSIVE SUBSTITUTION
function pressb = pressbub_purecomp_ss(pressb_ini, temp, pressc, tempc, acentric, tol, maxiter)

%pressb = pressbubest_purecomp(pressc, tempc, acentric, temp);
pressb = pressb_ini;

for loop = 1:maxiter
    
    [fugcoef_vap, ~] = fugacitycoef_purecomp_vapor(pressb, temp, pressc, tempc, acentric);
    [fugcoef_liq, ~] = fugacitycoef_purecomp_liquid(pressb, temp, pressc, tempc, acentric);
    
    fug_vap = pressb*fugcoef_vap;
    fug_liq = pressb*fugcoef_liq;
    
    eps = abs((fug_vap - fug_liq)/fug_vap);
    
    if eps < tol
        break;
    end
    
    pressb = pressb*fug_liq/fug_vap;
    
end

end