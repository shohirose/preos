function pressd = pressdew_purecomp(temp, pressc, tempc, acentric, tol, maxiter)

pressd = pressdest(pressc, tempc, acentric, temp);

for loop = 1:maxiter
    
    [fugcoef_vap, ~] = fugacitycoef_purecomp_vapor(pressd, temp, pressc, tempc, acentric);
    [fugcoef_liq, ~] = fugacitycoef_purecomp_liquid(pressd, temp, pressc, tempc, acentric);
    
    fug_vap = pressd*fugcoef_vap;
    fug_liq = pressd*fugcoef_liq;
    
    eps = abs((fug_vap - fug_liq)/fug_vap);
    
    if eps < tol
        break;
    end
    
    pressd = pressd*fug_liq/fug_vap;
    
end

end

function pressb = pressdest(pressc, tempc, acentric, temp)

pressb = pressc*exp(5.373*(1 + acentric)*(1 - tempc/temp));

end