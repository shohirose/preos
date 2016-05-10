function pressb = pressbub_purecomp(temp, pressc, tempc, acentric, tol, maxiter)

pressb = pressbest(pressc, tempc, acentric, temp);

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

function pressb = pressbest(pressc, tempc, acentric, temp)

pressb = pressc*exp(5.373*(1 + acentric)*(1 - tempc/temp));

end