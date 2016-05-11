%% CALCULATE BUBBLE POINT PRESSURE BY USING SUCCESSIVE SUBSTITUTION
function pressb = pressbub_purecomp_newton(pressb_ini, temp, pressc, tempc, acentric, tol, maxiter)

%pressb = pressbubest_purecomp(pressc, tempc, acentric, temp);
pressb = pressb_ini;

fugcoef_vap = @(press) fugacitycoef_purecomp_vapor(press, temp, pressc, tempc, acentric);
fugcoef_liq = @(press) fugacitycoef_purecomp_liquid(press, temp, pressc, tempc, acentric);
objfun = @(press) fugcoef_vap(press) - fugcoef_liq(press);

perturbp = 1e-3;

for loop = 1:maxiter
    
    % Calculate an update direction by using Newton method.
    f0 = objfun(pressb);
    f1 = objfun(pressb + perturbp);
    dfdp = (f1 - f0)/perturbp;
    dp = -f0/dfdp;
    
    % Update bubble point prssure.
    pressb = pressb + dp;
    
    eps = abs(f0);
    
    % Check the convergence.
    if eps < tol
        break;
    end
    
end

if (loop >= maxiter)
    fprintf('Iterations reached at the maximum limit in pressbub_purecomp_newton().');
end

end

