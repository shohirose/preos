function [pressb, comp_vap] = pressbub_multicomp_ss(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_liq,1);

% Input initial values.
%pressb = pressbubest_multicomp(comp_liq, temp, pressc, tempc, acentric);
pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);

for loop = 1:maxiter
    
    [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP);
    
    % Update K and bubble point pressure.
    K = Knew;
    pressb = pressbnew;
    
    % Check convergence.
    eps = abs(f);
    if eps < tol
        break;
    end
    
end

% Echo a message if the loop did not converge.
if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp() did not converge. eps = %E\n', eps);
else
    fprintf('Iteration = %d\n', loop);
end

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

end

function comp_vap = calccompvap(K, comp_liq)

ncomp = size(comp_liq, 1);

comp_vap = zeros(ncomp, 1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_liq(i);
end

end

function [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP);

% Calculate new K and bubble point pressure.
fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);
pressbnew = 0;
for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressb;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressb;
    % Calculate new K by successive substitution method.
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    % Calculate new bubble point pressure.
    pressbnew = pressbnew + fug_liq(i)/fugcoef_vap(i);
    
end

% Calculate objective function f.
f = comp_liq'*Knew - 1;

end