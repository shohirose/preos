function [K, comp_vap, comp_liq, phasefrac] = vaporliquideq(press, temp, comp_overall, pressc, tempc, acentric, BIP, tol, maxiter)

% initial estimate of equilibrium constant
K = wilsoneq(press, temp, pressc, tempc, acentric);

fun = @(x) objfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP);
grad = @(x) jacobfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP);

for i = 1:maxiter
    
    f    = fun(K);
    
    eps = max(abs(f));
    if eps < tol
        break;
    end
    
    dfdK = grad(K);
    dK   = -dfdK\f;
    K  = K + dK;
    
end

[phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
comp_vap = comp(:, 1);
comp_liq = comp(:, 2);

end

function f = objfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP)

ncomp = size(comp_overall, 1);
tol = 1e-5;
maxiter = 20;
[~, comp] = phasefraction(K, comp_overall, tol, maxiter);
comp_vap = comp(:, 1);
comp_liq = comp(:, 2);

[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, press, temp, pressc, tempc, acentric, BIP);
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, press, temp, pressc, tempc, acentric, BIP);

f = zeros(ncomp, 1);
for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end

%f(ncomp) = sum(comp_vap) - 1;

end

function J = jacobfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP)

ncomp = size(comp_overall, 1);

fun = @(x) objfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP);

f0 = fun(K);

perturb_K = 1e-6;
J = zeros(ncomp, ncomp);

for i = 1:ncomp
        
    dK = zeros(ncomp, 1);
    dK(i) = perturb_K;
    K1 = K + dK;
    
    f1 = fun(K1);
    J(:, i) = (f1 - f0)/perturb_K;
    
end

end