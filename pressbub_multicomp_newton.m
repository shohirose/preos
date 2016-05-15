function [pressb, comp_vap] = pressbub_multicomp_newton(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_liq,1);

% Input initial values.
%pressb = pressbubest_multicomp(comp_liq, temp, pressc, tempc, acentric);
pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);
K = updatekss(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP);

% m : model parameters
fun = @(m) objfun(m(1:ncomp, 1), comp_liq, m(ncomp + 1, 1), temp, pressc, tempc, acentric, BIP);
updatek = @(m) updatekss(m(1:ncomp, 1), comp_liq, m(ncomp + 1, 1), temp, pressc, tempc, acentric, BIP);
x = [K; pressb];

for loop = 1:maxiter
    
    % Calculate update direction by using Newton-Raphson method.
    f = fun(x);
    J = jacob(f, x, fun);
    dx = -J\f;
    
    % Update x.
    x = x + dx;
    
    % Check convergence.
    eps = abs(max(f));
    if eps < tol
        break;
    end
    
end

% Echo a message if the loop did not converge.
if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp_newton() did not converge: eps = %e\n', eps);
else
    fprintf('iter = %d, objfun = [ ', loop);
    for i = 1:ncomp+1
        fprintf('%1.3e ', f(i));
    end
    fprintf(']\n');
end

% Update K and pressb.
K = x(1:ncomp, 1);
pressb = x(ncomp + 1, 1);

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

%% Update K by using successive substitution
function K = updatekss(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP);

% Calculate equilibrium constant.
K = zeros(ncomp, 1);
for i = 1:ncomp
    K(i) = fugcoef_liq(i)/fugcoef_vap(i);
end

end

%% Objective function
% $f_i = K_i - \hat{\phi}_i^L/\hat{\phi}_i^V,  i = 1,\dots,N_c$
% $f_{N_c + 1} = \sum_i z_i K_i - 1$
function f = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP);

% Calculate the objective function.
f = zeros(ncomp + 1, 1);
for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end
f(ncomp + 1) = sum(comp_vap) - 1;

end

%% Jacobian matrix
% $J_{ij} = \frac{\partial f_i}{\partial x_j}$
function J = jacob(f0, x, fun)
N = size(x, 1);
J = zeros(N, N);
perturb_x = 1e-6;
for i = 1:N
    dx = zeros(N, 1);
    dx(i) = perturb_x;
    f1 = fun(x + dx);
    J(:, i) = (f1 - f0)/perturb_x;
end
end