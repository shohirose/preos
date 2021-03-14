%% Computes dew point pressure in a multi-component system
%
% pressd     : Dewpoint pressure
% comp_liq   : Liquid composition
% comp_vap   : Vapor composition
% pressd_ini : Initial estimate of vapor pressure
% temp       : Temperature
% pressc     : Critical pressure
% tempc      : Critical temperature
% acentric   : Accentric factor
% BIP        : Binary interaction parameters
% tol        : Iteraction tolerance
% maxiter    : Maximum number of iterations
function [pressd, comp_liq] = pressdew_multicomp_newton(comp_vap, pressd_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_vap,1);

% Input initial values.
%pressd = pressdewest_multicomp(comp_vap, temp, pressc, tempc, acentric);
pressd = pressd_ini;
K = wilsoneq(pressd, temp, pressc, tempc, acentric);
K = updatek(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP);

fun = @(x) objfun(x(1:ncomp), comp_vap, x(ncomp + 1), temp, pressc, tempc, acentric, BIP);

% model parameters.
m = [K; pressd];

for loop = 1:maxiter
    
    f = fun(m);
    J = jacobfun(f, m, fun);
    dm = -J\f;
    
    % Update m.
    m = m + dm;
    
    % Check convergence.
    eps = abs(f);
    if eps < tol
        break;
    end
    
end

% Echo a message if the iteration did not converge.
if loop >= maxiter
    fprintf('The iteration in pressdew_multicomp_newton() did not converge: eps = %e\n', eps);
else
    fprintf('iter = %d, objfun = [ ', loop);
    for i = 1:ncomp+1
        fprintf('%1.3e ', f(i));
    end
    fprintf(']\n');
end

K = m(1:ncomp);
pressd = m(ncomp + 1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

end

function comp_liq = calccompliq(K, comp_vap)

ncomp = size(comp_vap, 1);

comp_liq = zeros(ncomp, 1);
for i = 1:ncomp
    comp_liq(i) = comp_vap(i)/K(i);
end

end

function Knew = updatek(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP)

ncomp = size(comp_vap, 1);
Knew = zeros(ncomp, 1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

% Calculate fuacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP);

for i = 1:ncomp
   Knew(i) = fugcoef_liq(i)/fugcoef_vap(i); 
end

end

function f = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

% Calculate fuacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP);

% Calculate objective function.
f = zeros(ncomp + 1, 1);

for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end
f(ncomp + 1) = sum(comp_liq) - 1;

end

%% Jacobian matrix
% $J_{ij} = \frac{\partial f_i}{\partial x_j}$
function J = jacobfun(f0, x, fun)
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