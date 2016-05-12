function [pressb, comp_vap] = pressbub_multicomp_newton(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_liq,1);

% Input initial values.
%pressb = pressbubest_multicomp(comp_liq, temp, pressc, tempc, acentric);
pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);

% m : model parameters
fun = @(m) objfun(m(1:ncomp,1), comp_liq, m(ncomp + 1, 1), temp, pressc, tempc, acentric, BIP);

x = [K; pressb];

for loop = 1:maxiter
    
    % Calculate update direction by using Newton-Raphson method.
    f = fun(x);
    dfdx = 0;
    dx = 0;
    if loop <= 2
        % update only K by using successive substitution.
        dx = -f;
        dx(ncomp + 1) = 0;
    else
        dfdx = grad(f, x, fun);
        %dx = bicg(dfdx, -f, 1e-10, 20);
        dx = -dfdx\f;
    end
    
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
    fprintf('The iteration in pressbub_multicomp() did not converge. eps = %E\n', eps);
else
    fprintf('Iteration = %d, Objective Function f = [ ', loop);
    for i = 1:ncomp+1
        fprintf('%1.3e ', f(i));
    end
    fprintf(']\n');
end

% Update K and pressb.
K = x(1:ncomp, 1);
pressb = x(ncomp + 1, 1);

% Calculate vapor composition.
comp_vap = zeros(ncomp,1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_liq(i);
end

end

function f = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = zeros(ncomp,1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_liq(i);
end

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP);

% Calculate the objective function.
% $F_i = K_i - \hat{\phi}_i^L/\hat{\phi}_i^V,  i = 1,\dots,N_c$
% $F_{N_c + 1} = \sum_i z_i K_i - 1$
f = zeros(ncomp + 1, 1);
for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end
f(ncomp + 1) = sum(comp_vap) - 1;

end

function dfdx = grad(f0, x, fun)
N = size(x, 1);
dfdx = zeros(N, N);
perturb_x = 1e-6;
for i = 1:N
    dx = zeros(N, 1);
    dx(i) = perturb_x;
    f1 = fun(x + dx);
    dfdx(:, i) = (f1 - f0)/perturb_x;
end
end