function [pressd, comp_liq] = pressdew_multicomp(comp_vap, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_vap,1);

% Input initial values.
pressd = pressdewest_multicomp(comp_vap, temp, pressc, tempc, acentric);
K = wilsoneq(pressd, temp, pressc, tempc, acentric);
%K = updatek(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP);

for loop = 1:maxiter
    
    [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP);
    
    % Update K and dew point pressure.
    K = Knew;
    pressd = pressdnew;
    
    % Check convergence.
    eps = abs(f);
    if eps < tol
        break;
    end
    
end

% Echo a message if the iteration did not converge.
if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp() did not converge. eps = %E\n', eps);
else
    fprintf('Iteration = %d\n', loop);
end

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

function [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

% Calculate fuacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP);

% Calculate new K and dew point pressure.
fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);
f = 0;
pressdnew = 0;

for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressd;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressd;
    % Calculate new K.
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    % Calculate new dew point pressure.
    pressdnew = pressdnew + fug_vap(i)/fugcoef_liq(i);
    % Calculate the objective funciton.
    f = f + comp_vap(i)/Knew(i);
    
end

f = f - 1;

end