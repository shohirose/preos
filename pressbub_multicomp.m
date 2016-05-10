function [pressb, comp_vap] = pressbub_multicomp(comp_liq, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_liq,1);

pressb = pressbest(comp_liq, temp, pressc, tempc, acentric);
K = wilsoneq(pressb, temp, pressc, tempc, acentric);

for loop = 1:maxiter
    
    [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP);
    
    K = Knew;
    pressb = pressbnew;
    
    eps = abs(f);
    
    if eps < tol
        break;
    end
    
end

if loop >= maxiter
    
    fprintf('The iteration in pressbub_multicomp() did not converge.\n');
    
end

comp_vap = zeros(ncomp,1);

for i = 1:ncomp
    
    comp_vap(i) = K(i)*comp_liq(i);
    
end

end

function pressb = pressbest(comp_liq, temp, pressc, tempc, acentric)

ncomp = size(comp_liq,1); 
pressb = 0;

for i = 1:ncomp
    
    pressb = pressb + comp_liq(i)*pressc(i)*...
        exp(5.373*(1 + acentric(i))*(1 - tempc(i)/temp));
    
end

end

function [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

comp_vap = zeros(ncomp,1);

for i = 1:ncomp
    
    comp_vap(i) = K(i)*comp_liq(i);
    
end

[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacityCoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP);

fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);

f = 0;
pressbnew = 0;

for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressb;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressb;
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    pressbnew = pressbnew + fug_liq(i)/fugcoef_vap(i);
    
end

f = comp_liq'*Knew - 1;

end