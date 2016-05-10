function [pressd, comp_liq] = pressdew_multicomp(comp_vap, temp, pressc, tempc, acentric, BIP, tol, maxiter)

ncomp = size(comp_vap,1);

pressd = pressdest(comp_vap, temp, pressc, tempc, acentric);
K = wilsoneq(pressd, temp, pressc, tempc, acentric);

iter = 0;

for loop = 1:maxiter
    
    [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP);
    
    K = Knew;
    pressd = pressdnew;
    
    eps = abs(f);
    
    if eps < tol
        break;
    end
    
end

if loop >= maxiter
    
    fprintf('The iteration in pressbub_multicomp() did not converge.\n');
    
end

comp_liq = zeros(ncomp,1);

for i = 1:ncomp
    
    comp_liq(i) = comp_vap(i)/K(i);
    
end

end

function pressd = pressdest(comp_vap, temp, pressc, tempc, acentric)

ncomp = size(comp_vap,1); 
pressd = 0;

for i = 1:ncomp
    
    pressd = pressd + comp_vap(i)/...
        (pressc(i)*exp(5.373*(1 + acentric(i))*(1 - tempc(i)/temp)));
    
end

pressd = 1/pressd;

end

function [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP)

ncomp = size(K,1);

comp_liq = zeros(ncomp,1);

for i = 1:ncomp
    
    comp_liq(i) = comp_vap(i)/K(i);
    
end

[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP);

fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);

f = 0;
pressdnew = 0;

for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressd;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressd;
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    pressdnew = pressdnew + fug_vap(i)/fugcoef_liq(i);
    f = f + comp_vap(i)/Knew(i);
    
end

f = f - 1;

end