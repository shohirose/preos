%% Tangend Plane Distance analysis by using successive substitution
%
% twophase : 
% comp_overall : Overall composition
% press        : Pressure
% temp         : Temperature
% pressc       : Critical pressure
% tempc        : Critical temperature
% accentric    : Accentric factor
% BIP          : Binary interaction parameters
function twophase = tpdss(comp_overall, press, temp, pressc, tempc, accentric, BIP)

% The initial estimate of equilibrium constant
K = wilsoneq(press, temp, pressc, tempc, accentric);

% Vapor-like phase
comp_vap = compvapest(K, comp_overall);
vaporphase = checkphasesplit(comp_vap, comp_overall, press, temp, pressc, tempc, accentric, BIP);

% Liquid-like phase
comp_liq = compliqest(K, comp_overall);
liquidphase = checkphasesplit(comp_liq, comp_overall, press, temp, pressc, tempc, accentric, BIP);

if (vaporphase == true) || (liquidphase == true)
    twophase = true;
else
    twophase = false;
end
end

function phasesplit = checkphasesplit(comp_check, comp_overall, press, temp, pressc, tempc, accentric, BIP)

[fugcoef_overall, ~] = fugacitycoef_multicomp(comp_overall, press, temp, pressc, tempc, accentric, BIP);

ncomp = size(comp_overall, 1);
x = comp_check;
tol = 1e-6;
eps = 1;
maxiter = 20;

for iter = 1:maxiter
    
    xsum = sum(x);
    x = x/xsum;
    
    if max(abs(x - comp_overall)) < tol
        phasesplit = false;
        return;
    end
    
    if eps < tol
        break;
    end
    
    [fugcoef_x, ~] = fugacitycoef_multicomp(x, press, temp, pressc, tempc, accentric, BIP);
    eps = calcrsd(x, fugcoef_x, comp_overall, fugcoef_overall);
    x = updatex(comp_overall, fugcoef_overall, fugcoef_x);
    
end

if iter >= maxiter
    fprintf('The iteration in checkphasesplit() did not converge. One phase is assumed.\n');
%     fprintf('P = %1.3e, T = %4.2f, ', press, temp);
%     % plot comp_check
%     fprintf('comp_check = [ ');
%     for i = 1:ncomp
%         fprintf('%1.3f ', comp_check(i));
%     end
%     fprintf('], ');
%     % plot x
%     fprintf('x = [ ');
%     for i = 1:ncomp
%         fprintf('%1.3f ', x(i));
%     end
%     fprintf(']\n');
end

if xsum > 1
    phasesplit = true;
else
    phasesplit = false;
end

end

function comp_vap = compvapest(K, comp_overall)

ncomp = size(comp_overall,1);

comp_vap = zeros(ncomp,1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_overall(i);
end

end

function comp_liq = compliqest(K, comp_overall)

ncomp = size(comp_overall,1);

comp_liq = zeros(ncomp,1);
for i = 1:ncomp
    comp_liq(i) = comp_overall(i)/K(i);
end

end

function rsd = calcrsd(x, fugcoef_x, comp_overall, fugcoef_overall)

ncomp = size(x,1);

rsd = 0;
for i = 1:ncomp
    rsd = rsd + ( log(x(i)/comp_overall(i)*fugcoef_x(i)/fugcoef_overall(i)) )^2;
end

end

function x = updatex(comp_overall, fugcoef_overall, fugcoef_x)

ncomp = size(comp_overall,1);

x = zeros(ncomp,1);
for i = 1:ncomp
    x(i) = comp_overall(i)*fugcoef_overall(i)/fugcoef_x(i);
end

end