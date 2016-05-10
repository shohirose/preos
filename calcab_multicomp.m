%% CALCULATE DIMENSIONLESS ATTRACTION AND COVOLUME, A & B
% The Definition of Variables.
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% ncomp   : the number of components

function [A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric)

ncomp = size(pressc,1);
m = zeros(ncomp,1);
alpha = zeros(ncomp,1);
A = zeros(ncomp,1);
B = zeros(ncomp,1);
omegaa = 0.45724;
omegab = 0.0778;

for i = 1:ncomp
    
    % Calculate m.
    if acentric(i) > 0.49
        m(i) = 0.379642 + 1.48503*acentric(i) - 0.164423*acentric(i)^2 + 0.016666*acentric(i)^3;
    else
        m(i) = 0.37464 + 1.54226*acentric(i) - 0.26992*acentric(i)^2;
    end
    
    % Calculate reduced pressure and temperature.
    pressr = press/pressc(i);
    tempr = temp/tempc(i);
    
    % Calculate alpha.
    alpha(i) = ( 1 + m(i)*(1 - sqrt(tempr)) )^2;
    
    A(i) = omegaa*alpha(i)*pressr/tempr^2;
    B(i) = omegab*pressr/tempr;
    
end

end
