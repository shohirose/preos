%% CALCULATE DIMENSIONLESS ATTRACTION AND COVOLUME, A & B
% The Definition of Variables.
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
function [A, B] = calcab_purecomp(press, temp, pressc, tempc, acentric)

omegaa = 0.45724;
omegab = 0.0778;

pressr = press/pressc;  % reduced pressure
tempr = temp/tempc;     % reduced temperature

if acentric > 0.49
    m = 0.379642 + 1.48503*acentric - 0.164423*acentric^2 + 0.016666*acentric^3;
else
    m = 0.37464 + 1.54226*acentric - 0.26992*acentric^2;
end

alpha = ( 1 + m*(1 - sqrt(tempr)) )^2;

A = omegaa*alpha*pressr/tempr^2;
B = omegab*pressr/tempr;

end