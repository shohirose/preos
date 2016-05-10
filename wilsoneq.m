function K = wilsoneq(press, temp, pressc, tempc, acentric)

ncomp = size(pressc,1);   % the number of components
K  = zeros(ncomp,1);      % equilibrium constant

% estimate K by using Wilson equation.
for i = 1:ncomp
    
    K(i) = pressc(i)/press*exp(5.373*(1 + acentric(i))*(1 - tempc(i)/temp));
    
end

end