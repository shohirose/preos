function [Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP)

ncomp = size(comp,1);
Aij = zeros(ncomp, ncomp);

for i = 1:ncomp

    for j = 1:ncomp
        
        Aij(i,j) = sqrt(A(i)*A(j))*(1 - BIP(i,j));
        
    end
    
end

% Amix and Bmix are scalars.
Amix = comp'*Aij*comp;
Bmix = comp'*B;
% Amix2 is a vector.
Amix2 = Aij*comp;

end