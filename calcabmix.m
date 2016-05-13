%% Calculate A and B for mixtures by using the classical mixing rule.

% $A_{mix} = \sum_{ij} x_i x_j A_{ij}$
% $A_{ij} = \sqrt{A_i A_j}(1-k_{ij})$
% $B_{mix} = \sum_{i}x_i B_i$
% $A_{mix2,j} = \sum_{i} A_{ij} x_i$

% $k_{ij}$ is binary interaction parameter.

function [Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP)

ncomp = size(comp,1);
% Calculate 
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