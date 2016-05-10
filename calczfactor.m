function zfactor = calczfactor(a, b)

% Calculate the coefficients of cubic equation.
a1 = 1;
a2 = b - 1;
a3 = a - 2*b - 3*b^2;
a4 = -a*b + b^2 + b^3;

% Solve the cubic equation.
zroots = roots([a1 a2 a3 a4]);

% Choose the real roots.
zfactor = [];
for i = 1:3
    if imag( zroots(i) ) == 0
        z = real(zroots(i));
        zfactor = cat(1, zfactor, z);
    end
end

zfactor = sort(zfactor);
end
