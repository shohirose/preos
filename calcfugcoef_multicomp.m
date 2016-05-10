function fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2)

ncomp = size(A,1);
fugcoef = zeros(ncomp, 1);

c0 = 2*sqrt(2);
c1 = 1 + sqrt(2);
c2 = 1 - sqrt(2);

for i = 1:ncomp
    
    if zfactor < Bmix
        error('Z-factor must be larger than Bmix.\n');
    end
    
    fugcoef(i) = B(i)/Bmix*(zfactor - 1) - log(zfactor - Bmix) ...
        - Amix/(c0*Bmix) * (2*Amix2(i)/Amix - B(i)/Bmix) ...
        *log((zfactor + c1*Bmix)/(zfactor + c2*Bmix));
    
    fugcoef(i) = exp(fugcoef(i));
    
end

end