function pressb = pressbubest_multicomp(comp_liq, temp, pressc, tempc, acentric)

ncomp = size(comp_liq,1); 
pressb = 0;

for i = 1:ncomp
    pressb = pressb + comp_liq(i)*pressc(i)*...
        exp(5.373*(1 + acentric(i))*(1 - tempc(i)/temp));
end

end