function pressd = pressdewest_multicomp(comp_vap, temp, pressc, tempc, acentric)

ncomp = size(comp_vap,1); 
pressd = 0;
for i = 1:ncomp
    pressd = pressd + comp_vap(i)/...
        (pressc(i)*exp(5.373*(1 + acentric(i))*(1 - tempc(i)/temp)));
end
pressd = 1/pressd;

end