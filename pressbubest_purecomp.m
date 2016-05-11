function pressb = pressbubest_purecomp(pressc, tempc, acentric, temp)

pressb = pressc*exp(5.373*(1 + acentric)*(1 - tempc/temp));

end