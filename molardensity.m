%% Calculate molar density

function rho = molardensity(press, temp, zfactor)

gconst = 8.3144598; % [J/(K-mol)]

rho = press/(zfactor*gconst*temp); % [mol/m3]

end