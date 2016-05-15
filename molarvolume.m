%% Calculate molar volume

function v = molarvolume(press, temp, zfactor)

gconst = 8.3144598; % [J/(K-mol)]

v = zfactor*gconst*temp/press; % [m3/mol]

end