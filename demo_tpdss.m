%% Demo for tpdss function (2)

% Input data
temp   = 270; % [K]
name = {'CH4' 'C3H8'};
pressc   = [4.600, 4.246]'*1e6; % [Pa]
tempc    = [190.6, 369.8]'; % [K]
acentric = [0.008, 0.152]';
BIP = [
    0    0.09
    0.09 0    ];

ncomp = 2;
tol = 1e-6;
maxiter = 20;
comp_overall = [0.4, 0.6]';

% The low and high limit for calculation
press_low = 1e5; % 0.5 MPa
press_high = 1e7; % 10 MPa

% Calculate dew point pressure for phase identification.
pressdew_ini = pressdewest_multicomp(comp_overall, temp, pressc, tempc, acentric);
pressdew = pressdew_multicomp_newton(comp_overall, pressdew_ini, temp, pressc, tempc, acentric, BIP, 1e-6, 30);

% The number of calculation points
N = 101;

pdata = zeros(N,1);
vdata = zeros(N,2);

for i = 1:N
    
    press = press_low + (press_high - press_low)*(i - 1)/(N - 1);
    pdata(i,1) = press;
    
    % Check if phase will split.
    phasesplit = tpdss(comp_overall, press, temp, pressc, tempc, acentric, BIP);
    
    if phasesplit == true
        % For two phase.
        % Calculate phase composition and phase mole fraction.
        [~, comp_vap, comp_liq, phasefrac] = vaporliquideq(press, temp, comp_overall, pressc, tempc, acentric, BIP, tol, maxiter);
        % Calculate z-factor.
        [~, zfactor_vap] = fugacitycoef_multicomp_vapor(comp_vap, press, temp, pressc, tempc, acentric, BIP);
        [~, zfactor_liq] = fugacitycoef_multicomp_liquid(comp_liq, press, temp, pressc, tempc, acentric, BIP);
        % Calculate molar volume.
        v_vap = molarvolume(press, temp, zfactor_vap);  % vapor phase
        v_liq = molarvolume(press, temp, zfactor_liq);  % liquid phase
        v = v_vap*phasefrac + v_liq*(1 - phasefrac);    % overall
        vdata(i,2) = phasefrac;
    else
        % For one phase.
        % Calculate z-factor.
        [~, zfactor] = fugacitycoef_multicomp(comp_overall, press, temp, pressc, tempc, acentric, BIP);
        % Calculate molar volume.
        v = molarvolume(press, temp, zfactor);
        % Identify the phase by using dew point pressure.
        if press < pressdew
            vdata(i,2) = 1;
        end
    end
    vdata(i,1) = v;
    
end

figure;
subplot(1,2,1);
plot(vdata(:,1), pdata*1e-6);
ylabel('Pressure [MPa]');
xlabel('Molar Volume [m3/mol]');
%ax = gca;
%ax.XScale = 'log';

subplot(1,2,2);
plot(vdata(:,2), pdata*1e-6);
ylabel('Pressure [MPa]');
xlabel('Vapor Phase Mole Fraction');
axis([0,1,0,Inf]);
