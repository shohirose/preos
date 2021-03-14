%% Demo for pressbub_purecomp_newton()

% The Definition of Variables
% temp : temperature [K]
% name : component's name
% pressc : critical pressure
% tempc  : critical temperature
% acentric : acentric factor
% ngas : the number of gas components
% tol : toleratnce for convergence
% maxiter : the maximum iteration

% Input data
temp   = 270; % [K]
name = {'CH4', 'C2H6' 'C3H8' 'n-C4H10' 'n-C5H12' 'CO2'};
pressc   = [4.600, 4.884, 4.246, 3.800, 3.374, 7.376]'*1e6; % [Pa]
tempc    = [190.6, 305.4, 369.8, 425.2, 469.6, 304.2]'; % [K]
acentric = [0.008, 0.098, 0.152, 0.193, 0.251, 0.225]';
ngas = size(pressc, 1);

tol = 1e-8;
maxiter = 20;

% Calculate bubble point pressure for each component.
fprintf('Calculating vapor pressure of pure component systems...\n');
pressb_calc = [];
for i = 1:ngas;
    if (temp >= tempc(i))
        fprintf('Gas %s has no bubble point at the temperature %4.2f.\\n', char(name(i)), temp);
        pressb = 0;
    else
        % Calculate the initial estimate of bubble point pressure.
        pressb_int = pressbubest_purecomp(pressc(i), tempc(i), acentric(i), temp);
        % Calculate bubble point pressure.
        pressb = pressbub_purecomp_newton(pressb_int, temp, pressc(i), tempc(i), acentric(i), tol, maxiter);
    end
    % Store the calculated value.
    pressb_calc = cat(1, pressb_calc, pressb);
end

% Plot the calculation results.
figure;
plot([1:ngas], pressb_calc*1e-6, '-o');
%legend(char(name));
title('VAPOR PRESSURE AT T = 270 K');
xlabel('GAS COMPONENT');
ylabel('PRESSURE [MPa]');
%ax = gca;
%ax.XTick = [1:ngas];
%ax.XTickLabel = char(name);

%% Demo for pressbub_multicomp_newton() and pressdew_multicomp_newton()
% Methane - Propane system
% comp1 = CH4, comp2 = C3H8
x = [0:0.01:0.65];
comp = [x; 1 - x];
pressc1 = [pressc(1); pressc(3)];
tempc1  = [tempc(1); tempc(3)];
acentric1 = [acentric(1); acentric(3)];
BIP = [0, 0.09; 0.09, 0];

tol = 1e-8;
maxiter = 50;

ndata = size(x, 2);
comp1 = zeros(ndata, 4);
press_calc = zeros(ndata, 4);

fprintf('Calculating bubble and dew point pressure curves for CH4-C3H8 system...\n');

for i = 1:ndata
    
    % Composition in liquid phase
    x0 = comp(:,i);
    
    % Estimate initial bubble ponit pressure
    if i == 1
        pressb_ini = pressbubest_multicomp(x0, temp, pressc1, tempc1, acentric1);
        pressd_ini = pressdewest_multicomp(x0, temp, pressc1, tempc1, acentric1);
    else
        pressb_ini = pressb;
        pressd_ini = pressd;
    end
    
    % Calculate bubble and dew point pressure
    [pressb, yb] = pressbub_multicomp_newton(x0, pressb_ini, temp, pressc1, tempc1, acentric1, BIP, tol, maxiter);
    [pressd, yd] = pressdew_multicomp_newton(x0, pressd_ini, temp, pressc1, tempc1, acentric1, BIP, tol, maxiter);
    
    % Save the values.
    press_calc(i, 1) = pressb;
    press_calc(i, 2) = pressb;
    press_calc(i, 3) = pressd;
    press_calc(i, 4) = pressd;
    
    comp1(i, 1) = x0(1);
    comp1(i, 2) = yb(1);
    comp1(i, 3) = yd(1);
    comp1(i, 4) = x0(1);
    
end

% Plot the result.
figure;
subplot(1,2,1);
plot(comp1(:,1:2), press_calc(:,1:2)*1e-6);
title('CH_4-C_3H_8 SYSTEM AT T = 270 K');
xlabel('Methane Mole Composition');
ylabel('Pressure [MPa]');
legend('Pb (pressbub)', 'Pd (pressbub)', 'Location','northwest');
axis([0,1,0,12]);

subplot(1,2,2);
plot(comp1(:,3:4), press_calc(:,3:4)*1e-6);
title('CH_4-C_3H_8 SYSTEM AT T = 270 K');
xlabel('Methane Mole Composition');
ylabel('Pressure [MPa]');
legend('Pb (pressdew)', 'Pd (pressdew)', 'Location','northwest');
axis([0,1,0,12]);