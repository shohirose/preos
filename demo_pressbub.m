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
    if (temp >= tempc)
        error('Gas %s has no bubble point at the temperature %f.\\n', name(i), temp);
    end
    % Calculate the initial estimate of bubble point pressure.
    pressb_int = pressbubest_purecomp(pressc(i), tempc(i), acentric(i), temp);
    % Calculate bubble point pressure.
    pressb = pressbub_purecomp_newton(pressb_int, temp, pressc(i), tempc(i), acentric(i), tol, maxiter);
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
ax = gca;
ax.XTick = [1:ngas];
ax.XTickLabel = char(name);

%% DEMO FOR PRESSBUB_MULTICOMP
% Methane - Propane system
% comp1 = CH4, comp2 = C3H8
x = [0:0.01:0.67];
comp_liq = [x; 1 - x];
pressc1 = [pressc(1); pressc(3)];
tempc1  = [tempc(1); tempc(3)];
acentric1 = [acentric(1); acentric(3)];
BIP = [0, 0.09; 0.09, 0];

tol = 1e-8;
maxiter = 2000;

pressb_calc = [];
comp_vap = [];
fprintf('Calculating bubble and dew point pressure curves for CH4-C3H8 system...\n');
for i = 1:size(x,2);
    % Composition in liquid phase
    x0 = comp_liq(:,i);
    % Estimate initial bubble ponit pressure
    if i == 1
        press_ini = pressbubest_multicomp(x0, temp, pressc1, tempc1, acentric1);
    else
        press_ini = pressb;
    end
    % Calculate bubble point pressure
    [pressb, y] = pressbub_multicomp_newton(x0, press_ini, temp, pressc1, tempc1, acentric1, BIP, tol, maxiter);
    % Save the values.
    pressb_calc = cat(2, pressb_calc, pressb);
    comp_vap = cat(2, comp_vap, y);
end

% Plot the result.
figure;
plot(x, pressb_calc*1e-6, comp_vap(1,:), pressb_calc*1e-6);
title('BUBBLE POINT PRESSURE OF CH4-C3H8 SYSTEM AT T = 270 K');
xlabel('METHANE MOLE FRACTION');
ylabel('PRESSURE [MPa]');
legend('Bubble point pressure', 'Dew point pressure', 'Location','northwest');