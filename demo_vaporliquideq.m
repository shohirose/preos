%% Demo for vaporliquideq function

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

tol = 1e-8;
maxiter = 50;

comp_overall = [0.4, 0.6]';

pressbub_ini = pressbubest_multicomp(comp_overall, temp, pressc, tempc, acentric);
pressdew_ini = pressdewest_multicomp(comp_overall, temp, pressc, tempc, acentric);
[pressbub, ~] = pressbub_multicomp_newton(comp_overall, pressbub_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter);
[pressdew, ~] = pressdew_multicomp_newton(comp_overall, pressdew_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter);

N = 41;
data = zeros(N, 4);

for i = 1:N
    press = pressdew + (pressbub - pressdew)*(i - 1)/(N - 1);
    [K, comp_vap, comp_liq, phasefrac] = vaporliquideq(press, temp, comp_overall, pressc, tempc, acentric, BIP, tol, maxiter);
    data(i, 1) = press;
    data(i, 2) = comp_vap(1);
    data(i, 3) = comp_liq(1);
    data(i, 4) = phasefrac;
end

figure;
plot(data(:,1)*1e-6, data(:,2:4));
xlabel('Pressure [MPa]');
ylabel('n_V, x_{CH_4}^V, x_{CH_4}^L');
axis([-Inf,Inf,0,1]);
legend('x_{CH_4}^V', 'x_{CH_4}^L', 'n_V','Orientation', 'horizontal');

%% Demo for tpdss function

press = pressdew*0.8;
phasesplit = tpdss(comp_overall, press, temp, pressc, tempc, acentric, BIP);
if phasesplit == true
    fprintf('V and L phase exist at P = %1.3f MPa\n', press*1e-6);
else
    fprintf('Either V or L phase coexists at P = %1.3f MPa\n', press*1e-6);
end

press = (pressdew + pressbub)*0.5;
phasesplit = tpdss(comp_overall, press, temp, pressc, tempc, acentric, BIP);
if phasesplit == true
    fprintf('V and L phase exist at P = %1.3f MPa\n', press*1e-6);
else
    fprintf('Either V or L phase coexists at P = %1.3f MPa\n', press*1e-6);
end

press = pressbub*1.2;
phasesplit = tpdss(comp_overall, press, temp, pressc, tempc, acentric, BIP);
if phasesplit == true
    fprintf('V and L phase exist at P = %1.3f MPa\n', press*1e-6);
else
    fprintf('Either V or L phase coexists at P = %1.3f MPa\n', press*1e-6);
end