%% DEMO FOR PRESSBUB_PURECOMP

temp   = 230; % [K]
name = {'C2H6' 'C3H8' 'n-C4H10' 'n-C5H12' 'CO2'};
pressc   = [4.884, 4.246, 3.800, 3.374, 7.376]'*1e6; % [Pa]
tempc    = [305.4, 369.8, 425.2, 469.6, 304.2]'; % [K]
acentric = [0.098, 0.152, 0.193, 0.251, 0.225]';
ngas = size(pressc, 1);

tol = 1e-6;
maxiter = 20;

pressb_calc = [];
for i = 1:ngas;
    if (temp >= tempc)
        error('Gas %s has no bubble point at the temperature %f.\\n', name(i), temp);
    end
    pressb_int = pressbubest_purecomp(pressc(i), tempc(i), acentric(i), temp);
    pressb = pressbub_purecomp_newton(pressb_int, temp, pressc(i), tempc(i), acentric(i), tol, maxiter);
    pressb_calc = cat(1, pressb_calc, pressb);
end

figure;
plot([1:ngas], pressb_calc*1e-6, '-o');
%legend(char(name));
title('BUBBLE POINT PRESSURE AT T = 270 K');
xlabel('GAS COMPONENT');
ylabel('PRESSURE [MPa]');
ax = gca;
ax.XTick = [1:ngas];
ax.XTickLabel = char(name);