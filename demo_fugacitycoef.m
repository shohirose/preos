%% Demo for fugacitycoef_purecomp()

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
press = 3e6; % [Pa]
name = {'CH_4', 'C_2H_6' 'C_3H_8' 'n-C_4H_{10}' 'n-C_5H_{12}' 'CO_2'};
pressc   = [4.600, 4.884, 4.246, 3.800, 3.374, 7.376]'*1e6; % [Pa]
tempc    = [190.6, 305.4, 369.8, 425.2, 469.6, 304.2]'; % [K]
acentric = [0.008, 0.098, 0.152, 0.193, 0.251, 0.225]';
ngas = size(pressc, 1);

temp = [300:1:450]; % [K]

phi_all = [];
z_all = [];

for i = 1:ngas
    
    phi = [];
    z = [];
    
    for j = 1:size(temp, 2)
        [fugacoef, zfactor] = fugacitycoef_purecomp(press, temp(j), pressc(i), tempc(i), acentric(i));
        phi = cat(2, phi, fugacoef);
        z = cat(2, z, zfactor);
    end
    
    phi_all = cat(1, phi_all, phi);
    z_all = cat(1, z_all, z);
    
end

figure;
hold on;

subplot(1, 2, 1);
plot(temp, phi_all);
xlabel('Temperature, K');
ylabel('Fugacity Coefficient');
legend(char(name),'Location','eastoutside');

subplot(1, 2, 2);
plot(temp, z_all);
xlabel('Temperature, K');
ylabel('Z factor');
legend(char(name),'Location','eastoutside');

ax = gca;
axis([-Inf,Inf,0,Inf]);


%% Demo for fugacitycoef_multicomp()