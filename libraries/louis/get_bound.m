clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN          = 10000;
N_values    = round(linspace(5, 1000, NN));
a           = linspace(0, 1, NN);
a_opt       = zeros(size(N_values));

for k = 1:length(N_values)
    N           = N_values(k);
    Bound       = get_Bound(N,a);
    E           = abs(Bound - 0.5);
    [~, i_opt]  = min(E);    
    a_opt(k)    = a(i_opt);
end

% figure
% loglog(N_values,abs(a_opt-0.5))
% xlabel("number of symbols that interfere")
% ylabel("intergral parameter - 1/2")
% set(gca,'fontsize',20)
% grid on

x = N_values;
y = a_opt-0.5;

[fitResult,gof] = fit_rat(x, y);

step = 50;
M = [x(1:step:end)',y(1:step:end)'];
T = array2table(M, 'VariableNames', {'N','a'});
fname = sprintf(sprintf('upper_bound_a_%d.txt',step));
writetable(T, fname, 'Delimiter', ',')

NN          = 100;
N_values    = round(linspace(5, 1000, NN));
A           = 0.5041+0.82./N_values;


Bounds = get_Bound(N_values,A);
figure
semilogx(N_values,Bounds-0.5)
% ylim([0,1])

step = 50;
M = [N_values',Bounds'];
T = array2table(M, 'VariableNames', {'N','B'});
fname = sprintf(sprintf('upper_bound_B_%d.txt',step));
writetable(T, fname, 'Delimiter', ',')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = get_Bound(N,a)
    B = sqrt(pi) * (erf((N - a)/2) - erf(a/2)) - (1 - 1./N) * exp(-1/4);
end

function [fitResult,gof] = fit_pow(x, y)
    fitType             = 'power1';
    fitOptions          = fitoptions(fitType);
    [fitResult, gof]    = fit(x', y', fitType, fitOptions);
end

function [fitResult, gof] = fit_rat(x, y)
    fitType             = 'rat11';           % Modèle: (p1*x + p2) / (x + q1)
    fitOptions          = fitoptions(fitType);
    [fitResult, gof]    = fit(x', y', fitType, fitOptions);
end