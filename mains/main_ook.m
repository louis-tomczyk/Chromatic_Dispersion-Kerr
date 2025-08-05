clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function of baudrate @20 [ps/nm/km], 1500 [nm]

% Rs              = [10,15,20,25,30];
% DiffRs_10       = [225,375,525,625,850,950,1100,1250,1425,1525,1675,1850];
% DiffRs_15       = [165,225,295,370,440,500,535,630,700];
% DiffRs_20       = [54,95,125,175,200];
% DiffRs_25       = [35,59,83,100,125];
% DiffRs_30       = [25,40,55,74,91];
% 
% DiffRs_10_mean  = mean(DiffRs_10);
% DiffRs_15_mean  = mean(DiffRs_15);
% DiffRs_20_mean  = mean(DiffRs_20);
% DiffRs_25_mean  = mean(DiffRs_25);
% DiffRs_30_mean  = mean(DiffRs_30);
% 
% DiffRs          = [DiffRs_10_mean,DiffRs_15_mean,DiffRs_20_mean,DiffRs_25_mean,DiffRs_30_mean];
% 
% DiffRs_10_std   = std(DiffRs_10);
% DiffRs_15_std   = std(DiffRs_15);
% DiffRs_20_std   = std(DiffRs_20);
% DiffRs_25_std   = std(DiffRs_25);
% DiffRs_30_std   = std(DiffRs_30);
% StdRs            = [DiffRs_10_std,DiffRs_15_std,DiffRs_20_std,DiffRs_25_std,DiffRs_30_std];
% 
% x = Rs;
% y = DiffRs;
% % [fitResult,gof] = fit_pow(x, y);
% [fitResult,gof] = fit_pow_fixed_b(x, y,-2);
% xfit = 1:31;
% yfit = fitResult.a*xfit.^(fitResult.b);
% 
% 
% 
% subplot(1,3,1)
% hold on
% errorbar(Rs,DiffRs,StdRs)
% plot(xfit,yfit,'--')
% title(sprintf('y = ax^b, a = %.2e, b = %d, R^2 = %.2f',fitResult.a,fitResult.b,gof.rsquare))
% 
% xlim([0,35])
% ylim([0,2000])
% xlabel('Rs [GBd]')
% ylabel('half pseudo period')
% legend('data','fit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function of dispersion @25 [GBd], 1500 [nm]

D               = [2,5,10,15,20,25];%,30];
DiffD_2         = [295,510,695,905,1110,1295,1495,1710,1910,2115];
DiffD_5         = [117,202,282,362,437,522,607,677,765,842,922];
DiffD_10        = [137,180,220,257,300,340,382,412,455,505,545,585,620,655,700,740];
DiffD_15        = [120,150,170,205,230,255,295,310,330,365,390,410,435];
DiffD_20        = [30,50,70,90,112,128,154,170];
DiffD_25        = [24,40,56,72,88,104,120,136,154,170,184];
DiffD_30        = [20,34,46,62,72,88,100,112,128,140,152,168,184];

Diff_2          = diff(DiffD_2);
Diff_5          = diff(DiffD_5);
Diff_10         = diff(DiffD_10);
Diff_15         = diff(DiffD_15);
Diff_20         = diff(DiffD_20);
Diff_25         = diff(DiffD_25);
Diff_30         = diff(DiffD_30);

DiffD_2_mean    = mean(DiffD_2);
DiffD_5_mean    = mean(DiffD_5);
DiffD_10_mean   = mean(DiffD_10);
DiffD_15_mean   = mean(DiffD_15);
DiffD_20_mean   = mean(DiffD_20);
DiffD_25_mean   = mean(DiffD_25);
DiffD_30_mean   = mean(DiffD_30);
DiffD           = [DiffD_2_mean,DiffD_5_mean,DiffD_10_mean,DiffD_15_mean,...
                    DiffD_20_mean,DiffD_25_mean];%,DiffD_30_mean];

DiffD_2_std     = std(DiffD_2);
DiffD_5_std     = std(DiffD_5);
DiffD_10_std    = std(DiffD_10);
DiffD_15_std    = std(DiffD_15);
DiffD_20_std    = std(DiffD_20);
DiffD_25_std    = std(DiffD_25);
DiffD_30_std    = std(DiffD_30);
StdD            = [DiffD_2_std,DiffD_5_std,DiffD_10_std,DiffD_15_std,...
                    DiffD_20_std,DiffD_25_std];%,DiffD_30_std];

x = D;
y = DiffD;
% [fitResult,gof] = fit_pow(x, y);
[fitResult,gof] = fit_pow_fixed_b(x, y,-1);
% [a_opt,b_opt,R2] = fit_power_manual(x, y,-1);
xfit = 1:31;
yfit = fitResult.a*xfit.^(fitResult.b);
% yfit = a_opt*(xfit-b_opt).^(-1);

subplot(1,3,2)
hold on
errorbar(D,DiffD,StdD)
plot(xfit,yfit,'--')
title(sprintf('y = ax^b, a = %.2e, b = %d, R^2 = %.2f',fitResult.a,fitResult.b,gof.rsquare))
xlim([0,35])
ylim([0,2000])
xlabel('D [ps/nm/km]')
ylabel('half pseudo period')
legend('data','fit')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function of wavelength @20 [ps/nm/km], 25 [GBd]
% 
% Lam             = [1000,1100,1200,1300,1400,1500,1600]%,1700];
% DiffLam1000     = [290,510,700,910,1100,1290,1500,1690,1910,2150,2310,2550,2690];
% DiffLam1100     = [160,270,390,500,610,710,830,950];
% DiffLam1200     = [70,120,170,220,260,310,360,420,460,510,550,600];
% DiffLam1300     = [50,80,110,150,190,210,240,280,310,340,380,410];
% DiffLam1400     = [38,65,88,115,139,167,192];
% DiffLam1500     = [32,54,76,95,120,140,161,180];
% DiffLam1600     = [23,40,54,72,87,104,116,135,148,165,181,192];
% DiffLam1700     = [21,37,52,66,82,91,110,124,139,154,169,180,197];
% 
% Diff_1000       = diff(DiffLam1000);
% Diff_1100       = diff(DiffLam1100);
% Diff_1200       = diff(DiffLam1200);
% Diff_1300       = diff(DiffLam1300);
% Diff_1400       = diff(DiffLam1400);
% Diff_1500       = diff(DiffLam1500);
% Diff_1600       = diff(DiffLam1600);
% Diff_1700       = diff(DiffLam1700);
% 
% DiffLam1000_mean= mean(DiffLam1000);
% DiffLam1100_mean= mean(DiffLam1100);
% DiffLam1200_mean= mean(DiffLam1200);
% DiffLam1300_mean= mean(DiffLam1300);
% DiffLam1400_mean= mean(DiffLam1400);
% DiffLam1500_mean= mean(DiffLam1500);
% DiffLam1600_mean= mean(DiffLam1600);
% DiffLam1700_mean= mean(DiffLam1700);
% % DiffLam1800_mean= mean(DiffLam1800);
% % DiffLam1900_mean= mean(DiffLam1900);
% DiffLam         = [DiffLam1000_mean,DiffLam1100_mean,DiffLam1200_mean,DiffLam1300_mean,...
%                     DiffLam1400_mean,DiffLam1500_mean,DiffLam1600_mean];%,DiffLam1700_mean];
% 
% DiffLam1000_std = std(DiffLam1000);
% DiffLam1100_std = std(DiffLam1100);
% DiffLam1200_std = std(DiffLam1200);
% DiffLam1300_std = std(DiffLam1300);
% DiffLam1400_std = std(DiffLam1400);
% DiffLam1500_std = std(DiffLam1500);
% DiffLam1600_std = std(DiffLam1600);
% DiffLam1700_std = std(DiffLam1700);
% 
% StdLam            = [DiffLam1000_std,DiffLam1100_std,DiffLam1200_std,DiffLam1300_std,...
%                     DiffLam1400_std,DiffLam1500_std,DiffLam1600_std,DiffLam1700_std];
% 
% 
% X = Lam;
% Y = DiffLam;
% 
% % z = StdLam;
% 
% % M = [x',y',z'];
% % writematrix(M,'M_diff_lambda.txt')
% % [fitResult,gof] = fit_pow(x, y);
% % [fitResult,gof] = fit_pow_fixed_c(x, y,-2);
% xfit = linspace(min(x),max(x),100);
% Xfit = linspace(min(X),max(X),100);
% % yfit = 10^(7.7)*(xfit-800).^(-2);
% 
% 
% [a_opt, b_opt, R2] = fit_power2_manual(X, Y);
% Yfit = a_opt*(Xfit-b_opt).^(-2);
% 
% 
% [a_opt_n, b_opt_n, R2_n] = fit_power2_manual(x, y);
% yfit = a_opt_n*(xfit-b_opt_n).^(-2);
% 
% subplot(1,2,1)
% hold on
% plot(Lam,DiffLam)
% plot(Xfit,Yfit,'--')
% % title(sprintf('y = ax^b, a = %.2e, b = %d, R^2 = %.2f',fitResult.a,fitResult.b,gof.rsquare))
% % xlim([0,35])
% % ylim([0,2000])
% xlabel('D [ps/nm/km]')
% ylabel('half pseudo period')
% % legend('data','fit')
% set(gca,'yscale','log')
% 
% subplot(1,2,2)
% hold on
% plot(x,y)
% plot(xfit,yfit,'--')
% % title(sprintf('y = ax^b, a = %.2e, b = %d, R^2 = %.2f',fitResult.a,fitResult.b,gof.rsquare))
% % xlim([0,35])
% % ylim([0,2000])
% xlabel('D [ps/nm/km]')
% ylabel('half pseudo period')
% % legend('data','fit')
% % set(gca,'yscale','log')
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a_opt, b_opt, R2] = fit_power_manual(x, y,c)
    % Fonction de coût : erreur quadratique
    cost_fun = @(p) sum((y - p(1) * (x - p(2)).^(-2)).^2);

    % Estimation initiale : [a, b]
    a0 = max(y) * min(x)^2;
    b0 = min(x) - 100;

    % Optimisation : recherche de [a, b]
    opts = optimset('Display','off');
    params_opt = fminsearch(cost_fun, [a0, b0], opts);

    % Paramètres optimisés
    a_opt = params_opt(1);
    b_opt = params_opt(2);

    % Valeurs ajustées
    yfit = a_opt * (x - b_opt).^(c);

    SS_res = sum((y - yfit).^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = 1 - SS_res / SS_tot;
end



function [fitResult, gof] = fit_pow(x, y)
    % Définir le modèle : tous les paramètres sont libres
    fitType = fittype('a * (x - b).^c', ...
                      'independent', 'x', ...
                      'coefficients', {'a', 'b', 'c'});

    % Définir les options de fit
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares');
    fitOptions.StartPoint = [1, 0, 1];  % Valeurs initiales : [a, b, c]

    % Ajustement
    [fitResult, gof] = fit(x(:), y(:), fitType, fitOptions);
end

function [fitResult, gof] = fit_pow_fixed_b(x, y, b_fixed)
    % Définir le modèle avec b comme paramètre "problème" (fixe)
    fitType = fittype('a * x.^b', ...
                      'independent', 'x', ...
                      'coefficients', 'a', ...
                      'problem', 'b');

    % Définir les options de fit
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares');
    % Initialiser a (par exemple à 1)
    fitOptions.StartPoint = 10;

    % Lancer l'ajustement en passant b_fixed comme problème
    [fitResult, gof] = fit(x(:), y(:), fitType, fitOptions, 'problem', b_fixed);
end

function [fitResult, gof] = fit_pow_fixed_c(x, y, c_fixed)
    % Définir le modèle avec c comme paramètre "problème" (fixe)
    fitType = fittype('a * (x - b).^c', ...
                      'independent', 'x', ...
                      'coefficients', {'a', 'b'}, ...
                      'problem', 'c');

    % Définir les options de fit
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares');
    % Valeurs initiales pour a et b (à ajuster au besoin)
    fitOptions.StartPoint = [0,0];

    % Ajustement avec c fixé
    [fitResult, gof] = fit(x(:), y(:), fitType, fitOptions, 'problem', c_fixed);
end
