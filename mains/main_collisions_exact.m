clear
close all
clc

lambda  = 1550e-9;
c       = 299792458;
dnu     = 12.5e9;
nu      = c/lambda;
dlambda = dnu/nu*lambda;

D       = 17e-6;
% Rtest   = 90;
% Rs      = (10:1:2000)*1e9;
Rs      = logspace(0,3,10)*1e9;
Ls      = [(1:1:9),(10:1:100),(125:25:1000),(1000:100:8000)]*1e3;
NLs     = length(Ls);

NRs     = length(Rs);

M       = zeros(NLs,NRs);
dt_isi_l= zeros(1,NLs);
Sisi    = zeros(NLs,NRs);
Snisi   = zeros(NLs,NRs);
for i = 1:NLs
    L           = Ls(i);

    dt_isi      = D * L * dlambda;

    for j = 1:NRs

        N_isi   = ceil(dt_isi * Rs(j));
        r_isi_l = zeros(1,N_isi);

        % utilisé pour le manuscrit
        % for m = 1:N_isi-1
        %     r_isi_l(m) = (1 - m/N_isi *exp(-m^2/4);%* cos(2*pi*nu/Rs(j)*m);
        %     Sisi(i,j) = Sisi(i,j) + (N_isi+1-m)* exp(-m^2/4);
        % end


        for m = 1:N_isi
            r_isi_l(m) = 2*(1 - m/(N_isi+1)) *exp(-m^2/4);%* cos(2*pi*nu/Rs(j)*m);
            Sisi(i,j) = Sisi(i,j) + (N_isi+1-m)* exp(-m^2/4);
        end

        r_isi   = sum(r_isi_l);
        M(i,j)  = r_isi-r_isi_l(1);
        % M(i,j)
        Snisi(i,j) = N_isi+1;
    end
end

% Diff_S1 = abs(Sisi-Snisi)./Snisi*100;
% Diff_S2 = abs(Sisi-Snisi)./Sisi*100;
% 
% subplot(2,2,1)
%     plot(Rs*1e-9,mean(Diff_S1,1))
%     xlabel('Rs [Gbd]')
%     set(gca,'Xscale','log')
% 
% subplot(2,2,2)
%     plot(Ls*1e-3,mean(Diff_S1,2))
%     set(gca,'Xscale','log')
%     xlabel('L [km]')
% 
% subplot(2,2,3)
%     plot(Rs*1e-9,mean(Diff_S2,1))
%     set(gca,'Xscale','log')
%     xlabel('Rs [Gbd]')
% 
% subplot(2,2,4)
%     plot(Ls*1e-3,mean(Diff_S2,2))
%     set(gca,'Xscale','log')
%     xlabel('L [km]')

line_styles = {'-', '--', ':', '-.'};
markers     = {'o', 's', 'd', '^', 'v', 'x', '+', '*', '>'};
Nstyles     = length(line_styles);
Nmarkers    = length(markers);


figure
hold on

for j = 1:NRs
    % Index pour style et marqueur
    ls = line_styles{mod(j-1, Nstyles)+1};
    mk = markers{mod(j-1, Nmarkers)+1};

    % Épaisseur inversement proportionnelle au débit
    lw = 2.5 - 1.5 * (j-1)/(NRs-1);  % entre 2.5 et 1

    plot(Ls*1e-3, M(:,j), ...
        'LineStyle', ls, ...
        'Marker', mk, ...
        'LineWidth', lw, ...
        'MarkerIndices', 1:10:NLs, ...
        'DisplayName', sprintf('%3d GBd', round(Rs(j)*1e-9)))
end

hold off
xlim([1,8e3])
xlabel('Fibre length (km)')
ylabel('ISI metric')
% legend('show','Location','southeast','numcolumns',2)
title('ISI metric vs Fibre Length for different Symbol Rates')
set(gca,'fontsize',20,'Xscale','log')
grid on


% T = array2table(M, 'VariableNames', strcat('G',string(round(Rs*1e-9))));
% T = addvars(T, (Ls')*1e-3, 'Before', 1, 'NewVariableNames', 'Lkm');
% 
% fname = sprintf('isi_table_by_baudrate_th.txt');
% writetable(T, fname, 'Delimiter', ',')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to plot the theoretical limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nmax    = 50;
% res     = zeros(1,Nmax-2);
% 
% for n = 3:Nmax
%     N       = n;
%     ArgErf  = (N-1)/2;
% 
%     res(n-2)  = sqrt(pi)*erf(ArgErf)+2/N*exp(-((N-1)/2)^2)-2/N-1;
% end
% 
% u_1 = (1-1/N)*exp(-0.25);
% 
% figure
% plot(3:Nmax,res-u_1)
% 
% M = [(3:Nmax)',res'];
% T = array2table(M);
% 
% T.Properties.VariableNames = {'nPulses', 'Bound'};
% writetable(T, 'upper_bound.txt', 'Delimiter', ',')

% N = 10;
% m = 2:N;
% a = 1;
% vm = exp(-(m/(2*a)).^2);
% v1 = exp(-(1/(2*a)).^2);
% wm = v1-vm;
% rISI_bound = sqrt(pi)*(erfc(1/4)-erfc(N-1/4))-(1-1./N).*exp(-1/4);
% figure
% plot(N,rISI_bound)
% set(gca,'Xscale','log')
% 
% M = [N',rISI_bound'];
% T = array2table(M);
% 
% T.Properties.VariableNames = {'nPulses', 'Bound'};
% writetable(T, 'upper_bound_rI
% wm2 = exp(-(m.^2+1)/(8*a^2)).*sinh((1-m.^2)./(8*a^2));
% 
% um = (1-m/N).*vm;
% u1 = (1-1/N).*v1;
% wmm = u1-um;
% wmm2 = (m-1)/N.*exp(-(m/(2*a)).^2);
% 
% figure
% hold on
% plot(wmm)
% plot(wmm2)

% N = 2:1000;
% rISI_bound = sqrt(pi)*(erfc(1/4)-erfc(N-1/4))-(1-1./N).*exp(-1/4);
% figure
% plot(N,rISI_bound)
% set(gca,'Xscale','log')
% 
% M = [N',rISI_bound'];
% T = array2table(M);
% 
% T.Properties.VariableNames = {'nPulses', 'Bound'};
% writetable(T, 'upper_bound_rISI.txt', 'Delimiter', ',')