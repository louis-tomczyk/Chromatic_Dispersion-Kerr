clear
close all
clc
addpath '/home/louis/Documents/6_TélécomParis/3_Codes/0_louis/1_ClassicDSP/ISI_Kerr/libraries'/optilux/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global parameters
Nsymb   = 2^8;               % number of symbols
Nt      = 2;                 % number of discrete points per symbol
Nrea    = 2;                % number of realizations

Ls      = (1:0.5:1e3)*1e3;    % distance range [m]
Nlength = length(Ls);

Rs      = 25;                 % symbol rate [Gbaud]
symbrate= Rs;                 % alias for compatibility

modfors = {'qpsk'};            % modulation formats
disps   = [20];                % dispersion values [ps/nm/km]
Ndisps  = length(disps);
lams    = [1550];%,1400,1500,1600,1700,1800,1900];
Nlam    = length(lams);

pds     = 0;%[340,680,1020,1360,1700,3400,5100];
Npds    = length(pds);

n0s     = [-100,-90,-80,-70,-60,-50];
Nn0s    = length(n0s);

Rfs     = [0,0.05];%,0.1,0.2,0.3,0.4,0.5];
NRfs    = length(Rfs);

dnus     = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6];
Ndnus    = length(dnus);

% Preallocate results matrix
M       = zeros(Nrea,Nlength,NRfs);
Msave   = zeros(NRfs,Nlength);
% t = zeros(1,Nlength);

%% Main simulation loop
for i = 1:NRfs
    sprintf('roll off = %.2f',round(Rfs(i),2))

    Nsamp           = Nsymb*Nt;        % total number of samples
    fs              = symbrate*Nt;     % sampling rate [GHz]

    modfor          = modfors{1};      % modulation format

    PdBm            = 0;               % power [dBm]
    lam             = lams(1);            % carrier wavelength [nm]
    Plin            = 10^(PdBm/10);    % power [mW]

    tx.rolloff      = Rfs(i);               % pulse roll-off
    tx.emph         = 'asin';          % digital-premphasis type

    ft.lambda       = 1550;            % wavelength [nm]
    ft.alphadB      = 0;               % attenuation [dB/km]
    ft.disp         = disps(1);        % dispersion [ps/nm/km]    
    ft.slope        = 0;               % slope [ps/nm^2/km]
    ft.aeff         = 50;              % effective area [um^2]
    ft.n2           = 0;               % nonlinear index [m^2/W]
    
    fpd             = ft;
    fpd.length      = 1e3;
    fpd.disp        = pds(1);

    rx.modformat    = modfor;          % modulation format
    rx.sync.type    = 'da';            % time-recovery method
    rx.oftype       = 'gauss';         % optical filter type
    rx.obw          = Inf;             % optical filter bandwidth
    rx.eftype       = 'rootrc';      % electrical filter type % costails
    rx.ebw          = 0.5;             % electrical filter bandwidth
    rx.epar         = tx.rolloff;
    rx.type         = 'bin';           % binary pattern

    for k = 1:Nlength
        tempM = zeros(Nrea,1);
        ft.length      = Ls(k);           % length [m]
        % tic
        for j = 1:Nrea
            inigstate(Nsamp,fs);            % initialize global variables
            E            = lasersource(Plin,lam,struct('pol','single'));
            % E            = lasersource(Plin,lam,struct('pol','single'));%,'n0',n0s(i)));
            % E            = lasersource(Plin,lam,struct('pol','single','linewidth',dnus(i)));


            [patx, patbinx] = datapattern(Nsymb,'rand',struct('format',modfor));
            [sigx, normx]   = digitalmod(patbinx,modfor,symbrate,rx.eftype,tx);

            %% Channel
            Ein             = iqmodulator(E, sigx,struct('norm',normx));
            Epd             = fiber(Ein,fpd);
            Eout            = fiber(Epd,ft);

            dapr_in         = std(empower(Ein))/mean(empower(Ein));
            dapr_out        = std(empower(Eout))/mean(empower(Eout));
            tempM(j)        = dapr_out-dapr_in;
        end
        % t(k) = toc;

        M(:,k,i) = tempM;
    end

    % Plot results for current dispersion
    hold on
    % plot(Ls*1e-3,mean(M(:,:,i)),'DisplayName',sprintf('dnu=%.1e GHz',dnus(i)));
    % plot(Ls*1e-3,mean(M(:,:,i)),'DisplayName',sprintf('n0=%.1e dB/GHz',n0s(i)));
    % plot(Ls*1e-3,mean(M(:,:,i)),'DisplayName',sprintf('pd=%.1e ps/nm',pds(i)));
    plot(Ls*1e-3,mean(M(:,:,i)),'DisplayName',sprintf('rolloff=%.2f',Rfs(i)));
    Msave(i,:) = mean(M(:,:,i));
end


% for i = 1:Nn0s
%     hold on
%     plot(Ls*1e-3,mean(M(:,:,i)),'DisplayName',sprintf('n0=%d dB/GHz',n0s(i)));
% end

%% Plot formatting
set(gca,'Xscale','log')
xlabel('Distance de propagation [km]')
ylabel('dapr(out) - dapr(in)')
title(sprintf('%s, att = %.1f dB/km, gamma = 0',modfor,ft.alphadB))
set(gca,'fontsize',20)
legend show
grid on

% T = array2table(Msave', 'VariableNames',strcat('G',string(round(Rs))));
% T = array2table(Msave', 'VariableNames',modfors);
% T = array2table(Msave', 'VariableNames',strcat('L',string(round(lams))))
% T = array2table(Msave', 'VariableNames',strcat(["Dnu1G","Dnu100M","Dnu10M","Dnu1M","Dnu100k","Dnu10k","Dnu1k"]);
% T = array2table(Msave', 'VariableNames',strcat('d',string(abs(n0s))));
% T = array2table(Msave', 'VariableNames',strcat('pd',string(pds)));
% T = addvars(T, (Ls')*1e-3, 'Before', 1, 'NewVariableNames', 'Lkm');

T = array2table(Msave', 'VariableNames',strcat('rolloff',string(Rfs)));
T = addvars(T, (Ls')*1e-3, 'Before', 1, 'NewVariableNames', 'Lkm');
% 
% % writetable(T, 'isi_table_by_baudrate_simu_modfor.txt', 'Delimiter', ',')
% writetable(T, 'isi_table_by_baudrate_simu_ook_lambda.txt', 'Delimiter', ',');
% writetable(T, 'isi_table_by_baudrate_simu_qpsk_dnu.txt', 'Delimiter', ',');
% writetable(T, 'isi_table_qpsk_by_baudrate_simu_pd_01_1.txt', 'Delimiter', ',');

writetable(T, 'isi_table_qpsk_by_baudrate_simu_rolloff_0_1_test.txt', 'Delimiter', ',');




