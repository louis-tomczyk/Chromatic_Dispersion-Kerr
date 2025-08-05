clear
close all
clc
addpath '/home/louis/Documents/6_TélécomParis/3_Codes/0_louis/1_ClassicDSP/ISI_Kerr/libraries'/optilux/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global parameters

Nsymb   = 2^8;                 % number of symbols
Nt      = 2;                   % number of discrete points per symbol


Nrea    = 2;

Ls      = (1e2:5:1e4)*1e3;
Nlength = length(Ls);

Rs      = 25;
NRs     = length(Rs);
M       = zeros(Nrea,Nlength,NRs);
Msave   = zeros(NRs,Nlength);
modfors = {'ook'};
Nmodfor = 3;

disps   = [17];
Ndisps  = length(disps);


% for i = 1:Nmodfor
% for i = 1:NRs
for i = 1:Ndisps

    % sprintf('Rs=%i',Rs(i))
    sprintf('D = %i',disps(i))
    
    for j = 1:Nrea
        for k = 1:Nlength
        
        
            %% Tx parameters
            symbrate        = Rs;        % symbol rate [Gbaud].

            tx.rolloff      = 0;            % pulse roll-off
            tx.emph         = 'asin';       % digital-premphasis type
            % modfor          = modfors{i};       % modulation format
            modfor          = modfors{1};       % modulation format
            PdBm            = 0;            % power [dBm]
            lam             = 1550;         % carrier wavelength [nm]
            
            %% Channel parameters            
            ft.length       = Ls(k);        % length [m]
            ft.lambda       = 1550;         % wavelength [nm] of fiber parameters
            ft.alphadB      = 0;          % attenuation [dB/km]
            ft.disp         = disps(i);           % dispersion [ps/nm/km] @ ft.lambda
            ft.slope        = 0;            % slope [ps/nm^2/km] @ ft.lambda
            ft.aeff         = 50;           % effective area [um^2]
            
            gamma2          = 0.002;%2*pi/ft.lambda*ft.n2/ft.aeff*1e12*1e9*1000; % [1/W/km]
            n2              = 0;%(ft.lambda*1e-9)/2/pi*(ft.aeff*1e-12)*gamma2;
            ft.n2           = n2;%2.6*1e-20;    % nonlinear index [m^2/W]

            fpd = ft;
            fpd.length      = 1e3;
            fpd.disp        = 1700;

            %% Rx parameters
            rx.modformat    = modfor;           % modulation format
            rx.sync.type    = 'da';             % time-recovery method
            rx.oftype       = 'gauss';          % optical filter type
            rx.obw          = Inf;              % optical filter bandwidth normalized to symbrate
            rx.eftype       = 'costails';         % electrical filter type
            rx.ebw          = 0.5;              % electrical filter bandwidth normalized to symbrate
            rx.epar         = tx.rolloff;
            rx.type         = 'bin';            % binary pattern
            
            %% Init
            Nsamp           = Nsymb*Nt;         % overall number of samples
            fs              = symbrate*Nt;      % sampling rate [GHz]
            inigstate(Nsamp,fs);                % initialize global variables: Nsamp and fs.
            
            %% Tx side
            
            Plin            = 10.^(PdBm/10);   % [mW]
            E               = lasersource(Plin,lam,struct('pol','single'));%,'linewidth',1e-4));  % y-pol does not exist
            
            [patx, patbinx] = datapattern(Nsymb,'rand',struct('format',modfor));
            [sigx, normx]   = digitalmod(patbinx,modfor,symbrate,'rootrc',tx);
            
            %% Channel
            
            Ein             = iqmodulator(E, sigx,struct('norm',normx));
            Epd            = fiber(Ein,fpd);
            Eout            = fiber(Ein,ft);
            
            dapr_in         = std(empower(Ein))/mean(empower(Ein));
            dapr_out        = std(empower(Eout))/mean(empower(Eout));
            M(j,k,i)        = dapr_out-dapr_in;
        
        end
    end
    hold on
    plot(Ls*1e-3,mean(M(:,:,i)),DisplayName=sprintf('%d [GBd]',symbrate));
    % plot(Ls*1e-3,mean(M(:,:,i)),DisplayName=sprintf('%s',modfors{i}));
    Msave(i,:) = mean(M(:,:,i));
end

set(gca,'Xscale','log')
xlabel('distance of propagation [km]')
ylabel('dapr(out) - dapr(in)')
t=title(sprintf('%s, att = %.1f [dB/km], disp = %d [ps/nm/km], gamma = %.1f',modfor,ft.alphadB,ft.disp,gamma2*1e3));
set(gca,'fontsize',25)
% ylim([0,0.5])
set(gca,'fontsize',20)
% T = array2table(Msave', 'VariableNames',strcat('G',string(round(Rs))));
% % T = array2table(Msave', 'VariableNames',modfors);
% T = addvars(T, (Ls')*1e-3, 'Before', 1, 'NewVariableNames', 'Lkm');
% 
% % writetable(T, 'isi_table_by_baudrate_simu_modfor.txt', 'Delimiter', ',')
% writetable(T, 'isi_table_by_baudrate_simu_ook.txt', 'Delimiter', ',');
