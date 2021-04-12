
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   Code to check the impact of estmation error of Bob's channel performed
%   at Alice. Due to the error, the precoding does not lead to a perfect
%   matched filter at Bob, reducing the SR performances. 
%   Plot of the SR as a function of alpha for different BOR and different
%   decoding structures at Eve.
%
%   Code Started:   19.11.2020
%   Last Update:    15.03.2021
%
%
%  © SIDNEY GOLSTEIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%close all;

set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',32)
set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'DefaultLineMarkerSize',15);
set(0, 'defaultFigurePosition',  [-1267  44   1256    872])

h = waitbar(0,'Simulation Progression...');

%% Parameters
% Simulation parameters
nb_run = 1000;              % number of experiments
alpha_step = 4;             % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;  

sigma_tilde = [0:0.02:0.5];                % sigma_tilde := percentage of CSI error --> = 0 : perfect CSI @Alice , = 1: 100% error
sigma_tilde_db = 10*log10(sigma_tilde);     % CSI error in dB.

% Communication parameters
Q = 16;
U = [2 4 8 16];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
snr_b  = [10]; %EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = [20]; %EbN0_e + 10*log10(k);  % SNR @Eve
sigma_b = 1./U.'./10.^(snr_b./10);    % expected noise energy @Bob
sigma_e = 1./U.'./10.^(snr_e./10);    % expected noise energy @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance



%% Matrices preallocation
e_an_TX             = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));

e_sym_decod1_b      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_noise_decod1_b    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_an_decod1_b       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_denom_decod1_b    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));

e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
% 
% e_sym_decod3_e      = zeros(nb_run,length(alpha),length(U));
% e_noise_decod3_e    = zeros(nb_run,length(alpha),length(U));
% e_an_decod3_e       = zeros(nb_run,length(alpha),length(U));
% e_denom_decod3_e    = zeros(nb_run,length(alpha),length(U));
% 
% e_sym_decod4_e      = zeros(nb_run,length(alpha),length(U));
% e_noise_decod4_e    = zeros(nb_run,length(alpha),length(U));
% e_an_decod4_e       = zeros(nb_run,length(alpha),length(U));
% e_denom_decod4_e    = zeros(nb_run,length(alpha),length(U));
% 
e_sym_decod5_e      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_noise_decod5_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_an_decod5_e       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));
e_denom_decod5_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde),length(snr_b),length(snr_e));


%% Mainloop
for iter = 1:nb_run
for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run




% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));

% Channel generation
% Rayleigh channel - independant subcarrier fading - unit power per
% subcarrier : CN(0,I_Q) - circulary complex random normal variable.

% Bob channel
Hb_TX = channelRayleigh(Q, mu , sigma);                                     % Correct Bob's channel 
Hb_delta_TX = channelRayleigh(Q, mu , sigma);                               % Estimation error
Hb_delta_RX = ctranspose(Hb_delta_TX.');
Hb_RX = ctranspose(Hb_TX.');

% Eve channel
He_TX = channelRayleigh(Q, mu , sigma);                                     % Eve channel 
He_RX = ctranspose(He_TX.');

for ss = 1:length(sigma_tilde)
% Bob channel estimated by Alice (imperfect CSI)
Hb_tilde_TX = sqrt(1-sigma_tilde(ss))*Hb_TX + sqrt(sigma_tilde(ss))*Hb_delta_TX;    % Estimate Bob CSI with imperfect CSI.
Hb_tilde_RX = ctranspose(Hb_tilde_TX.');

%% Encoder
% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % Qx1
sym_precoded = Hb_tilde_TX*sym_spread; % Qx1, not weighted 

% AN generation
an = generateAN(Hb_tilde_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd"); % Qx1, not weighted

for aa = 1:length(alpha)
% Weighting
sym_precoded_TX = sqrt(alpha(aa))*sym_precoded;     % weighted
an_TX = sqrt(1-alpha(aa))*an;                       % weighted
e_an_TX(iter,aa) = energy(an_TX);

%% Receiver
% Whole stream
sym_RX_b = Hb_RX*(sym_precoded_TX + an_TX);
sym_RX_e = He_RX*(sym_precoded_TX + an_TX);

% Useful symbol
sym_b = Hb_RX*sym_precoded_TX; % Qx1
sym_e = He_RX*sym_precoded_TX;



for nb = 1:length(snr_b)
for ne = 1:length(snr_e)
% Noise symbol
[noise_b, ~ ] = addNoise(sym_b , snr_b(nb), energy(sym_precoded_TX+an_TX));     %addNoise(sym_b , snr_b, energy(sym_RX_b)); 
[noise_e, e_noise_e ] = addNoise(sym_e , snr_e(ne), energy(sym_precoded_TX+an_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   

% AN symbol
an_b = Hb_RX*an_TX; % @Bob : not 0 since imperfect CSI
an_e = He_RX*an_TX; % @Eve



%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_tilde_RX*He_TX;                 % matched filter : Eve perfectly estimate the channel Hb_tilde_RX He_TX since bad precoding @A
% decod3 = matrix_despread*(Hb_tilde/He_RX);                % AN killer
% gamma_E = (He_RX*Hb_TX)*matrix_spread;
% gamma_EH = ctranspose(gamma_E);
% decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an) + e_noise_e*eye(Q) );   % LMMSE
decod5 = matrix_despread*He_TX;                                   % Only He known by Eve


sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_decod2_e = decod2*sym_e;
% sym_decod3_e = decod3*sym_e;
% sym_decod4_e = decod4*sym_e;
sym_decod5_e = decod5*sym_e;

noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;
noise_decod2_e = decod2*noise_e;
% noise_decod3_e = decod3*noise_e;
% noise_decod4_e = decod4*noise_e;
noise_decod5_e = decod5*noise_e;

an_decod1_b = decod1*an_b;
an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;
% an_decod3_e = decod3*an_e;
% an_decod4_e = decod4*an_e;
an_decod5_e = decod5*an_e;
% 

% sinr_no_decod_e(iter,aa,:) = U*diag((ctranspose(He_equiv)*He_equiv)/(ctranspose(an_decod2_e+noise_decod2_e)*(an_decod2_e+noise_decod2_e)));


%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,aa,bb,ss,nb,ne)      = energy(sym_decod1_b);
e_noise_decod1_b(iter,aa,bb,ss,nb,ne)    = energy(noise_decod1_b);
e_an_decod1_b(iter,aa,bb,ss,nb,ne)       = energy(an_decod1_b);
e_denom_decod1_b(iter,aa,bb,ss,nb,ne)    = energy(noise_decod1_b + an_decod1_b);        % energy of the sinr denominator for decoder 1 @Eve


% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,ss,nb,ne)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,ss,nb,ne)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,ss,nb,ne)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,ss,nb,ne)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,aa,bb,ss,nb,ne)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb,ss,nb,ne)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb,ss,nb,ne)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb,ss,nb,ne)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve
% 
% % @ Eve : decod 3
% e_sym_decod3_e(iter,aa,bb)     = energy(sym_decod3_e);
% e_noise_decod3_e(iter,aa,bb)   = energy(noise_decod3_e);
% e_an_decod3_e(iter,aa,bb)      = energy(an_decod3_e);
% e_denom_decod3_e(iter,aa,bb)   = energy(noise_decod3_e + an_decod3_e);        % energy of the sinr denominator for decoder 3 @Eve
% 
% % @ Eve : decod 4
% e_sym_decod4_e(iter,aa,bb)     = energy(sym_decod4_e);
% e_noise_decod4_e(iter,aa,bb)   = energy(noise_decod4_e);
% e_an_decod4_e(iter,aa,bb)      = energy(an_decod4_e);
% e_denom_decod4_e(iter,aa,bb)   = energy(noise_decod4_e + an_decod4_e);        % energy of the sinr denominator for decoder 4 @Eve
% 
% @ Eve : decod 5
e_sym_decod5_e(iter,aa,bb,ss,nb,ne)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,aa,bb,ss,nb,ne)   = energy(noise_decod5_e);
e_an_decod5_e(iter,aa,bb,ss,nb,ne)      = energy(an_decod5_e);
e_denom_decod5_e(iter,aa,bb,ss,nb,ne)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve

end
end
end
waitbar(iter / nb_run)
end
end
end
close(h)

%% Post processing
% Energy averaging over all realisations
e_avg_an                = squeeze(mean(e_an_TX,1));

e_avg_sym_decod1_b      = squeeze(mean(e_sym_decod1_b,1));  
e_avg_noise_decod1_b    = squeeze(mean(e_noise_decod1_b,1));
e_avg_an_decod1_b       = squeeze(mean(e_an_decod1_b,1));
e_avg_denom_decod1_b    = squeeze(mean(e_denom_decod1_b,1));


e_avg_sym_decod1_e      = squeeze(mean(e_sym_decod1_e,1));
e_avg_noise_decod1_e    = squeeze(mean(e_noise_decod1_e,1));
e_avg_an_decod1_e       = squeeze(mean(e_an_decod1_e,1));
e_avg_denom_decod1_e    = squeeze(mean(e_denom_decod1_e,1));

e_avg_sym_decod2_e      = squeeze(mean(e_sym_decod2_e,1));
e_avg_noise_decod2_e    = squeeze(mean(e_noise_decod2_e,1));
e_avg_an_decod2_e       = squeeze(mean(e_an_decod2_e,1));
e_avg_denom_decod2_e    = squeeze(mean(e_denom_decod2_e,1));
% 
% e_avg_sym_decod3_e      = squeeze(mean(e_sym_decod3_e,1));
% e_avg_noise_decod3_e    = squeeze(mean(e_noise_decod3_e,1));
% e_avg_an_decod3_e       = squeeze(mean(e_an_decod3_e,1));
% e_avg_denom_decod3_e    = squeeze(mean(e_denom_decod3_e,1));
% 
% e_avg_sym_decod4_e      = squeeze(mean(e_sym_decod4_e,1));
% e_avg_noise_decod4_e    = squeeze(mean(e_noise_decod4_e,1));
% e_avg_an_decod4_e       = squeeze(mean(e_an_decod4_e,1));
% e_avg_denom_decod4_e    = squeeze(mean(e_denom_decod4_e,1));
% 
e_avg_sym_decod5_e      = squeeze(mean(e_sym_decod5_e,1));
e_avg_noise_decod5_e    = squeeze(mean(e_noise_decod5_e,1));
e_avg_an_decod5_e       = squeeze(mean(e_an_decod5_e,1));
e_avg_denom_decod5_e    = squeeze(mean(e_denom_decod5_e,1));

% instantaneous SINRs
sinr1_b = e_sym_decod1_b./e_denom_decod1_b;
sinr1_e = e_sym_decod1_e./e_denom_decod1_e;
sinr2_e = e_sym_decod2_e./e_denom_decod2_e;
% sinr3_e = e_sym_decod3_e./e_denom_decod3_e;
% sinr4_e = e_sym_decod4_e./e_denom_decod4_e;
sinr5_e = e_sym_decod5_e./e_denom_decod5_e;


% ergodic SINR
sinr1_b_avg = squeeze(mean(sinr1_b,1));
sinr1_e_avg = squeeze(mean(sinr1_e,1));
sinr2_e_avg = squeeze(mean(sinr2_e,1));
sinr5_e_avg = squeeze(mean(sinr5_e,1));


% instantaneous Secrecy capacity
sr1 = secrecyCapacity(sinr1_b,sinr1_e);
sr2 = secrecyCapacity(sinr1_b,sinr2_e);
% sr3 = secrecyCapacity(sinr1_b,sinr3_e);
% sr4 = secrecyCapacity(sinr1_b,sinr4_e);
sr5 = secrecyCapacity(sinr1_b,sinr5_e);
% 

% ergodic Secrecy capacity
sr1_avg = squeeze(mean(sr1,1));
sr2_avg = squeeze(mean(sr2,1));
% sr3_avg = squeeze(mean(sr3,1));
% sr4_avg = squeeze(mean(sr4,1));
sr5_avg = squeeze(mean(sr5,1));





%% SINR modeling: % PCSI : perfect CSI / ICSI : imperfect CSI

% Matrix allocation

sinr1_ICSI_model_b = zeros(length(alpha),length(U),length(sigma_tilde),length(snr_b));
sinr1_ICSI_model_e = zeros(length(alpha),length(U),length(sigma_tilde),length(snr_e));
sinr2_ICSI_model_e = zeros(length(alpha),length(U),length(sigma_tilde),length(snr_e));
sinr5_ICSI_model_e = zeros(length(alpha),length(U),length(sigma_tilde),length(snr_e));


for bb = 1:length(U)
for ss = 1:length(sigma_tilde)           
for nb = 1:length(snr_b)
sinr1_ICSI_model_b(:,bb,ss,nb) = sinrModelingICSI(alpha,U(bb),snr_b(nb),NaN,sigma_tilde(ss),"bob_decod1");
end
for ne = 1:length(snr_e)
%sinr1_PCSI_model_b(:,bb,nb,ne) = sinrModeling(alpha,U(bb),snr_b(nb),snr_e(ne),1,1,1,"bob_decod1");              
sinr1_ICSI_model_e(:,bb,ss,ne) = sinrModelingICSI(alpha,U(bb),NaN,snr_e(ne),sigma_tilde(ss),"eve_decod1");
sinr2_ICSI_model_e(:,bb,ss,ne) = sinrModelingICSI(alpha,U(bb),NaN,snr_e(ne),sigma_tilde(ss),"eve_decod2");
sinr5_ICSI_model_e(:,bb,ss,ne) = sinrModelingICSI(alpha,U(bb),NaN,snr_e(ne),sigma_tilde(ss),"eve_decod5");
end
end
end
% 
%sr1_PCSI_model = secrecyCapacity(sinr1_PCSI_model_b,sinr1_ICSI_model_e);
for bb = 1:length(U)
for ss = 1:length(sigma_tilde)           
for nb = 1:length(snr_b)
for ne = 1:length(snr_e)
sr1_ICSI_model(:,bb,ss,nb,ne) = secrecyCapacity(sinr1_ICSI_model_b(:,bb,ss,nb),sinr1_ICSI_model_e(:,bb,ss,ne));
sr2_ICSI_model(:,bb,ss,nb,ne) = secrecyCapacity(sinr1_ICSI_model_b(:,bb,ss,nb),sinr2_ICSI_model_e(:,bb,ss,ne));
sr5_ICSI_model(:,bb,ss,nb,ne) = secrecyCapacity(sinr1_ICSI_model_b(:,bb,ss,nb),sinr5_ICSI_model_e(:,bb,ss,ne));
end
end
end
end


%% PLOT Section


% 1. SR as a function of AN for a given ICSI: model vs simu
nb = find(snr_b == 10);
ne = find(snr_e == 20);
bb = find(U == 4);
ss = find(sigma_tilde == 0.1);
set(0,'defaultAxesFontSize',26)


figure;
plot(100*(1-alpha),sr1_avg(:,bb,ss,nb,ne),'Marker','o'); hold on;
plot(100*(1-alpha),sr1_ICSI_model(:,bb,ss,nb,ne),'Marker','<'); hold on; 
plot(100*(1-alpha),sr2_avg(:,bb,ss,nb,ne),'Marker','square','color','m'); hold on;
plot(100*(1-alpha),sr2_ICSI_model(:,bb,ss,nb,ne),'Marker','square','color','g'); hold on;
plot(100*(1-alpha),sr5_avg(:,bb,ss,nb,ne),'Marker','diamond','color','c'); hold on;
plot(100*(1-alpha),sr5_ICSI_model(:,bb,ss,nb,ne),'Marker','diamond','color','y'); hold on;
box on; grid on;
xlabel('Percentage of radiated AN energy (\%)')
ylabel('Secrecy rate (bit/channel use)')
%ylim([min(sr2_avg(:,bb,ss,nb,ne)-1,max(sr1_ICSI_model(:,bb,ss,nb,ne))+1])
legend('SDS - Simulation' , 'SDS - Modeling' , 'MF - Simulation' , 'MF - Modeling', 'OC - Simulation', 'OC - Modeling', 'location', 'bestoutside')
title(['Imperfect CSI = ', num2str(sigma_tilde_db(ss)), 'dB, U = ', num2str(U(bb)), ' SNR @B = ', num2str(snr_b(nb)), 'dB, SNR @E = ', num2str(snr_b(nb)), 'dB'])



%% 2. SR as a function of AN for different ICSI and a given model : model

% 2.1 For SDS
figure;
plot(100*(1-alpha),squeeze(sr1_ICSI_model(:,bb,1:4:end,nb,ne)),'Marker','o'); hold on; 
legend(cellstr(num2str(sigma_tilde_db(1:4:end).', 'Imperfect CSI = %-0.f dB')), 'location', 'bestoutside')
box on; grid on;
xlabel('Percentage of radiated AN energy (\%)')
ylabel('Secrecy rate (bit/channel use)')
title(['SDS Decoder, U = ', num2str(U(bb)), ', SNR @B = ', num2str(snr_b(nb)), 'dB, SNR @E = ', num2str(snr_b(nb)), 'dB'])


% 2.2 For MF
figure;
plot(100*(1-alpha),squeeze(sr2_ICSI_model(:,bb,1:4:end,nb,ne)),'Marker','o'); hold on; 
legend(cellstr(num2str(sigma_tilde_db(1:4:end).', 'Imperfect CSI = %-0.f dB')), 'location', 'bestoutside')
box on; grid on;
xlabel('Percentage of radiated AN energy (\%)')
ylabel('Secrecy rate (bit/channel use)')
title(['MF Decoder, U = ', num2str(U(bb)), ', SNR @B = ', num2str(snr_b(nb)), 'dB, SNR @E = ', num2str(snr_b(nb)), 'dB'])



% 2.3 For OC
figure;
plot(100*(1-alpha),squeeze(sr5_ICSI_model(:,bb,1:4:end,nb,ne)),'Marker','o'); hold on; 
legend(cellstr(num2str(sigma_tilde_db(1:4:end).', 'Imperfect CSI = %-0.f dB')), 'location', 'bestoutside')
box on; grid on;
xlabel('Percentage of radiated AN energy (\%)')
ylabel('Secrecy rate (bit/channel use)')
title(['OC Decoder, U = ', num2str(U(bb)), ', SNR @B = ', num2str(snr_b(nb)), 'dB, SNR @E = ', num2str(snr_b(nb)), 'dB'])



%% 3. Best alpha to inject
% 3.1 For SDS
alpha_opt1 = optimalAlphaICSI(U,snr_b,snr_e,sigma_tilde,"decod1");
figure; 
plot(sigma_tilde_db,100*(1-alpha_opt1),'-o');
legend(cellstr(num2str(U.', 'U = %-0.f')), 'location', 'bestoutside')
box on; grid on;
xlabel('CSI error (dB)')
ylabel('Optimal AN energy to inject (\%)')
title(['SDS Decoder, SNR @B = ', num2str(snr_b(nb)), 'dB, SNR @E = ', num2str(snr_b(nb)), 'dB'])

% 3.2 For MF
% Not determined since condition on sigma to have positive SR

% 3.3 For OC
alpha_opt5 = optimalAlphaICSI(U,snr_b,snr_e,sigma_tilde,"decod5");
figure; 
plot(sigma_tilde_db,100*(1-alpha_opt5),'-o');
legend(cellstr(num2str(U.', 'U = %-0.f')), 'location', 'bestoutside')
box on; grid on;
xlabel('CSI error (dB)')
ylabel('Optimal AN energy to inject (\%)')
title(['SDS Decoder, SNR @B = ', num2str(snr_b(nb)), 'dB, SNR @E = ', num2str(snr_b(nb)), 'dB'])


%% 4. Sigma max in order to have SR > 0 when Eve is noisy
snr_b  = [0:25]; %EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = [0:25 200];

sigma_b = 1./U.'./10.^(snr_b./10);    % expected noise energy @Bob
sigma_e = 1./U.'./10.^(snr_e./10);    % expected noise energy @Eve

% 4.1 SDS Decoder

sigma_max1 = ICSImaxSigmaEveNoisy(U,snr_b,snr_e,alpha,"decod1");
sigma_max1_db = 10*log10(sigma_max1);

% 4.1 MF Decoder

sigma_max2 = ICSImaxSigmaEveNoisy(U,snr_b,snr_e,alpha,"decod2");
sigma_max2_db = 10*log10(sigma_max2);

% 4.1 OC Decoder

sigma_max5 = ICSImaxSigmaEveNoisy(U,snr_b,snr_e,alpha,"decod5");
sigma_max5_db = 10*log10(sigma_max5);


%% 5. Guaranteeing SR = ∆ bit/channel use
SR_target = [0.2:.1:4];
SR_target_lin = 2.^SR_target.';

% 5.1 Max allowed sigma to ensure SR = ∆ when sigma_e = 0 (Eve noiseless)
% 5.1.1 SDS Decoder
sigma_max_noiseless1 = ICSImaxSigmaEveNoiseless(U,SR_target_lin,"decod1");
sigma_max_noiseless1_db = 10*log10(sigma_max_noiseless1);
figure;
plot(SR_target,sigma_max_noiseless1_db,'-o');
legend(cellstr(num2str(U.', 'U = %-0.f')), 'location', 'bestoutside')
box on; grid on;
ylabel('CSI error (dB)')
xlabel('Targeted SR (bit/channel use)')
xlim([min(SR_target) max(SR_target)])
title(['SDS Decoder, maximum allowed CSI error that can be made @A'])



% 5.1.2 MF Decoder
sigma_max_noiseless2 = ICSImaxSigmaEveNoiseless(U,SR_target_lin,"decod2");
sigma_max_noiseless2_db = 10*log10(sigma_max_noiseless2);
figure;
plot(SR_target,sigma_max_noiseless2_db,'-o');
legend(cellstr(num2str(U.', 'U = %-0.f')), 'location', 'bestoutside')
box on; grid on;
ylabel('CSI error (dB)')
xlabel('Targeted SR (bit/channel use)')
xlim([min(SR_target) max(SR_target)])
title(['MF Decoder, maximum allowed CSI error that can be made @A'])


% 5.1.3 OC Decoder
sigma_max_noiseless5 = ICSImaxSigmaEveNoiseless(U,SR_target_lin,"decod5");
sigma_max_noiseless5_db = 10*log10(sigma_max_noiseless5);
figure;
plot(SR_target,sigma_max_noiseless5_db,'-o');
legend(cellstr(num2str(U.', 'U = %-0.f')), 'location', 'bestoutside')
box on; grid on;
ylabel('CSI error (dB)')
xlabel('Targeted SR (bit/channel use)')
xlim([min(SR_target) max(SR_target)])
title(['OC Decoder, maximum allowed CSI error that can be made @A'])


% 5.2 Optimal alpha to inject to ensure SR = ∆ when sigma_e = 0 
% If Alice is aware of the error of CSI estimation, the real alpha_opt is
% given by:

nb_sigma = 36; % Number of CSI error from 0 to sigma_max
sigma_tilde1 = zeros(length(U),nb_sigma,size(sigma_max_noiseless1,1));
sigma_tilde2 = zeros(length(U),nb_sigma,size(sigma_max_noiseless1,1));
sigma_tilde5 = zeros(length(U),nb_sigma,size(sigma_max_noiseless1,1));

for aa = 1:size(sigma_max_noiseless1,1)
    for bb = 1:length(U)
        sigma_tilde1(bb,:,aa) = linspace(0,sigma_max_noiseless1(aa,bb),nb_sigma);
        sigma_tilde2(bb,:,aa) = linspace(0,sigma_max_noiseless2(aa,bb),nb_sigma);
        sigma_tilde5(bb,:,aa) = linspace(0,sigma_max_noiseless5(aa,bb),nb_sigma);
    end
end
    
sigma_tilde1_db = 10*log10(sigma_tilde1);
sigma_tilde2_db = 10*log10(sigma_tilde2);
sigma_tilde5_db = 10*log10(sigma_tilde5);
alpha_opt1_inf = zeros(length(U), nb_sigma-1,length(SR_target_lin));
alpha_opt2_inf = zeros(length(U), nb_sigma-1,length(SR_target_lin));

for aa = 1:length(U)
for bb = 1:nb_sigma-1
for cc = 1:length(SR_target_lin)
    A = (U(aa)+1)*(1-sigma_tilde1(aa,bb,cc));
    B = sigma_tilde1(aa,bb,cc)*(SR_target_lin(cc)-1);
    C = A+B;
    delta = 4*A^2*(SR_target_lin(cc)-1)^2 - 4*A*(-(SR_target_lin(cc)-1)*C-B);
    alpha_opt1_inf(aa,bb,cc) = (-2*A*(SR_target_lin(cc)-1) + sqrt(delta))/(2*A);
    if alpha_opt1_inf(aa,bb,cc) > 1
        alpha_opt1_inf(aa,bb,cc) = NaN;
    end
    
    BB = U(aa)^2+3*U(aa)+3;
    X = sigma_tilde2(aa,bb,cc)*(SR_target_lin(cc)*BB +U(aa)*(U(aa)+1)) - U(aa)*(U(aa)+1);
    Y = sigma_tilde2(aa,bb,cc)*(SR_target_lin(cc)*U(aa) - SR_target_lin(cc)*BB - U(aa)*(U(aa)+2)) + U(aa)*(U(aa)+1);
    Z = U(aa)+SR_target_lin(cc)*BB;
    T = U(aa)*(1-SR_target_lin(cc));
    alpha_opt2_inf(aa,bb,cc) = (X*T-sqrt(X^2*T^2+Z*X*T*(sigma_tilde2(aa,bb,cc)*Z+Y) ))/(Z*X);
    if alpha_opt2_inf(aa,bb,cc) > 1
        alpha_opt2_inf(aa,bb,cc) = NaN;
    end
end
end
end
alpha_opt5_inf  = alpha_opt1_inf;

% Since Alice is not aware of the error shes doing, she wants to inject
% alpha_opt considering that there is no CSI error. Therefore, she will
% inject another amount of AN energy. 

% % SDS Decoder
% A1 = SR_target_lin - 1;
% alpha_opt1_estimated = sqrt(A1.^2+A1)-A1;                    % alpha to inject in order to minimize Bob SNR expression when sigma_e = infinite and model = model1
% SNRb1 = 10*log10(sqrt(A1.^2+A1)./((U+1).*(-2*A1.^2-2*A1+sqrt(A1.^2+A1).*(2*A1+1))));
% 
% % MF Decoder
% B2 = U.*(sr_lin-1);
% A2 = sr_lin.*(U+1).*(U+3)-B2;
% alpha_opt2_estimated = (-2*B2 + sqrt(4.*A2.*B2+4.*B2.^2))./(2.*A2);                     % alpha to inject in order to minimize Bob SNR expression when sigma_e = infinite and model = model2
% SNRb2 = 10*log10((alpha2.*A2 + B2)./((-alpha2.^2+alpha2).*(U+1).*U));
% 
% % OC Decoder
% alpha_opt5_estimated = alpha_opt1_estimated;                % alpha to inject in order to minimize Bob SNR expression when sigma_e = infinite and model = model5
% SNRb5 = SNRb1;
% 
% 






sr1_tmp = find(SR_target == 2);
sr2_tmp = find(SR_target == 4);
sr3_tmp = find(SR_target == 6);

% SDS plot
figure; 
plot(squeeze(sigma_tilde1_db(1,1:end-1,sr1_tmp)),100*(1-squeeze(alpha_opt1_inf(1,:,sr1_tmp))),'-o'); hold on;
plot(squeeze(sigma_tilde1_db(2,1:end-1,sr1_tmp)),100*(1-squeeze(alpha_opt1_inf(2,:,sr1_tmp))),'-o'); hold on;
legendCell1 = cellstr(num2str(U.', 'U = %-0.f, SR = 2 bits/ch use'));

plot(squeeze(sigma_tilde1_db(1,1:end-1,sr2_tmp)),100*(1-squeeze(alpha_opt1_inf(1,:,sr2_tmp))),'marker','diamond'); hold on;
plot(squeeze(sigma_tilde1_db(2,1:end-1,sr2_tmp)),100*(1-squeeze(alpha_opt1_inf(2,:,sr2_tmp))),'marker','diamond'); hold on;
legendCell2 = cellstr(num2str(U.', 'U = %-0.f, SR = 4 bits/ch use'));

plot(squeeze(sigma_tilde1_db(1,1:end-1,sr3)),100*(1-squeeze(alpha_opt1_inf(1,:,sr3))),'marker','square'); hold on;
plot(squeeze(sigma_tilde1_db(2,1:end-1,sr3)),100*(1-squeeze(alpha_opt1_inf(2,:,sr3))),'marker','square');
legendCell3 = cellstr(num2str(U.', 'U = %-0.f, SR = 6 bits/ch use'));

legendCell =[legendCell1;legendCell2;legendCell3];
legend(legendCell, 'location', 'bestoutside')
box on; grid on;
xlabel('CSI error (dB)')
ylabel('Optimal AN energy to inject (\%)')
title(['SDS Decoder, Eve noiseless'])

% MF plot
figure; 
plot(squeeze(sigma_tilde2_db(1,1:end-1,sr1_tmp)),100*(1-squeeze(alpha_opt2_inf(1,:,sr1_tmp))),'-o'); hold on;
plot(squeeze(sigma_tilde2_db(2,1:end-1,sr1_tmp)),100*(1-squeeze(alpha_opt2_inf(2,:,sr1_tmp))),'-o'); hold on;
legendCell1 = cellstr(num2str(U.', 'U = %-0.f, SR = 2 bits/ch use'));

plot(squeeze(sigma_tilde2_db(1,1:end-1,sr2_tmp)),100*(1-squeeze(alpha_opt2_inf(1,:,sr2_tmp))),'marker','diamond'); hold on;
plot(squeeze(sigma_tilde2_db(2,1:end-1,sr2_tmp)),100*(1-squeeze(alpha_opt2_inf(2,:,sr2_tmp))),'marker','diamond'); hold on;
legendCell2 = cellstr(num2str(U.', 'U = %-0.f, SR = 4 bits/ch use'));

plot(squeeze(sigma_tilde2_db(1,1:end-1,sr3_tmp)),100*(1-squeeze(alpha_opt2_inf(1,:,sr3_tmp))),'marker','square'); hold on;
plot(squeeze(sigma_tilde2_db(2,1:end-1,sr3_tmp)),100*(1-squeeze(alpha_opt2_inf(2,:,sr3_tmp))),'marker','square');
legendCell3 = cellstr(num2str(U.', 'U = %-0.f, SR = 6 bits/ch use'));

legendCell =[legendCell1;legendCell2;legendCell3];
legend(legendCell, 'location', 'bestoutside')
box on; grid on;
xlabel('CSI error (dB)')
ylabel('Optimal AN energy to inject (\%)')
title(['MF Decoder, Eve noiseless'])


% OC plot
figure; 
plot(squeeze(sigma_tilde5_db(1,1:end-1,sr1)),100*(1-squeeze(alpha_opt5_inf(1,:,sr1))),'-o'); hold on;
plot(squeeze(sigma_tilde5_db(2,1:end-1,sr1)),100*(1-squeeze(alpha_opt5_inf(2,:,sr1))),'-o'); hold on;
legendCell1 = cellstr(num2str(U.', 'U = %-0.f, SR = 2 bits/ch use'));

plot(squeeze(sigma_tilde5_db(1,1:end-1,sr2)),100*(1-squeeze(alpha_opt5_inf(1,:,sr2))),'marker','diamond'); hold on;
plot(squeeze(sigma_tilde5_db(2,1:end-1,sr2)),100*(1-squeeze(alpha_opt5_inf(2,:,sr2))),'marker','diamond'); hold on;
legendCell2 = cellstr(num2str(U.', 'U = %-0.f, SR = 4 bits/ch use'));

plot(squeeze(sigma_tilde5_db(1,1:end-1,sr3)),100*(1-squeeze(alpha_opt5_inf(1,:,sr3))),'marker','square'); hold on;
plot(squeeze(sigma_tilde5_db(2,1:end-1,sr3)),100*(1-squeeze(alpha_opt5_inf(2,:,sr3))),'marker','square');
legendCell3 = cellstr(num2str(U.', 'U = %-0.f, SR = 6 bits/ch use'));

legendCell =[legendCell1;legendCell2;legendCell3];
legend(legendCell, 'location', 'bestoutside')
box on; grid on;
xlabel('CSI error (dB)')
ylabel('Optimal AN energy to inject (\%)')
title(['OC Decoder, Eve noiseless'])


% 5.3 Required SNR at Bob to ensure SR = ∆ when sigma_e = 0 
% 5.3.1 SDS Decoder cfr recto p16
alpha1 = alpha_opt1_inf;
snr_b_inf1 = zeros(length(U), nb_sigma-1,length(SR_target_lin));
for aa = 1:length(U)
for bb = 1:nb_sigma-1
for cc = 1:length(SR_target_lin)
    snr_b_inf1(aa,bb,cc) = 10*log10((alpha1(aa,bb,cc) + (SR_target_lin(cc)-1) ) / ...
                           (-alpha1(aa,bb,cc)^2*(U(aa)+1)*(1-sigma_tilde1(aa,bb,cc)) + ...
                           alpha1(aa,bb,cc)*( (U(aa) + 1)*(1-sigma_tilde1(aa,bb,cc)) + SR_target_lin(cc)*sigma_tilde1(aa,bb,cc) -sigma_tilde1(aa,bb,cc) ) + ...
                           sigma_tilde1(aa,bb,cc)*(1-SR_target_lin(cc)) )) ;                        
end
end
end

% 5.3.2 MF Decoder
alpha2 = alpha_opt2_inf;
snr_b_inf2 = zeros(length(U), nb_sigma-1,length(SR_target_lin));
for aa = 1:length(U)
for bb = 1:nb_sigma-1
for cc = 1:length(SR_target_lin)
    B = U(aa)^2 + 3*U(aa) + 3;
    snr_b_inf2(aa,bb,cc) = 10*log10((alpha2(aa,bb,cc)*(U(aa) + SR_target_lin(cc)*B ) + U(aa)*(SR_target_lin(cc)-1) ) / ...
                           ( alpha2(aa,bb,cc)^2*( SR_target_lin(cc)*B*sigma_tilde2(aa,bb,cc) -  ...
                           U(aa)*(U(aa)+1)*(1-sigma_tilde2(aa,bb,cc)))  + ...
                           alpha2(aa,bb,cc)*(SR_target_lin(cc)*U(aa)*sigma_tilde2(aa,bb,cc) - SR_target_lin(cc)*sigma_tilde2(aa,bb,cc)*B + ...
                           U(aa)*(U(aa)+1)*(1 - sigma_tilde2(aa,bb,cc)) - U(aa)*sigma_tilde2(aa,bb,cc)) + ...
                           U(aa)*sigma_tilde2(aa,bb,cc)*(1-SR_target_lin(cc)) )) ;                        
end
end
end




% 5.3.3 OC Decoder
snr_b_inf5 = snr_b_inf1;



figure;
plot(SR_target,squeeze(snr_b_inf1(:,1,:)),'--'); hold on;
legendCell1 = cellstr(num2str(U.', 'U = %-0.f, $\sigma$ = -\infty dB'));
plot(SR_target,squeeze(snr_b_inf1(:,5,:)),'-o');
legendCell2 = cellstr(num2str(U.', 'U = %-0.f, $\sigma$ = -\infty dB'));


figure;
plot(SR_target,squeeze(snr_b_inf2(:,1,:)),'--'); hold on;
legendCell1 = cellstr(num2str(U.', 'U = %-0.f, $\sigma$ = -\infty dB'));
plot(SR_target,squeeze(snr_b_inf2(:,30,:)),'-o');


%% 6 Outage probability
s_outage1 = sort(squeeze(max(sr1,[],2)),1); % Max over AN dimension (i.e. the 2nd dimension of the sr1 matrix) --> dim of tmp1 = nb_iter x nb_BOR x nb_ICSI
figure;
plot(s_outage1(:,:,1),linspace(0,1,nb_run))

%% FAIRE DES SWITCH CASE POUR LES PARTIES DU CODE: 1 SIMU ; 2 MODELING & POST PROCESS






