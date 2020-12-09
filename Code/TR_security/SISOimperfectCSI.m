% LAST MODIF: 
% should uncomment plot section, change the legends and include the new
% model (model 5) where Eve knows her own channel He. 
% DATE: 28.07.2020


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
alpha_step = 2;             % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;  
sigma_tilde = 0.3868;    % sigma_tilde := percentage of CSI error --> = 0 : perfect CSI @Alice , = 1: 100% error

% Communication parameters
Q = 128;
U = [4];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
snr_b  = 20; %EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = 20; %EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance


if (length(snr_b) == 1 && length(snr_e) ==1)
%% Matrices preallocation
e_an_TX             = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));

e_sym_decod1_b      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_noise_decod1_b    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_an_decod1_b       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_denom_decod1_b    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));

e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
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
e_sym_decod5_e      = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_noise_decod5_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_an_decod5_e       = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));
e_denom_decod5_e    = zeros(nb_run,length(alpha),length(U),length(sigma_tilde));


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




% Noise symbol
[noise_b, ~ ] = addNoise(sym_b , snr_b, energy(sym_precoded_TX+an_TX));     %addNoise(sym_b , snr_b, energy(sym_RX_b)); 
[noise_e, e_noise_e ] = addNoise(sym_e , snr_e, energy(sym_precoded_TX+an_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   

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
e_sym_decod1_b(iter,aa,bb,ss)      = energy(sym_decod1_b);
e_noise_decod1_b(iter,aa,bb,ss)    = energy(noise_decod1_b);
e_an_decod1_b(iter,aa,bb,ss)       = energy(an_decod1_b);
e_denom_decod1_b(iter,aa,bb,ss)    = energy(noise_decod1_b + an_decod1_b);        % energy of the sinr denominator for decoder 1 @Eve


% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,ss)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,ss)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,ss)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,ss)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,aa,bb,ss)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb,ss)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb,ss)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb,ss)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve
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
e_sym_decod5_e(iter,aa,bb,ss)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,aa,bb,ss)   = energy(noise_decod5_e);
e_an_decod5_e(iter,aa,bb,ss)      = energy(an_decod5_e);
e_denom_decod5_e(iter,aa,bb,ss)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve


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

sinr1_PCSI_model_b = sinrModeling(alpha,U,snr_b,snr_e,1,1,1,"bob_decod1");
sinr1_ICSI_model_b = sinrModelingICSI(alpha,U,snr_b,snr_e,sigma_tilde,"bob_decod1");
sinr1_ICSI_model_e = sinrModelingICSI(alpha,U,snr_b,snr_e,sigma_tilde,"eve_decod1");
sinr2_ICSI_model_e = sinrModelingICSI(alpha,U,snr_b,snr_e,sigma_tilde,"eve_decod2");
sinr5_ICSI_model_e = sinrModelingICSI(alpha,U,snr_b,snr_e,sigma_tilde,"eve_decod5");

sr1_PCSI_model = secrecyCapacity(sinr1_PCSI_model_b,sinr1_ICSI_model_e);
sr1_ICSI_model = secrecyCapacity(sinr1_ICSI_model_b,sinr1_ICSI_model_e);
sr2_ICSI_model = secrecyCapacity(sinr1_ICSI_model_b,sinr2_ICSI_model_e);
sr5_ICSI_model = secrecyCapacity(sinr1_ICSI_model_b,sinr5_ICSI_model_e);
% 

figure;
plot(100*(1-alpha),sr2_avg,'Marker','o'); hold on;
plot(100*(1-alpha),sr2_ICSI_model,'Marker','<'); hold on; 
% plot(100*(1-alpha),sr2_avg,'Marker','square','color','m'); hold on;
% plot(100*(1-alpha),sr2_model,'Marker','square','color','g'); hold on;
% plot(100*(1-alpha),sr5_avg,'Marker','diamond','color','c'); hold on;
% plot(100*(1-alpha),sr5_model,'Marker','diamond','color','y'); hold on;
box on; grid on;
xlabel('Percentage of radiated AN energy (\%)')
ylabel('Secrecy rate (bit/channel use)')
ylim([-2 3])
legend('Same decoder - Simulation' , 'Same decoder - Modeling' , 'Matched filtering - Simulation' , 'Matched filtering - Modeling', 'Own channel - Simulation', 'Own channel - Modeling', 'location', 'best')


end



%% FAIRE DES SWITCH CASE POUR LES PARTIES DU CODE: 1 SIMU ; 2 MODELING & POST PROCESS






