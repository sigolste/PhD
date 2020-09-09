% LAST MODIF: 
% should uncomment plot section, change the legends and include the new
% model (model 5) where Eve knows her own channel He. 
% DATE: 28.07.2020


clear all;
close all;

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
alpha_step = .25;             % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

% Communication parameters
Q = 64;
U = [2 4 8 16];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 10; % energy per bit over noise psd @Bob - dB
EbN0_e = [5]; % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance


if (length(snr_b) == 1 && length(snr_e) ==1)
%% Matrices preallocation
e_an_TX             = zeros(nb_run,length(alpha),length(U));

e_sym_decod1_b      = zeros(nb_run,length(alpha),length(U));
e_noise_decod1_b    = zeros(nb_run,length(alpha),length(U));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U));

e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U));

e_sym_decod3_e      = zeros(nb_run,length(alpha),length(U));
e_noise_decod3_e    = zeros(nb_run,length(alpha),length(U));
e_an_decod3_e       = zeros(nb_run,length(alpha),length(U));
e_denom_decod3_e    = zeros(nb_run,length(alpha),length(U));

e_sym_decod4_e      = zeros(nb_run,length(alpha),length(U));
e_noise_decod4_e    = zeros(nb_run,length(alpha),length(U));
e_an_decod4_e       = zeros(nb_run,length(alpha),length(U));
e_denom_decod4_e    = zeros(nb_run,length(alpha),length(U));

e_sym_decod5_e      = zeros(nb_run,length(alpha),length(U));
e_noise_decod5_e    = zeros(nb_run,length(alpha),length(U));
e_an_decod5_e       = zeros(nb_run,length(alpha),length(U));
e_denom_decod5_e    = zeros(nb_run,length(alpha),length(U));


%% Mainloop
for iter = 1:nb_run
for bb =1:length(U)
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run

% Channel generation
% Rayleigh channel - independant subcarrier fading - unit power per
% subcarrier : CN(0,I_Q) - circulary complex random normal variable.
Hb_TX = channelRayleigh(Q, mu , sigma);  
He_TX = channelRayleigh(Q, mu , sigma);
Hb_RX = ctranspose(Hb_TX.');
He_RX = ctranspose(He_TX.');



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));

%% Encoder

% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % Qx1
sym_precoded = Hb_TX*sym_spread; % Qx1, not weighted 

% AN generation
an = generateAN(Hb_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd"); % Qx1, not weighted

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
an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob



%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_RX*He_TX;                       % matched filter
decod3 = matrix_despread*(Hb_RX/He_RX);                     % AN killer
gamma_E = (He_RX*Hb_TX)*matrix_spread;
gamma_EH = ctranspose(gamma_E);
decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an) + e_noise_e*eye(Q) );   % LMMSE
decod5 = matrix_despread*He_TX;                                   % Only He known by Eve


sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_decod2_e = decod2*sym_e;
sym_decod3_e = decod3*sym_e;
sym_decod4_e = decod4*sym_e;
sym_decod5_e = decod5*sym_e;

noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;
noise_decod2_e = decod2*noise_e;
noise_decod3_e = decod3*noise_e;
noise_decod4_e = decod4*noise_e;
noise_decod5_e = decod5*noise_e;

an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;
an_decod3_e = decod3*an_e;
an_decod4_e = decod4*an_e;
an_decod5_e = decod5*an_e;


% sinr_no_decod_e(iter,aa,:) = U*diag((ctranspose(He_equiv)*He_equiv)/(ctranspose(an_decod2_e+noise_decod2_e)*(an_decod2_e+noise_decod2_e)));


%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,aa,bb) = energy(sym_decod1_b);
e_noise_decod1_b(iter,aa,bb) = energy(noise_decod1_b);

% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,aa,bb)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

% @ Eve : decod 3
e_sym_decod3_e(iter,aa,bb)     = energy(sym_decod3_e);
e_noise_decod3_e(iter,aa,bb)   = energy(noise_decod3_e);
e_an_decod3_e(iter,aa,bb)      = energy(an_decod3_e);
e_denom_decod3_e(iter,aa,bb)   = energy(noise_decod3_e + an_decod3_e);        % energy of the sinr denominator for decoder 3 @Eve

% @ Eve : decod 4
e_sym_decod4_e(iter,aa,bb)     = energy(sym_decod4_e);
e_noise_decod4_e(iter,aa,bb)   = energy(noise_decod4_e);
e_an_decod4_e(iter,aa,bb)      = energy(an_decod4_e);
e_denom_decod4_e(iter,aa,bb)   = energy(noise_decod4_e + an_decod4_e);        % energy of the sinr denominator for decoder 4 @Eve

% @ Eve : decod 5
e_sym_decod5_e(iter,aa,bb)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,aa,bb)   = energy(noise_decod5_e);
e_an_decod5_e(iter,aa,bb)      = energy(an_decod5_e);
e_denom_decod5_e(iter,aa,bb)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve


end
waitbar(iter / nb_run)
end
end
close(h)

%% Post processing
% Energy averaging over all realisations
e_avg_an                = squeeze(mean(e_an_TX,1));

e_avg_sym_decod1_b      = squeeze(mean(e_sym_decod1_b,1));  
e_avg_noise_decod1_b    = squeeze(mean(e_noise_decod1_b,1));

e_avg_sym_decod1_e      = squeeze(mean(e_sym_decod1_e,1));
e_avg_noise_decod1_e    = squeeze(mean(e_noise_decod1_e,1));
e_avg_an_decod1_e       = squeeze(mean(e_an_decod1_e,1));
e_avg_denom_decod1_e    = squeeze(mean(e_denom_decod1_e,1));

e_avg_sym_decod2_e      = squeeze(mean(e_sym_decod2_e,1));
e_avg_noise_decod2_e    = squeeze(mean(e_noise_decod2_e,1));
e_avg_an_decod2_e       = squeeze(mean(e_an_decod2_e,1));
e_avg_denom_decod2_e    = squeeze(mean(e_denom_decod2_e,1));

e_avg_sym_decod3_e      = squeeze(mean(e_sym_decod3_e,1));
e_avg_noise_decod3_e    = squeeze(mean(e_noise_decod3_e,1));
e_avg_an_decod3_e       = squeeze(mean(e_an_decod3_e,1));
e_avg_denom_decod3_e    = squeeze(mean(e_denom_decod3_e,1));

e_avg_sym_decod4_e      = squeeze(mean(e_sym_decod4_e,1));
e_avg_noise_decod4_e    = squeeze(mean(e_noise_decod4_e,1));
e_avg_an_decod4_e       = squeeze(mean(e_an_decod4_e,1));
e_avg_denom_decod4_e    = squeeze(mean(e_denom_decod4_e,1));

e_avg_sym_decod5_e      = squeeze(mean(e_sym_decod5_e,1));
e_avg_noise_decod5_e    = squeeze(mean(e_noise_decod5_e,1));
e_avg_an_decod5_e       = squeeze(mean(e_an_decod5_e,1));
e_avg_denom_decod5_e    = squeeze(mean(e_denom_decod5_e,1));

% instantaneous SINRs
sinr1_b = e_sym_decod1_b./e_noise_decod1_b;
sinr1_e = e_sym_decod1_e./e_denom_decod1_e;
sinr2_e = e_sym_decod2_e./e_denom_decod2_e;
sinr3_e = e_sym_decod3_e./e_denom_decod3_e;
sinr4_e = e_sym_decod4_e./e_denom_decod4_e;
sinr5_e = e_sym_decod5_e./e_denom_decod5_e;


% instantaneous Secrecy capacity
sr1 = secrecyCapacity(sinr1_b,sinr1_e);
sr2 = secrecyCapacity(sinr1_b,sinr2_e);
sr3 = secrecyCapacity(sinr1_b,sinr3_e);
sr4 = secrecyCapacity(sinr1_b,sinr4_e);
sr5 = secrecyCapacity(sinr1_b,sinr5_e);


% ergodic Secrecy capacity
sr1_avg = squeeze(mean(sr1,1));
sr2_avg = squeeze(mean(sr2,1));
sr3_avg = squeeze(mean(sr3,1));
sr4_avg = squeeze(mean(sr4,1));
sr5_avg = squeeze(mean(sr5,1));


%% SINR modeling:
sinr1_model_b = sinrModeling(alpha,U,snr_b,snr_e,"bob_decod1");
sinr1_model_e = sinrModeling(alpha,U,snr_b,snr_e,"eve_decod1");
sinr2_model_e = sinrModeling(alpha,U,snr_b,snr_e,"eve_decod2");
sinr5_model_e = sinrModeling(alpha,U,snr_b,snr_e,"eve_decod5");


sr1_model = secrecyCapacity(sinr1_model_b,sinr1_model_e);
sr2_model = secrecyCapacity(sinr1_model_b,sinr2_model_e);
sr5_model = secrecyCapacity(sinr1_model_b,sinr5_model_e);


%% Optimal amount of data to inject and max SR

% Optimal alpha decod1
alpha1_model_opt = optimalAlpha(U,snr_b,snr_e,"model1");

% Optimal alpha decod2
alpha2_model_opt = optimalAlpha(U,snr_b,snr_e,"model2");

% Optimal alpha decod2
alpha5_model_opt = optimalAlpha(U,snr_b,snr_e,"model5");

% Max theoretical value of sr derived from optimal alpha
sr1_model_max = maxThSr(U,sr1_model,alpha1_model_opt,alpha_step);
sr2_model_max = maxThSr(U,sr2_model,alpha2_model_opt,alpha_step);
sr5_model_max = maxThSr(U,sr5_model,alpha5_model_opt,alpha_step);



% Max simu value of sr and corresponding optimal alpha
[sr1_simu_max , alpha1_simu_opt] = maxSimuSr(sr1_avg,alpha_step);
[sr2_simu_max , alpha2_simu_opt] = maxSimuSr(sr2_avg,alpha_step);
[sr5_simu_max , alpha5_simu_opt] = maxSimuSr(sr5_avg,alpha_step);


%% Plot section
if length(U) == 1
% Comparaison simu vs models - SR curves
figure;
plot(100*alpha,sr1_avg,'Marker','o','color','b'); hold on;
plot(100*alpha,sr2_avg,'Marker','diamond','color','r'); hold on;
plot(100*alpha,sr3_avg,'Marker','square','color','m'); hold on;
plot(100*alpha,sr4_avg,'Marker','v','color','g'); hold on;
plot(100*alpha,sr5_avg,'Marker','v','color','c'); hold on;

box on; grid on;
xlabel('Percentage of energy radiated dedicated for data (\%)')
ylabel('Secrecy rate (bit/channel use)')
legend('Same decoder - Simulation' , 'Matched filtering - Simulation' , 'AN killer - Simulation','LMMSE - Simulation', 'Own channel - Simulation', 'location', 'best')

% Comparaison models - SR curves
figure;
plot(100*alpha,sr1_avg,'Marker','o','color','b'); hold on;
plot(100*alpha,sr1_model,'Marker','o','color','r'); hold on;
plot(100*alpha,sr2_avg,'Marker','square','color','m'); hold on;
plot(100*alpha,sr2_model,'Marker','square','color','g'); hold on;
plot(100*alpha,sr5_avg,'Marker','diamond','color','c'); hold on;
plot(100*alpha,sr5_model,'Marker','diamond','color','y'); hold on;
box on; grid on;
xlabel('Percentage of energy radiated dedicated for data (\%)')
ylabel('Secrecy rate (bit/channel use)')
legend('Same decoder - Simulation' , 'Same decoder - Modeling' , 'Matched filtering - Simulation' , 'Matched filtering - Modeling', 'Own channel - Simulation', 'Own channel - Modeling', 'location', 'best')


else
% Comparaison simu vs models - maxima of SR
% figure;
% yyaxis left
% plot(U,100*alpha1_simu_opt,'Marker','square') ; hold on;
% plot(U,100*alpha1_model_opt,'Marker','square') ; hold on;
% plot(U,100*alpha2_simu_opt,'Marker','v') ; hold on;
% plot(U,100*alpha2_model_opt,'Marker','v') ; hold on;
% plot(U,100*alpha5_simu_opt,'Marker','>') ; hold on;
% plot(U,100*alpha5_model_opt,'Marker','>') ; 
% ylabel('Percentage of energy radiated dedicated for data (\%)')
% ylim([0 100])
% 
% yyaxis right
% plot(U,sr1_simu_max,'Marker','o') ; hold on;
% plot(U,sr1_model_max,'Marker','o') ; hold on;
% plot(U,sr2_simu_max,'Marker','<') ; hold on;
% plot(U,sr2_model_max,'Marker','<') ; hold on;
% plot(U,sr5_simu_max,'Marker','diamond') ; hold on;
% plot(U,sr5_model_max,'Marker','diamond') 
% ylabel('Secrecy rate (bit/channel use)')
% ylim([0 max(sr1_simu_max)])
% 
% xlabel('BOR')
% legend('Same decoder - Simulation','Same decoder - Modeling', ...
%     'Matched filtering - Simulation', 'Matched filtering - Modeling',...
%     'Own channel - Simulation', 'Own channel - Modeling',...
%     'Same decoder - Simulation','Same decoder - Modeling', ...
%     'Matched filtering - Simulation' , 'Matched filtering - Modeling',...
%     'Own channel - Simulation', 'Own channel - Modeling')      
%         


figure;
subplot(1,2,1)
plot(U,100*alpha1_simu_opt,'Marker','square') ; hold on;
plot(U,100*alpha1_model_opt,'Marker','square') ; hold on;
plot(U,100*alpha2_simu_opt,'Marker','v') ; hold on;
plot(U,100*alpha2_model_opt,'Marker','v') ; hold on;
plot(U,100*alpha5_simu_opt,'Marker','>') ; hold on;
plot(U,100*alpha5_model_opt,'Marker','>') ; 
ylabel('Percentage of energy radiated dedicated for data (\%)')
xlabel('BOR')
box on; grid on;
ylim([min(100*alpha2_model_opt) max(100*alpha1_simu_opt)])
legend('Same decoder - Simulation','Same decoder - Modeling', ...
    'Matched filtering - Simulation' , 'Matched filtering - Modeling',...
    'Own channel - Simulation', 'Own channel - Modeling')


subplot(1,2,2)
plot(U,sr1_simu_max,'Marker','o') ; hold on;
plot(U,sr1_model_max,'Marker','o') ; hold on;
plot(U,sr2_simu_max,'Marker','<') ; hold on;
plot(U,sr2_model_max,'Marker','<') ; hold on;
plot(U,sr5_simu_max,'Marker','diamond') ; hold on;
plot(U,sr5_model_max,'Marker','diamond') 
ylabel('Secrecy rate (bit/channel use)')
ylim([min(sr2_simu_max) max(sr1_model_max)])
xlabel('BOR')
box on; grid on;
legend('Same decoder - Simulation','Same decoder - Modeling', ...
    'Matched filtering - Simulation', 'Matched filtering - Modeling',...
    'Own channel - Simulation', 'Own channel - Modeling')
          
        


end        




else
SRtargetted = 0.1;
% SNR calibration @ Bob for a targetted SR and different SNR levels @ Eve
[snr_b_calib_model1 , alpha_calib_model1] = snrBobCalibration(SRtargetted,U,snr_e,alpha,alpha_step,"model1");
[snr_b_calib_model2 , alpha_calib_model2] = snrBobCalibration(SRtargetted,U,snr_e,alpha,alpha_step,"model2");
[snr_b_calib_model5 , alpha_calib_model5] = snrBobCalibration(SRtargetted,U,snr_e,alpha,alpha_step,"model5");

end


%% FAIRE DES SWITCH CASE POUR LES PARTIES DU CODE: 1 SIMU ; 2 MODELING & POST PROCESS






