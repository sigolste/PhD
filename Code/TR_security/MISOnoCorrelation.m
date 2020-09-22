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
alpha_step = 1;             % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

Na = 1:5;           % Number of antennas @Alice
Nb = 1;             % Number of antennas @Bob
Ne = 1;             % Number of antennas @Eve
% Communication parameters
Q = 8;
U = [4];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 15;        % energy per bit over noise psd @Bob - dB
EbN0_e = [15];      % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance


%% Matrices preallocation
e_an_TX             = zeros(nb_run,length(alpha),length(U),length(Na));

e_sym_decod1_b      = zeros(nb_run,length(alpha),length(U),length(Na));
e_noise_decod1_b    = zeros(nb_run,length(alpha),length(U),length(Na));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),length(Na));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),length(Na));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),length(Na));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),length(Na));

for iter = 1:nb_run
for bb =1:length(U)
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run

% Channel generation
% Rayleigh channel - independant subcarrier fading - unit power per
% subcarrier : CN(0,I_Q) - circulary complex random normal variable.
for nt = 1:length(Na)
Hb_TX = channelRayleighMISO(Na(nt),Q, mu , sigma);  
He_TX = channelRayleighMISO(Na(nt),Q, mu , sigma);
Hb_RX = conj(Hb_TX);
He_RX = conj(He_TX);


% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrixMISO(Q,Na(nt),N(bb),U(bb)); % MISO --> size of spreading matrix: Q.Na x N
                                                                             % beamforming gain: S^H S = Na*I

%% Encoder

% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % QxQ
sym_precoded = Hb_TX*sym_spread; % QxQxNa, not weighted 

% AN generation
an = generateANMISO(Na(nt),Hb_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd"); % QxNa, not weighted

for aa = 1:length(alpha)
sym_precoded_TX = sqrt(alpha(aa))*sym_precoded;     % weighted
an_TX = sqrt(1-alpha(aa))*an;                       % weighted
%e_an_TX(iter,aa) = energy(an_TX);

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

sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;

noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;

an_decod1_e = decod1*an_e;


%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,aa,bb,nt) = energy(sym_decod1_b);
e_noise_decod1_b(iter,aa,bb,nt) = energy(noise_decod1_b);

% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,nt)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,nt)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,nt)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,nt)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve



end
end
waitbar(iter / nb_run)
end
end

%% Post processing
% Energy averaging over all realisations
e_avg_an                = squeeze(mean(e_an_TX,1));

e_avg_sym_decod1_b      = squeeze(mean(e_sym_decod1_b,1));  
e_avg_noise_decod1_b    = squeeze(mean(e_noise_decod1_b,1));

e_avg_sym_decod1_e      = squeeze(mean(e_sym_decod1_e,1));
e_avg_noise_decod1_e    = squeeze(mean(e_noise_decod1_e,1));
e_avg_an_decod1_e       = squeeze(mean(e_an_decod1_e,1));
e_avg_denom_decod1_e    = squeeze(mean(e_denom_decod1_e,1));


% instantaneous SINRs
sinr1_b = e_sym_decod1_b./e_noise_decod1_b;
sinr1_e = e_sym_decod1_e./e_denom_decod1_e;


% instantaneous Secrecy capacity
sr1 = secrecyCapacity(sinr1_b,sinr1_e);


% ergodic Secrecy capacity
sr1_avg = squeeze(mean(sr1,1));


sinr1_model_b = sinrModeling(alpha,U,snr_b,snr_e,Na,Nb,Ne,"bob_MISO_decorrelated");
sinr1_model_e = sinrModeling(alpha,U,snr_b,snr_e,1,1,1,"eve_decod1");

sr1_model = secrecyCapacity(sinr1_model_b,sinr1_model_e);
% Comparaison simu vs models - SR curves
figure;
plot(100*alpha,sr1_avg,'Marker','o','color','b'); hold on;
plot(100*alpha,sr1_model,'Marker','o','color','r'); hold on;
close(h);