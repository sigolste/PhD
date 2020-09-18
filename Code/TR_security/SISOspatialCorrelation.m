%**************************************************************************
%
%   This script simulates an implementation of a secure Time Reversal
%   scheme in the frequency domain using OFDM. 
%
%   SR computed after despreading
%
%   Frequency domain secure TR in a MISO configuration    
%   Complex Envelope (CE) simulations
%
%   OFDM to mitigate MPC effects
%
%   Artificial Noise added to secure the data transmission at Bob
%
%   Simulation of the effect of spatial correlation between Eve and Bob on
%   the SR. 
%   
%	Example code to call the function corr_spatial to generates 
%   nb_subca independent Rayleigh channel exhibiting spatial correlation 
%   assuming a three-dimensional Rayleigh channel
%
%   The spatial dependence of the correlation is given by eq. 24
%   from Ph. De Doncker, "Spatial Correlation Functions for Fields in Three
%   -Dimensional Rayleigh Channels," PIER, Vol. 40, 55-69, 2003
%
%   Inputs:
%   Q: number of sub-carrier (or frequency samples within the
%   frequency bandwidth)
%   U: Back off rate
%   b_subca: bandwidth of each subcarrier (of frequency separation between 
%   frequency samples)
%   r: vector of the considered distances away from Bob
%   fc: carrier frequency
%   nb_run: number of samples of each sub-channel 
%
%   Outputs:
%   H: channel matrix whose size is (nb_run x Q x nb_distance)

%
%
%   Code started: 9 September 2020
%   Last update: 9 September 2020
%
%
%   by Sidney Golstein & Julien Sarrazin
%
%**************************************************************************




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
nb_run = 10;                               % number of experiments
nb_pos = 80;                                % Number of spatial positions between Bob and Eve
fc = 2e9 ;                                  % Carrier frequency
c = 3e8;                                    % Light speed
b_subcar = .1e6 ;                           % Bandwidth of each subcarrier
r = linspace( 0 , 2*c/fc , nb_pos ) ;       % Eve's distance from Bob . Distance from 0 to 1 with 40 points
dist = r./(c/fc);                           % Distance between subsequent position of Eve

alpha_step = 2;                           % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

% Communication parameters
Q = 64;
U = 4;
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



%% Matrices preallocation
%Bob matrices of size: nb_run x length(alpha) x length(U)
%Eve matrices of size: nb_run x length(alpha) x length(U) x nb_pos


e_an_TX             = zeros(nb_run,length(alpha),length(U));

e_sym_decod1_b      = zeros(nb_run,length(alpha),length(U));
e_noise_decod1_b    = zeros(nb_run,length(alpha),length(U));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),nb_pos);
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),nb_pos);
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),nb_pos);
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),nb_pos);

e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),nb_pos);
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),nb_pos);
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),nb_pos);
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),nb_pos);

e_sym_decod3_e      = zeros(nb_run,length(alpha),length(U),nb_pos);
e_noise_decod3_e    = zeros(nb_run,length(alpha),length(U),nb_pos);
e_an_decod3_e       = zeros(nb_run,length(alpha),length(U),nb_pos);
e_denom_decod3_e    = zeros(nb_run,length(alpha),length(U),nb_pos);

e_sym_decod4_e      = zeros(nb_run,length(alpha),length(U),nb_pos);
e_noise_decod4_e    = zeros(nb_run,length(alpha),length(U),nb_pos);
e_an_decod4_e       = zeros(nb_run,length(alpha),length(U),nb_pos);
e_denom_decod4_e    = zeros(nb_run,length(alpha),length(U),nb_pos);

e_sym_decod5_e      = zeros(nb_run,length(alpha),length(U),nb_pos);
e_noise_decod5_e    = zeros(nb_run,length(alpha),length(U),nb_pos);
e_an_decod5_e       = zeros(nb_run,length(alpha),length(U),nb_pos);
e_denom_decod5_e    = zeros(nb_run,length(alpha),length(U),nb_pos);



%% Random Variable Calculation
[H , RHO_SPATIAL] = corr_spatial( Q , b_subcar , r , fc , nb_run ) ;


for iter = 1:nb_run

Hb_TX = diag(squeeze(H(iter,:,1)).');
Hb_RX = ctranspose(Hb_TX);

for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                      % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');   % QAM modulation, can be changed to different modulation types, Nxnb_run



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


for dd = 1:nb_pos
%channel generation

He_TX = diag(squeeze(H(iter,:,dd)).');
He_RX = ctranspose(He_TX);


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
e_sym_decod1_b(iter,aa,bb,dd) = energy(sym_decod1_b);
e_noise_decod1_b(iter,aa,bb,dd) = energy(noise_decod1_b);

% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,dd)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,dd)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,dd)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,dd)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,aa,bb,dd)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb,dd)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

% @ Eve : decod 3
e_sym_decod3_e(iter,aa,bb,dd)     = energy(sym_decod3_e);
e_noise_decod3_e(iter,aa,bb,dd)   = energy(noise_decod3_e);
e_an_decod3_e(iter,aa,bb,dd)      = energy(an_decod3_e);
e_denom_decod3_e(iter,aa,bb,dd)   = energy(noise_decod3_e + an_decod3_e);        % energy of the sinr denominator for decoder 3 @Eve

% @ Eve : decod 4
e_sym_decod4_e(iter,aa,bb,dd)     = energy(sym_decod4_e);
e_noise_decod4_e(iter,aa,bb,dd)   = energy(noise_decod4_e);
e_an_decod4_e(iter,aa,bb,dd)      = energy(an_decod4_e);
e_denom_decod4_e(iter,aa,bb,dd)   = energy(noise_decod4_e + an_decod4_e);        % energy of the sinr denominator for decoder 4 @Eve

% @ Eve : decod 5
e_sym_decod5_e(iter,aa,bb,dd)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,aa,bb,dd)   = energy(noise_decod5_e);
e_an_decod5_e(iter,aa,bb,dd)      = energy(an_decod5_e);
e_denom_decod5_e(iter,aa,bb,dd)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve




end
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


figure;
plot(dist,sr1_avg(sr1_avg == max(sr1_avg))); hold on;
plot(dist,sr2_avg(sr2_avg == max(sr2_avg))); hold on;
plot(dist,sr5_avg(sr5_avg == max(sr5_avg))); hold on;
plot(dist,sr3_avg(sr3_avg == max(sr3_avg))); hold on;
plot(dist,sr4_avg(sr4_avg == max(sr4_avg))); hold on;





