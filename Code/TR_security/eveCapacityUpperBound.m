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
alpha_step = 5;             % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

% Communication parameters
Q = 16;
U = [2];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 20; % energy per bit over noise psd @Bob - dB
EbN0_e = EbN0_b; % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Mean channel
sigma = 1;      % Variance of channel



%% Matrices preallocation
e_denom_decod4_e    = zeros(nb_run,length(alpha),length(U));


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
sym_spread = matrix_spread*sym_TX;  % Qx1
sym_precoded = Hb_TX*sym_spread; % Qx1, not weighted

for aa = 1:length(alpha)

sym_precoded_TX = sqrt(alpha(aa))*sym_precoded;     % weighted
an_TX = sqrt(1-alpha(aa))*an;                       % weighted

sym_b = Hb_RX*sym_precoded_TX; % Qx1
sym_e = He_RX*sym_precoded_TX;


[noise_e, e_noise_e ] = addNoise(sym_e , snr_e, energy(sym_precoded_TX+an_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   
%mu_p = Hb_TX/norm(Hb_TX);
mu_p = Hb_TX;

a = sqrt(alpha(aa));
b = sqrt(1-alpha(aa));

term1(:,iter,aa,bb) = diag(log2(abs((a^2-b^2)*He_RX*mu_p*ctranspose(mu_p)*He_TX + b^2*He_TX*He_RX + e_noise_e*eye(Q))));
term2(:,iter,aa,bb) = diag(log2(abs(b^2*(eye(Q) - mu_p*ctranspose(mu_p))*He_TX + e_noise_e*eye(Q))));

end
waitbar(iter / nb_run)
end
end

term1_mean = squeeze(mean(term1,2));
term2_mean = squeeze(mean(term2,2));

term1_mean_per_subcar = mean(term1_mean,1);
term2_mean_per_subcar = mean(term2_mean,1);


upper_bound_per_subcar = term1_mean_per_subcar-term2_mean_per_subcar;

upper_bound_per_symb = U*upper_bound_per_subcar;
close(h)

plot(upper_bound_per_symb)

%% Post processing
% Energy averaging over all realisations




