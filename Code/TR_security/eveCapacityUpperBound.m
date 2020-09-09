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
Q = 4;
U = [4];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 20; % energy per bit over noise psd @Bob - dB
EbN0_e = EbN0_b; % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);  % SNR @Eve

sigma_e = 1./10^(snr_e./10);
sigma_b = 1./10^(snr_b./10);
% Channel parameters 
mu = 0;         % Mean channel
sigma = 1;      % Variance of channel




%% Mainloop
for iter = 1:nb_run
for bb =1:length(U)
% msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
% sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run
% 
% Channel generation
% Rayleigh channel - independant subcarrier fading - unit power per
% subcarrier : CN(0,I_Q) - circulary complex random normal variable.
Hb_TX = channelRayleigh(U(bb), mu , sigma);  
He_TX = channelRayleigh(U(bb), mu , sigma);
Hb_RX = ctranspose(Hb_TX.');
He_RX = ctranspose(He_TX.');


mu_p = Hb_TX/norm(Hb_TX);
%mu_p = Hb_TX;
Qw = eye(U(bb))-mu_p*ctranspose(mu_p);
% (De)-Spreading matrix creation
% [matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));

%% Encoder

% Spreading + TR precoding
% sym_spread = matrix_spread*sym_TX;  % Qx1
% sym_precoded = Hb_TX*sym_spread; % Qx1, not weighted 
% 
% % AN generation
% an = generateAN(Hb_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd"); % Qx1, not weighted
% sym_spread = matrix_spread*sym_TX;  % Qx1
% sym_precoded = Hb_TX*sym_spread; % Qx1, not weighted

for aa = 1:length(alpha)

% sym_precoded_TX = sqrt(alpha(aa))*sym_precoded;     % weighted
% an_TX = sqrt(1-alpha(aa))*an;                       % weighted
% 
% sym_b = Hb_RX*sym_precoded_TX; % Qx1
% sym_e = He_RX*sym_precoded_TX;
% 
% 
% [noise_e, e_noise_e ] = addNoise(sym_e , snr_e, energy(sym_precoded_TX+an_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   

a = sqrt(alpha(aa));
b = sqrt(1-alpha(aa));

term1(iter,aa,bb) = log2(det(a^2*He_RX*mu_p*ctranspose(mu_p)*He_TX + b^2*He_TX*Qw*He_RX + sigma_e*eye(U(bb))));
term2(iter,aa,bb) = log2(det(b^2*He_RX*Qw*He_TX + sigma_e*eye(U(bb))));

end
waitbar(iter / nb_run)
end
end

term1_mean = squeeze(mean(term1,1));
term2_mean = squeeze(mean(term2,1));



capa_max_eve = term1_mean-term2_mean;
capa_bob = log2(1 + alpha.*(U+1)./U./sigma_b); 

sr_min = capa_bob - capa_max_eve; 
figure;
plot(100*alpha,real(sr_min))
% upper_bound_per_symb = U*upper_bound_per_subcar;
close(h)


%% Post processing
% Energy averaging over all realisations




