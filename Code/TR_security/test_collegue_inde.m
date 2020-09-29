%%writefile sidney.m
sigma=1;

clear all;
close all;

%%%%SIBI: All Decorations Removed, as they wont run in google Colab.

%h = waitbar(0,'Simulation Progression...');

% function y = channelRayleigh(len, mu, sigma)
%   y = (1/sqrt(2)*(sigma*randn(len,1) + sigma*1j*randn(len,1)) + mu);
% end

%% Parameters
% Simulation parameters
nb_run = 2;                 % number of experiments
alpha_step = 1;             % Percentage between subsequent alpha values
alpha = (0:alpha_step/20:1).';         

% Communication parameters
U = [4];
% N = Q./U;

M = 4;
k = log2(M); %%%% SIBI Why M=4, we are not dealing with just QPSK.

% AWGN parameters
EbN0_b = 20;                                % energy per bit over noise psd @Bob - dB
EbN0_e = EbN0_b;                            % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);              % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);              % SNR @Eve
      
sigma_e = 1./U./10^(snr_e./10);
sigma_b = 1./U./10^(snr_b./10);
% Channel parameters 
mu = 0;         % Mean channel
sigma = 1;      % Variance of channel


%% Mainloop
for iter = 1:nb_run
for bb =1:length(U)

% Channel generation
% Rayleigh channel - independant subcarrier fading - unit power per
% subcarrier : CN(0,I_Q) - circulary complex random normal variable.
Hb_TX = diag(channelRayleigh(U(bb), mu , sigma));  
He_TX = diag(channelRayleigh(U(bb), mu , sigma));
Hb_RX = ctranspose(Hb_TX);
He_RX = ctranspose(He_TX);


mu_p = Hb_TX/norm(Hb_TX);
%mu_p = Hb_TX;
Qw = eye(U(bb))-mu_p*ctranspose(mu_p);

%% Encoder

for aa = 1:length(alpha)

a = sqrt(alpha(aa));
b = sqrt(1-alpha(aa));

term1e(iter,aa,bb) = log2(a^2*He_RX*mu_p*ctranspose(mu_p)*He_TX + b^2*He_RX*Qw*He_TX + sigma_e(bb)) ;
term2e(iter,aa,bb) =  log2(b^2*He_RX*Qw*He_TX + sigma_e(bb));

term1b(iter,aa,bb) = log2(a^2*Hb_RX*mu_p*ctranspose(mu_p)*Hb_TX + b^2*Hb_RX*Qw*Hb_TX + sigma_b(bb)) ;
term2b(iter,aa,bb) =  log2(b^2*Hb_RX*Qw*Hb_TX + sigma_b(bb));

%%%%SIBI Notice that in above Hb_RX*Qw*Hb_TX = 0, due to the choice of mu_p.


end
waitbar(iter / nb_run)
end
end

term1e_mean = squeeze(mean(term1e,1)); % Mean over all realizations
term2e_mean = squeeze(mean(term2e,1));

%%%%SIBI: Changed the lines below.
term1b_mean = squeeze(mean(term1b,1)); % Mean over all realizations
term2b_mean = squeeze(mean(term2b,1));


capa_max_eve = term1e_mean-term2e_mean;               % Upper bound capa @ Eve
capa_max_bob = term1b_mean-term2b_mean;               % Upper bound capa @ Eve

capa_bob = log2(1 + alpha.*((U+1)./U)./sigma_b).';      % Analytic capa @ Bob
sr_min = capa_bob - capa_max_eve;                   % minimum of secrecy rate

sr_min