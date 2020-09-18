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
alpha = (0:alpha_step/100:1).';         

% Communication parameters
U = [2 4 8 16];
% N = Q./U;

M = 4;
k = log2(M);

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
Hb_TX = channelRayleigh(U(bb), mu , sigma);  
He_TX = channelRayleigh(U(bb), mu , sigma);
Hb_RX = ctranspose(Hb_TX.');
He_RX = ctranspose(He_TX.');


mu_p = Hb_TX/norm(Hb_TX) ; 
%mu_p = Hb_TX;
Qw = eye(U(bb))-mu_p*ctranspose(mu_p);

%% Encoder

for aa = 1:length(alpha)

a = sqrt(alpha(aa));
b = sqrt(1-alpha(aa));

term1(iter,aa,bb) = log2(det(a^2*He_RX*mu_p*ctranspose(mu_p)*He_TX + b^2*He_TX*Qw*He_RX + sigma_e(bb)*eye(U(bb))));
term2(iter,aa,bb) = log2(det(b^2*He_RX*Qw*He_TX + sigma_e(bb)*eye(U(bb))));

end
waitbar(iter / nb_run)
end
end

term1_mean = squeeze(mean(term1,1)); % Mean over all realizations
term2_mean = squeeze(mean(term2,1));



capa_max_eve = term1_mean-term2_mean;               % Upper bound capa @ Eve
capa_bob = log2(1 + alpha.*(U+1)./U./sigma_b);      % Analytic capa @ Bob
sr_min = capa_bob - capa_max_eve;                   % minimum of secrecy rate

figure;
plot(100*alpha,real(sr_min)); hold on;
legend('BOR = 2', 'BOR = 4', 'BOR = 8', 'BOR = 16','Location', 'Best' )
box on ; grid on;
xlabel('Percentage of energy dedicated for data (\%)')
ylabel('Secrecy rate (bit/channel use)')
ylim([-20 5])

% upper_bound_per_symb = U*upper_bound_per_subcar;
close(h)






