% Code Started 21.10.2020
% Last Update: 
% By Sidney Golstein


clear all;
%close all;

set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',32)
set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'DefaultLineMarkerSize',15);
set(0, 'defaultFigurePosition',  [-1267  44   1256    872])


N_A = 2;       % nb of TX antennas
N_B = 3;     % nb of RX antennas @Bob
N_E = 4;     % nb of RX antennas @Eve


% Communication parameters
Q = 512;
U = [4];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% ELM parameters
c = 3e8;      
fc = 2.45e9 ;                                                               % Carrier frequency
lambda_c = c / fc ;                                                         % Carrier wavelength (m)
sigma_tau = 3e-6 ;                                                          % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                                        % Approximation of coherence bandwidth


% Spatial parameters
coef_space_A = 0.001;                                                         % Spacement between antennas in lambda_c         
coef_space_B = 10 ;
coef_space_E = 0.1;


% frequency parameters
coef_freq_B = N/6;                                                            % Subcarrier spacement in fraction of delta_fc for AB chanel
coef_freq_E = N;



[Hb_RX , Hb_space_RX , Hb_freq_RX , Hb_no_corr_RX] = channelMIMO(fc , c , lambda_c , sigma_tau, Q , N , N_A , N_B , coef_space_A , coef_space_B , coef_freq_B);
[He_RX , He_space_RX , He_freq_RX , He_no_corr_RX] = channelMIMO(fc , c , lambda_c , sigma_tau, Q , N , N_A , N_E , coef_space_A , coef_space_E , coef_freq_E);



figure;
set(gcf,'name','2 colocated TX, 3 RX')
subplot(4,1,1)
for ii = 1:N_B
plot(squeeze(abs(Hb_no_corr_RX(ii,1,:)))); hold on;
plot(squeeze(abs(Hb_no_corr_RX(ii,2,:)))); hold on;
end
title('No correlation')

subplot(4,1,2)
for ii = 1:N_B
plot(squeeze(abs(Hb_space_RX(ii,1,:)))); hold on;
plot(squeeze(abs(Hb_space_RX(ii,2,:)))); hold on;
end
title('Spatial correlation of TX antennas')

subplot(4,1,3)
for ii = 1:N_B
plot(squeeze(abs(Hb_freq_RX(ii,1,:)))); hold on;
plot(squeeze(abs(Hb_freq_RX(ii,2,:)))); hold on;
end
title('Frequency correlation of TX antennas')


subplot(4,1,4)
for ii = 1:N_B
plot(squeeze(abs(Hb_RX(ii,1,:)))); hold on;
plot(squeeze(abs(Hb_RX(ii,2,:)))); hold on;
end
title('Spatial correlation of TX antennas + frequency correlation')



