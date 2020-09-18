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

Na = 2;             % Number of antennas @Alice
Nb = 1;             % Number of antennas @Bob
Ne = 1;             % Number of antennas @Eve
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


%% Matrices preallocation
e_an_TX             = zeros(nb_run,length(alpha),length(U));

e_sym_decod1_b      = zeros(nb_run,length(alpha),length(U));
e_noise_decod1_b    = zeros(nb_run,length(alpha),length(U));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U));
