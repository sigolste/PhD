%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   Implementation of new approx for Bob capa and Eve capa with SDS/MF/OC 
%   decoders when  frequency correlation among Bob's subcarriers is  
%   introduced. 
%   Computation of the required SNR @B to ensure a communication SR = X 
%   To ensure: i.e. No AWGN @E
%
%
%   Code Started:   09.02.2021
%   Last Update:    11.02.2021
%
%  
%   Mathematical derivations: cfr FC35-36
%
%  © SIDNEY GOLSTEIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% PROBLEME APPARENT: MODELE NE CHANGE PAS TELLEMENT AVEC LA CORRELATION!
% TESTER CHAQUE TERME SEPAREMENT
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
nb_run = 1000;                               % number of experiments
fc = 2e9 ;                                  % Carrier frequency
c = 3e8;                                    % Light speed

alpha_step = 1;                           % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

% Communication parameters
Q = 32;
U = [2 4 8 16];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters

snr_b  = [5 10].';  % SNR @Bob
snr_e  = [200];  % SNR @Eve

% Noise energy
sigma_b = 1./U./10.^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10.^(snr_e./10);    % expected noise energy @Eve


%% Channel modeling

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance
sigma_tau = .5e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                        % Approximation of coherence bandwidth
coef_freq = [1:15].*N.'/6;                               % From min/6 ∆fc to max/6 ∆fc, depending on N
delta_f_n = coef_freq.*delta_f_c;   
b_subcar = delta_f_n./N.';                                    % Bandwidth of each subcarrier
x_axis  = b_subcar./delta_f_c;                              % number of coherence bandwidth where lies a subcarrier 
x_axis = x_axis(1,:);

%% Matrix instantiation
e_noise_e            = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_b            = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_an_TX              = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

e_sym_decod1_b       = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_decod1_b     = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));


e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

e_sym_decod3_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_decod3_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_an_decod3_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_denom_decod3_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

e_sym_decod4_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_decod4_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_an_decod4_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_denom_decod4_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

e_sym_decod5_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_noise_decod5_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_an_decod5_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
e_denom_decod5_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

sym_b_exp4          = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
sym_decod1_e_exp4   = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
sym_decod2_e_exp4   = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
sym_decod5_e_exp4   = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

noise_b_exp4        = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
denom_decod1_e_exp4 = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
denom_decod2_e_exp4 = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));
denom_decod5_e_exp4 = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b));

H       = zeros(nb_run,Q,length(U),size(b_subcar,2));
H1      = zeros(nb_run,Q,length(U),size(b_subcar,2));
T       = zeros(Q,Q,length(U),size(b_subcar,2));
abs_rho = zeros(Q,length(U),size(b_subcar,2));

% Channel Generation:
for bb = 1:length(U)
    for dd = 1:size(b_subcar,2)
        [H(:,:,bb,dd), H1(:,:,bb,dd), abs_rho(:,bb,dd), T(:,:,bb,dd)] = corr_frequency( Q , b_subcar(bb,dd) , sigma_tau , nb_run ) ;
    end
end


for iter = 1:nb_run

for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


for dd = 1:size(b_subcar,2)
%channel generation
Hb_TX = diag(squeeze(H1(iter,:,bb,dd)).');
Hw_TX = diag(squeeze(H(iter,:,bb,dd)).');

Hb_RX = ctranspose(Hb_TX);
Hw_RX = ctranspose(Hw_TX);

He_TX = channelRayleigh(Q, mu , sigma);             % Assumption: uncorrelated subcarrier for Eve channel
He_RX = ctranspose(He_TX);


%% Encoder

% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % Qx1


sym_precoded = Hb_TX*sym_spread;


% AN generation
%[an,V1_correl(:,:,iter),S_correl(:,:,iter)] = generateAN_TEST(Hb_RX,Q(bb),U(bb),matrix_despread,1/U(bb),"svd");% % Qx1, not weighted 
[an, ~ , ~ ] = generateAN_TEST(Hb_RX,Q,U(bb),matrix_despread,1/U(bb),"svd");% % Qx1, not weighted 


for aa = 1:length(alpha)
% Weighting
sym_precoded_TX = sqrt(alpha(aa))*sym_precoded;     % weighted

an_TX = sqrt(1-alpha(aa))*an;                % weighted
e_an_TX(iter,aa) = energy(an_TX);

%% Receiver
% Whole stream
sym_RX_b = Hb_RX*(sym_precoded_TX + an_TX);
sym_RX_e = He_RX*(sym_precoded_TX + an_TX);

% Useful symbol
sym_b = Hb_RX*sym_precoded_TX; % Qx1
sym_e = He_RX*sym_precoded_TX;



for nn = 1:length(snr_b)
% Noise symbol
[noise_b, e_noise_b(iter,aa,bb,dd,nn) ] = addNoise(sym_b , snr_b(nn), 1/U(bb));

[noise_e, e_noise_e(iter,aa,bb,dd,nn) ] = addNoise(sym_e , snr_e, 1/U(bb));

% AN symbol
an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob


%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_RX*He_TX;                 % matched filter
decod3 = matrix_despread*(Hb_RX/He_RX);                     % AN killer
% gamma_E = (He_RX*Hb_TX)*matrix_spread;
% gamma_EH = ctranspose(gamma_E);
% decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an));   % LMMSE
decod5 = matrix_despread*He_TX;                             % Only He known by Eve
% 
sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_decod2_e = decod2*sym_e;
sym_decod3_e = decod3*sym_e;
% sym_decod4_e = decod4*sym_e;
sym_decod5_e = decod5*sym_e;

sym_b_exp4(iter,aa,bb,dd,nn) = mean(abs(sym_decod1_b).^4); 
sym_decod1_e_exp4(iter,aa,bb,dd,nn) = mean(abs(sym_decod1_e).^4); 
sym_decod2_e_exp4(iter,aa,bb,dd,nn) = mean(abs(sym_decod2_e).^4); 
sym_decod5_e_exp4(iter,aa,bb,dd,nn) = mean(abs(sym_decod5_e).^4); 


% 
noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;
noise_b_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod1_b).^4); 

noise_decod2_e = decod2*noise_e;
noise_decod3_e = decod3*noise_e;
% noise_decod4_e = decod4*noise_e;
noise_decod5_e = decod5*noise_e;
% 
an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;
an_decod3_e = decod3*an_e;
% an_decod4_e = decod4*an_e;
an_decod5_e = decod5*an_e;

denom_decod1_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod1_e + an_decod1_e).^4);
denom_decod2_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod2_e + an_decod2_e).^4);
% denom_decod3_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod3_e + an_decod3_e).^4);
% denom_decod4_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod4_e + an_decod4_e).^4);
denom_decod5_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod5_e + an_decod5_e).^4);



%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,aa,bb,dd,nn)      = energy(sym_decod1_b); %  energy(matrix_despread*abs(Hb_RX).^2*matrix_spread);
e_noise_decod1_b(iter,aa,bb,dd,nn)    = energy(noise_decod1_b);



% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,dd,nn)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,dd,nn)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,dd,nn)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,dd,nn)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,aa,bb,dd,nn)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb,dd,nn)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb,dd,nn)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb,dd,nn)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

% @ Eve : decod 3
e_sym_decod3_e(iter,aa,bb,dd,nn)     = energy(sym_decod3_e);
e_noise_decod3_e(iter,aa,bb,dd,nn)   = energy(noise_decod3_e);
e_an_decod3_e(iter,aa,bb,dd,nn)      = energy(an_decod3_e);
e_denom_decod3_e(iter,aa,bb,dd,nn)   = energy(noise_decod3_e + an_decod3_e);        % energy of the sinr denominator for decoder 3 @Eve

% % @ Eve : decod 4
% e_sym_decod4_e(iter,aa,bb,dd,nn)     = energy(sym_decod4_e);
% e_noise_decod4_e(iter,aa,bb,dd,nn)   = energy(noise_decod4_e);
% e_an_decod4_e(iter,aa,bb,dd,nn)      = energy(an_decod4_e);
% e_denom_decod4_e(iter,aa,bb,dd,nn)   = energy(noise_decod4_e + an_decod4_e);        % energy of the sinr denominator for decoder 4 @Eve

% @ Eve : decod 5
e_sym_decod5_e(iter,aa,bb,dd,nn)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,aa,bb,dd,nn)   = energy(noise_decod5_e);
e_an_decod5_e(iter,aa,bb,dd,nn)      = energy(an_decod5_e);
e_denom_decod5_e(iter,aa,bb,dd,nn)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve


end
end
end
waitbar(iter / nb_run)
end
end

close(h)



% Ergodic energies
e_avg_sym_decod1_b      = squeeze(mean(e_sym_decod1_b,1)); 
e_avg_noise_decod1_b     = squeeze(mean(e_noise_decod1_b,1));


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

% e_avg_sym_decod4_e      = squeeze(mean(e_sym_decod4_e,1));
% e_avg_noise_decod4_e    = squeeze(mean(e_noise_decod4_e,1));
% e_avg_an_decod4_e       = squeeze(mean(e_an_decod4_e,1));
% e_avg_denom_decod4_e    = squeeze(mean(e_denom_decod4_e,1));

e_avg_sym_decod5_e      = squeeze(mean(e_sym_decod5_e,1));
e_avg_noise_decod5_e    = squeeze(mean(e_noise_decod5_e,1));
e_avg_an_decod5_e       = squeeze(mean(e_an_decod5_e,1));
e_avg_denom_decod5_e    = squeeze(mean(e_denom_decod5_e,1));

% Terms exponent 4 (needed for approximation with variance)
e_avg_sym_b_exp4 = squeeze(mean(sym_b_exp4,1));
e_avg_sym_decod1_e_exp4 = squeeze(mean(sym_decod1_e_exp4,1));
e_avg_sym_decod2_e_exp4 = squeeze(mean(sym_decod2_e_exp4,1));
e_avg_sym_decod5_e_exp4 = squeeze(mean(sym_decod5_e_exp4,1));

e_avg_noise_b_exp4 = squeeze(mean(noise_b_exp4,1));
e_avg_denom_decod1_e_exp4 = squeeze(mean(denom_decod1_e_exp4,1)); 
e_avg_denom_decod2_e_exp4 = squeeze(mean(denom_decod2_e_exp4,1)); 
e_avg_denom_decod5_e_exp4 = squeeze(mean(denom_decod5_e_exp4,1)); 


% Instantaneous SINR
sinr1_b = e_sym_decod1_b./e_noise_decod1_b;

sinr1_e = e_sym_decod1_e./e_denom_decod1_e;
sinr2_e = e_sym_decod2_e./e_denom_decod2_e;
sinr3_e = e_sym_decod3_e./e_denom_decod3_e;
% sinr4_e = e_sym_decod4_e./e_denom_decod4_e;
sinr5_e = e_sym_decod5_e./e_denom_decod5_e;



% Instantaneous capacity
capa1_b = log2(1+sinr1_b);

capa1_e = log2(1+sinr1_e);
capa2_e = log2(1+sinr2_e);
capa3_e = log2(1+sinr3_e);
% capa4_e = log2(1+sinr4_e);
capa5_e = log2(1+sinr5_e);

% Instantaneous SR
sr1 = secrecyCapacity(sinr1_b,sinr1_e);
sr2 = secrecyCapacity(sinr1_b,sinr2_e);
sr3 = secrecyCapacity(sinr1_b,sinr3_e);
% sr4 = secrecyCapacity(sinr1_b,sinr4_e);
sr5 = secrecyCapacity(sinr1_b,sinr5_e);


% Ergodic SINR
sinr1_b_avg      = squeeze(mean(sinr1_b,1)); 

sinr1_e_avg = squeeze(mean(sinr1_e,1)); 
sinr2_e_avg = squeeze(mean(sinr2_e,1)); 
sinr3_e_avg = squeeze(mean(sinr3_e,1)); 
% sinr4_e_avg = squeeze(mean(sinr4_e,1)); 
sinr5_e_avg = squeeze(mean(sinr5_e,1)); 


% Ergodic capacity
capa1_b_avg = squeeze(mean(capa1_b,1)); 

capa1_e_avg = squeeze(mean(capa1_e,1));
capa2_e_avg = squeeze(mean(capa2_e,1));
capa3_e_avg = squeeze(mean(capa3_e,1));
% capa4_e_avg = squeeze(mean(capa4_e,1));
capa5_e_avg = squeeze(mean(capa5_e,1));

% Ergodic Secrecy Rate
sr1_avg = squeeze(mean(sr1));%capa1_b_correl_avg - capa1_e_avg; 
sr2_avg = squeeze(mean(sr2));%capa1_b_correl_avg - capa2_e_avg;
sr3_avg = squeeze(mean(sr3));%capa1_b_correl_avg - capa3_e_avg;
% sr4_avg = squeeze(mean(sr4));%capa1_b_correl_avg - capa4_e_avg;
sr5_avg = squeeze(mean(sr5));%capa1_b_correl_avg - capa5_e_avg;

%% TMP SECTION : Numerical approximations for SDS/MF decoding structure @E

% Decoder 1

% sinr1_e_simu = e_avg_sym_decod1_e./e_avg_denom_decod1_e;
% sinr1_e_simu_square = sinr1_e_simu.^2;
% expected1_e_simu_square = e_avg_sym_decod1_e_exp4./e_avg_denom_decod1_e_exp4;
% var1_e_simu = expected1_e_simu_square - sinr1_e_simu_square;
% 
% capa1_e_ergodic = capa1_e_avg;
% capa1_e_jensen_simu = log2(1+sinr1_e_simu);
% capa1_e_new_approx_simu = log2(1+sinr1_e_simu) - 1/2*var1_e_simu./((1+sinr1_e_simu).^2);



% % Decoder 2: Matched Filter
% sinr2_e_simu = e_avg_sym_decod2_e./e_avg_denom_decod2_e;
% sinr2_e_simu_square = sinr2_e_simu.^2;
% expected2_e_simu_square = e_avg_sym_decod2_e_exp4./e_avg_denom_decod2_e_exp4;
% var2_e_simu = expected2_e_simu_square - sinr2_e_simu_square;
% 
% capa2_e_ergodic = capa2_e_avg;
% capa2_e_jensen_simu = log2(1+sinr2_e_simu);
% capa2_e_new_approx_simu = log2(1+sinr2_e_simu) - 1/2*var2_e_simu./((1+sinr2_e_simu).^2);
% 
% % Decoder 5 : Own channel knwoledge
% sinr5_e_simu = e_avg_sym_decod5_e./e_avg_denom_decod5_e;
% sinr5_e_simu_square = sinr5_e_simu.^2;
% expected5_e_simu_square = e_avg_sym_decod5_e_exp4./e_avg_denom_decod5_e_exp4;
% var5_e_simu = expected5_e_simu_square - sinr5_e_simu_square;
% 
% capa5_e_ergodic = capa5_e_avg;
% capa5_e_jensen_simu = log2(1+sinr5_e_simu);
% capa5_e_new_approx_simu = log2(1+sinr5_e_simu) - 1/2*var5_e_simu./((1+sinr5_e_simu).^2);


% h = waitbar(0,'Modeling SINR/VARIANCE/CAPA at Bob and Eve ...');
% tic;
% for bb = 1:length(U)
%     for dd = 1:size(b_subcar,2)
%         parfor nn = 1:length(snr_b)
%             sinr_b_model(:,bb,dd,nn)    = sinrModelingFrequencyCorrelation(alpha,U(bb),N(bb),T(:,:,bb,dd),snr_b(nn),snr_e,"bob_correl").';
%             var_b_model(:,bb,dd,nn)     = modelVariance(alpha,N(bb),U(bb),T(:,:,bb,dd),sinr_b_model(:,bb,dd,nn),snr_b(nn),snr_e,"bob_correl");
%             capa_b_jensen(:,bb,dd,nn)   = log2(1+sinr_b_model(:,bb,dd,nn));
%             capa_b_new_approx(:,bb,dd,nn)  = log2(1+sinr_b_model(:,bb,dd,nn)) - 1/2*var_b_model(:,bb,dd,nn)./((1+sinr_b_model(:,bb,dd,nn)).^2);
%             
%             sinr1_e_model(:,bb,dd,nn)    = sinrModelingFrequencyCorrelation(alpha,U(bb),N(bb),T(:,:,bb,dd),snr_b(nn),snr_e,"eve_decod1_correl").';
%             var1_e_model(:,bb,dd,nn)     = modelVariance(alpha,N(bb),U(bb),T(:,:,bb,dd),sinr1_e_model(:,bb,dd,nn),snr_b(nn),snr_e,"eve_decod1_correl");
%             capa1_e_jensen(:,bb,dd,nn)       = log2(1+sinr1_e_model(:,bb,dd,nn));
%             capa1_e_new_approx(:,bb,dd,nn)   = log2(1+sinr1_e_model(:,bb,dd,nn)) - 1/2*var1_e_model(:,bb,dd,nn)./((1+sinr1_e_model(:,bb,dd,nn)).^2);
%         end
%     end
%     waitbar(bb / length(U))
% end
% toc;
% close(h)



capa_b_ergodic      = capa1_b_avg;
% capa_b_jensen       = squeeze(capa_b_jensen);
% capa_b_new_approx   = squeeze(capa_b_new_approx);


capa1_e_ergodic      = capa1_e_avg;
% capa1_e_jensen       = squeeze(capa1_e_jensen);
% capa1_e_new_approx   = squeeze(capa1_e_new_approx);

% % ANALYTIC ERRORS
% error_1st_order_b = capa_b_ergodic - capa_b_jensen;
% error_2nd_order_b = capa_b_ergodic - capa_b_new_approx;
% 
% error_1st_order_1e = capa1_e_ergodic - capa1_e_jensen;
% error_2nd_order_1e = capa1_e_ergodic - capa1_e_new_approx;


% SIMULATION ERRORS
% error_1st_order_1e_simu = capa1_e_ergodic - capa1_e_jensen_simu;          % IDEM QUE MODEL
% error_2nd_order_1e_simu = capa1_e_ergodic - capa1_e_new_approx_simu;
% 
% error_1st_order_2e_simu = capa2_e_ergodic - capa2_e_jensen_simu;
% error_2nd_order_2e_simu = capa2_e_ergodic - capa2_e_new_approx_simu;
% 
% error_1st_order_5e_simu = capa5_e_ergodic - capa5_e_jensen_simu;
% error_2nd_order_5e_simu = capa5_e_ergodic - capa5_e_new_approx_simu;

% SR

% sr1_jensen = capa_b_jensen - capa1_e_jensen;
% sr1_new_approx = capa_b_new_approx - capa1_e_new_approx;
 sr1_ergodic = capa_b_ergodic - capa1_e_ergodic;


% sr2_jensen = capa_b_jensen - capa2_e_jensen_simu;
% sr2_new_approx = capa_b_new_approx - capa2_e_new_approx_simu;
sr2_ergodic = sr2_avg;

sr3_ergodic = sr3_avg;

% sr4_ergodic = sr4_avg;

% sr5_jensen = capa_b_jensen - capa5_e_jensen_simu;
% sr5_new_approx = capa_b_new_approx - capa5_e_new_approx_simu;
sr5_ergodic = sr5_avg;

%%
% %% Optimal alpha
% 
% % SDS Decoder:
% % Variable initialization
% T1B                  = zeros(length(U),size(b_subcar,2),length(snr_b));
% alpha_opt_1st_order1 = zeros(length(U),size(b_subcar,2),length(snr_b));
% 
% for bb = 1:length(U)
%     for dd = 1:size(b_subcar,2)
%         for nn = 1:length(snr_b)
%             H2H2 = modelCorrelH2H2(N(bb),U(bb),T(:,:,bb,dd));
%             H4 = modelCorrelH4(U(bb));
%             T1B(bb,dd,nn) = H2H2 + H4;
%             alpha_opt_1st_order1(bb,dd,nn) = (T1B(bb,dd,nn)*(U(bb)*sigma_e(bb)+1)-U(bb)^2*sigma_b(nn))/(2*T1B(bb,dd,nn)) ;
%         end
%     end
% end

% MF Decoder

%% Required Bob SNR in order to ensure SR = X bit/channel use

% Hypothesis:
% Noiseless Eve receiver
% No correlation among Eve's subcarriers

% Parameters:
SR_targeted = linspace(0.05,5,50);%[0.1:.2:10] ;  % Targetted SR from 0.5 bit/channel use to 5 bits/channel use
sr = 2.^SR_targeted ;       % Linear scale


snr_tmp_1b = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b),length(sr));   % With time dependency --> 6D matrix
snr_tmp_2b = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b),length(sr));
snr_tmp_3b = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b),length(sr));
% snr_tmp_4b = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b),length(sr));
snr_tmp_5b = zeros(nb_run,length(alpha),length(U),size(b_subcar,2),length(snr_b),length(sr));


for ss = 1:length(sr)
    for bb = 1:length(U)
    snr_tmp_1b(:,:,bb,:,:,ss) = 10*log10( ( sr(ss)*(e_an_decod1_e(:,:,bb,:,:) + e_sym_decod1_e(:,:,bb,:,:)) - e_an_decod1_e(:,:,bb,:,:) )  ... 
                            ./ (U(bb)*e_sym_decod1_b(:,:,bb,:,:).*e_an_decod1_e(:,:,bb,:,:)) );
    snr_tmp_2b(:,:,bb,:,:,ss) = 10*log10( ( sr(ss)*(e_an_decod2_e(:,:,bb,:,:) + e_sym_decod2_e(:,:,bb,:,:)) - e_an_decod2_e(:,:,bb,:,:) )  ... 
                            ./ (U(bb)*e_sym_decod1_b(:,:,bb,:,:).*e_an_decod2_e(:,:,bb,:,:)) );
    snr_tmp_3b(:,:,bb,:,:,ss) = 10*log10( ( sr(ss)*(e_an_decod3_e(:,:,bb,:,:) + e_sym_decod3_e(:,:,bb,:,:)) - e_an_decod3_e(:,:,bb,:,:) )  ... 
                            ./ (U(bb)*e_sym_decod1_b(:,:,bb,:,:).*e_an_decod3_e(:,:,bb,:,:)) );
%     snr_tmp_4b(:,:,bb,:,:,ss) = 10*log10( ( sr(ss)*(e_an_decod4_e(:,:,bb,:,:) + e_sym_decod4_e(:,:,bb,:,:)) - e_an_decod4_e(:,:,bb,:,:) )  ... 
%                             ./ (U(bb)*e_sym_decod1_b(:,:,bb,:,:).*e_an_decod4_e(:,:,bb,:,:)) );
    snr_tmp_5b(:,:,bb,:,:,ss) = 10*log10( ( sr(ss)*(e_an_decod5_e(:,:,bb,:,:) + e_sym_decod5_e(:,:,bb,:,:)) - e_an_decod5_e(:,:,bb,:,:) )  ... 
                            ./ (U(bb)*e_sym_decod1_b(:,:,bb,:,:).*e_an_decod5_e(:,:,bb,:,:)) );
    end
end

% snr_tmp_ergo_1b = squeeze(mean(squeeze(mean(snr_tmp_1b,1)),4));             
snr_tmp_ergo_1b = squeeze(mean(squeeze(mean(snr_tmp_1b,1)),4)); % Mean over the number of run and over the dimension on SNRb since it should not depend on Bob SNR . --> 4D matrix
snr_tmp_ergo_2b = squeeze(mean(squeeze(mean(snr_tmp_2b,1)),4));
snr_tmp_ergo_3b = squeeze(mean(squeeze(mean(snr_tmp_3b,1)),4));
% snr_tmp_ergo_4b = squeeze(mean(squeeze(mean(snr_tmp_4b,1)),4));
snr_tmp_ergo_5b = squeeze(mean(squeeze(mean(snr_tmp_5b,1)),4));

% snr_tmp_ergo_1/2/5b is a 4D matrix with the required SNR @B, depending on
% the decoding structure investigated, that allows to achieve a targetted
% SR = x bit/channel use, as a function of alpha, U, amount of correlation,
% and the targetted SR. 
% For a particular set of (U,correlation,SR), we obtain a SNR value for
% each alpha value. We now have to determine which value of alpha gives the
% lowest (i.e. the easiest to achieve) value of required SNR @B. In doing
% so, we will have the lowest SNR required and the optimal amount of alpha 
% to inject in order to target a SR = X for a particular triple set
% (U,correlation,SR)

% --> 3D matrix
[snr_required_1b , idx_1b] = min(snr_tmp_ergo_1b,[],1); % find the minimum required snr along the first dimension and the idx in the first dimension
[snr_required_2b , idx_2b] = min(snr_tmp_ergo_2b,[],1); % find the minimum required snr along the first dimension and the idx in the first dimension
[snr_required_3b , idx_3b] = min(snr_tmp_ergo_3b,[],1); % find the minimum required snr along the first dimension and the idx in the first dimension
% [snr_required_4b , idx_4b] = min(snr_tmp_ergo_4b,[],1); % find the minimum required snr along the first dimension and the idx in the first dimension
[snr_required_5b , idx_5b] = min(snr_tmp_ergo_5b,[],1); % find the minimum required snr along the first dimension and the idx in the first dimension

snr_required_1b     = squeeze(snr_required_1b); % value of required SNR @B for each triple (U,correlation,SR) 
idx_1b              = squeeze(idx_1b);
alpha_required_1b   = alpha(idx_1b);            % Data energy to inject for each triple (U,correlation,SR)         

snr_required_2b     = squeeze(snr_required_2b); % value of required SNR @B for each triple (U,correlation,SR) 
idx_2b              = squeeze(idx_2b);
alpha_required_2b   = alpha(idx_2b);            % Data energy to inject for each triple (U,correlation,SR) 

snr_required_3b     = squeeze(snr_required_3b); % value of required SNR @B for each triple (U,correlation,SR) 
idx_3b              = squeeze(idx_3b);
alpha_required_3b   = alpha(idx_3b);            % Data energy to inject for each triple (U,correlation,SR)   

% snr_required_4b     = squeeze(snr_required_4b); % value of required SNR @B for each triple (U,correlation,SR) 
% idx_4b              = squeeze(idx_4b);
% alpha_required_4b   = alpha(idx_4b);            % Data energy to inject for each triple (U,correlation,SR)         
% 
snr_required_5b     = squeeze(snr_required_5b); % value of required SNR @B for each triple (U,correlation,SR) 
idx_5b              = squeeze(idx_5b);
alpha_required_5b   = alpha(idx_5b);            % Data energy to inject for each triple (U,correlation,SR)         


%% PLOT SECTION 
% Default parameters
set(0,'defaultAxesFontSize',22)
set(0,'DefaultLineMarkerSize',12);
set(0,'defaultLineLineWidth',2);

% SNR/alpha en fonction SR targetted
figure;
subplot(2,1,1)
plot(SR_targeted,squeeze(snr_required_1b(1,1,:)).','Marker','o'); hold on;
plot(SR_targeted,squeeze(snr_required_2b(1,1,:)).','Marker','v'); hold on;
plot(SR_targeted,squeeze(snr_required_5b(1,1,:)).','Marker','diamond'); hold on;
title(['BOR = 2, $\Delta f_C / \Delta f_N$ = 1/6' ])
box on; grid on; 
xlabel('Targeted SR (bit/channel use)')
ylabel('Required SNR at Bob (dB)')
legend('SDS Decoder','MF Decoder', 'OC Decoder','location','bestoutside')


subplot(2,1,2)
plot(SR_targeted,100*(1-squeeze(alpha_required_1b(1,1,:)).'),'Marker','o'); hold on;
plot(SR_targeted,100*(1-squeeze(alpha_required_2b(1,1,:)).'),'Marker','v'); hold on;
plot(SR_targeted,100*(1-squeeze(alpha_required_5b(1,1,:)).'),'Marker','diamond'); hold on;
title(['BOR = 2, $\Delta f_C / \Delta f_N$ = 1/6' ])
box on; grid on; 
xlabel('Targeted SR (bit/channel use)')
ylabel('Required AN energy to inject ($\%$)')
legend('SDS Decoder','MF Decoder', 'OC Decoder','location','bestoutside')

% SNR/alpha en fonction du BOR
figure;
subplot(2,1,1)
plot(U,squeeze(snr_required_1b(:,1,end)).','Marker','o'); hold on;
plot(U,squeeze(snr_required_2b(:,1,end)).','Marker','v'); hold on;
plot(U,squeeze(snr_required_5b(:,1,end)).','Marker','diamond'); hold on;
title(['Targeted SR = ' num2str(SR_targeted(end)), ' bit/channel use, $\Delta f_C / \Delta f_N$ = 1/6' ])
box on; grid on; 
xlabel('BOR')
ylabel('Required SNR at Bob (dB)')
legend('SDS Decoder','MF Decoder', 'OC Decoder','location','bestoutside')


subplot(2,1,2)
plot(U,100*(1-squeeze(alpha_required_1b(:,1,end)).'),'Marker','o'); hold on;
plot(U,100*(1-squeeze(alpha_required_2b(:,1,end)).'),'Marker','v'); hold on;
plot(U,100*(1-squeeze(alpha_required_5b(:,1,end)).'),'Marker','diamond'); hold on;
title(['Targeted SR = ' num2str(SR_targeted(end)), ' bit/channel use, $\Delta f_C / \Delta f_N$ = 1/6' ])
box on; grid on; 
xlabel('BOR')
ylabel('Required AN energy to inject ($\%$)')
legend('SDS Decoder','MF Decoder', 'OC Decoder','location','bestoutside')



% SNR/alpha en fonction de la correlation
figure;
subplot(2,1,1)
plot(x_axis,squeeze(snr_required_1b(1,:,end)).','Marker','o'); hold on;
plot(x_axis,squeeze(snr_required_2b(1,:,end)).','Marker','v'); hold on;
plot(x_axis,squeeze(snr_required_5b(1,:,end)).','Marker','diamond'); hold on;
title(['Targeted SR = ' num2str(SR_targeted(end)), ' bit/channel use, $\Delta f_C / \Delta f_N$ = 1/6' ])
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('Required SNR at Bob (dB)')
legend('SDS Decoder','MF Decoder', 'OC Decoder','location','bestoutside')


subplot(2,1,2)
plot(x_axis,100*(1-squeeze(alpha_required_1b(1,:,end)).'),'Marker','o'); hold on;
plot(x_axis,100*(1-squeeze(alpha_required_2b(1,:,end)).'),'Marker','v'); hold on;
plot(x_axis,100*(1-squeeze(alpha_required_5b(1,:,end)).'),'Marker','diamond'); hold on;
title(['Targeted SR = ' num2str(SR_targeted(end)), ' bit/channel use, BOR = 2' ])
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('Required AN energy to inject ($\%$)')
legend('SDS Decoder','MF Decoder', 'OC Decoder','location','bestoutside')






% nn = find(snr_b == 10);

