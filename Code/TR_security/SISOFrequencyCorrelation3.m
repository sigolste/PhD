%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   Implementation of new approx for Bob capa and Eve capa with SDS/MF/OC 
%   decoders when  frequency correlation among Bob's subcarriers is  
%   introduced. Check the error matrices due to the approximations. These
%   matrices depend on all parameters (U, rho, snr_b, alpha,...). The idea
%   is to check if there is at least one parameter for which the
%   approximated SR (1st and/or 2nd order) converges to the ergodic one. 
%
%
%
%   Code Started:   09.02.2021
%   Last Update:    11.02.2021
%
%   © SIDNEY GOLSTEIN
% 
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
nb_run = 4000;                               % number of experiments
fc = 2e9 ;                                  % Carrier frequency
c = 3e8;                                    % Light speed

alpha_step = 10;                           % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

% Communication parameters
Q = 16;
U = [2 4 8 16];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters

snr_b  = [-5 0 5 10 15].';  % SNR @Bob
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
coef_freq = [1 3 6 9 12 15 18 21 24].*N.'/6;                  % From 1/6 ∆fc to 4 ∆fc, depending on N
delta_f_n = coef_freq.*delta_f_c;   
b_subcar = delta_f_n./N.';                                    % Bandwidth of each subcarrier
x_axis  = b_subcar./delta_f_c;                              % number of coherence bandwidth where lies a subcarrier 
for bb = 1:length(U)
    for dd = 1:size(b_subcar,2)
        [H(:,:,bb,dd), H1(:,:,bb,dd), abs_rho(:,bb,dd), T(:,:,bb,dd)] = corr_frequency( Q , b_subcar(bb,dd) , sigma_tau , nb_run ) ;
    end
end

%% Energy matrix instantiation
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
% decod3 = matrix_despread*(Hb_RX/He_RX);                     % AN killer
% gamma_E = (He_RX*Hb_TX)*matrix_spread;
% gamma_EH = ctranspose(gamma_E);
% decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an) + e_noise_e*eye(Q) );   % LMMSE
decod5 = matrix_despread*He_TX;                             % Only He known by Eve
% 
sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_decod2_e = decod2*sym_e;
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
% noise_decod3_e = decod3*noise_e;
% noise_decod4_e = decod4*noise_e;
noise_decod5_e = decod5*noise_e;
% 
an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;
an_decod5_e = decod5*an_e;

denom_decod1_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod1_e + an_decod1_e).^4);
denom_decod2_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod2_e + an_decod2_e).^4);
denom_decod5_e_exp4(iter,aa,bb,dd,nn) = mean(abs(noise_decod5_e + an_decod5_e).^4);


% an_decod3_e = decod3*an_e;
% an_decod4_e = decod4*an_e;

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

% % @ Eve : decod 3
% e_sym_decod3_e(iter,aa,bb,dd)     = energy(sym_decod3_e);
% e_noise_decod3_e(iter,aa,bb,dd)   = energy(noise_decod3_e);
% e_an_decod3_e(iter,aa,bb,dd)      = energy(an_decod3_e);
% e_denom_decod3_e(iter,aa,bb,dd)   = energy(noise_decod3_e + an_decod3_e);        % energy of the sinr denominator for decoder 3 @Eve
% 
% % @ Eve : decod 4
% e_sym_decod4_e(iter,aa,bb,dd)     = energy(sym_decod4_e);
% e_noise_decod4_e(iter,aa,bb,dd)   = energy(noise_decod4_e);
% e_an_decod4_e(iter,aa,bb,dd)      = energy(an_decod4_e);
% e_denom_decod4_e(iter,aa,bb,dd)   = energy(noise_decod4_e + an_decod4_e);        % energy of the sinr denominator for decoder 4 @Eve
% 
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
sinr5_e = e_sym_decod5_e./e_denom_decod5_e;



% Instantaneous capacity
capa1_b = log2(1+sinr1_b);

capa1_e = log2(1+sinr1_e);
capa2_e = log2(1+sinr2_e);
capa5_e = log2(1+sinr5_e);

% Instantaneous SR
sr1 = secrecyCapacity(sinr1_b,sinr1_e);
sr2 = secrecyCapacity(sinr1_b,sinr2_e);
sr5 = secrecyCapacity(sinr1_b,sinr5_e);


% Ergodic SINR
sinr1_b_avg      = squeeze(mean(sinr1_b,1)); 

sinr1_e_avg = squeeze(mean(sinr1_e,1)); 
sinr2_e_avg = squeeze(mean(sinr2_e,1)); 
sinr5_e_avg = squeeze(mean(sinr5_e,1)); 


% Ergodic capacity
capa1_b_avg = squeeze(mean(capa1_b,1)); 

capa1_e_avg = squeeze(mean(capa1_e,1));
capa2_e_avg = squeeze(mean(capa2_e,1));
capa5_e_avg = squeeze(mean(capa5_e,1));

% Ergodic Secrecy Rate
sr1_avg = squeeze(mean(sr1));%capa1_b_correl_avg - capa1_e_avg; 
sr2_avg = squeeze(mean(sr2));%capa1_b_correl_avg - capa2_e_avg;
sr5_avg = squeeze(mean(sr5));%capa1_b_correl_avg - capa2_e_avg;

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



% Decoder 2: Matched Filter
sinr2_e_simu = e_avg_sym_decod2_e./e_avg_denom_decod2_e;
sinr2_e_simu_square = sinr2_e_simu.^2;
expected2_e_simu_square = e_avg_sym_decod2_e_exp4./e_avg_denom_decod2_e_exp4;
var2_e_simu = expected2_e_simu_square - sinr2_e_simu_square;

capa2_e_ergodic = capa2_e_avg;
capa2_e_jensen_simu = log2(1+sinr2_e_simu);
capa2_e_new_approx_simu = log2(1+sinr2_e_simu) - 1/2*var2_e_simu./((1+sinr2_e_simu).^2);

% Decoder 5 : Own channel knwoledge
sinr5_e_simu = e_avg_sym_decod5_e./e_avg_denom_decod5_e;
sinr5_e_simu_square = sinr5_e_simu.^2;
expected5_e_simu_square = e_avg_sym_decod5_e_exp4./e_avg_denom_decod5_e_exp4;
var5_e_simu = expected5_e_simu_square - sinr5_e_simu_square;

capa5_e_ergodic = capa5_e_avg;
capa5_e_jensen_simu = log2(1+sinr5_e_simu);
capa5_e_new_approx_simu = log2(1+sinr5_e_simu) - 1/2*var5_e_simu./((1+sinr5_e_simu).^2);


h = waitbar(0,'Modeling SINR/VARIANCE/CAPA at Bob and Eve ...');
tic;
for bb = 1:length(U)
    for dd = 1:size(b_subcar,2)
        parfor nn = 1:length(snr_b)
            sinr_b_model(:,bb,dd,nn)    = sinrModelingFrequencyCorrelation(alpha,U(bb),N(bb),T(:,:,bb,dd),snr_b(nn),snr_e,"bob_correl").';
            var_b_model(:,bb,dd,nn)     = modelVariance(alpha,N(bb),U(bb),T(:,:,bb,dd),sinr_b_model(:,bb,dd,nn),snr_b(nn),snr_e,"bob_correl");
            capa_b_jensen(:,bb,dd,nn)   = log2(1+sinr_b_model(:,bb,dd,nn));
            capa_b_new_approx(:,bb,dd,nn)  = log2(1+sinr_b_model(:,bb,dd,nn)) - 1/2*var_b_model(:,bb,dd,nn)./((1+sinr_b_model(:,bb,dd,nn)).^2);
            
            sinr1_e_model(:,bb,dd,nn)    = sinrModelingFrequencyCorrelation(alpha,U(bb),N(bb),T(:,:,bb,dd),snr_b(nn),snr_e,"eve_decod1_correl").';
            var1_e_model(:,bb,dd,nn)     = modelVariance(alpha,N(bb),U(bb),T(:,:,bb,dd),sinr1_e_model(:,bb,dd,nn),snr_b(nn),snr_e,"eve_decod1_correl");
            capa1_e_jensen(:,bb,dd,nn)       = log2(1+sinr1_e_model(:,bb,dd,nn));
            capa1_e_new_approx(:,bb,dd,nn)   = log2(1+sinr1_e_model(:,bb,dd,nn)) - 1/2*var1_e_model(:,bb,dd,nn)./((1+sinr1_e_model(:,bb,dd,nn)).^2);
        end
    end
    waitbar(bb / length(U))
end
toc;
close(h)



capa_b_ergodic      = capa1_b_avg;
capa_b_jensen       = squeeze(capa_b_jensen);
capa_b_new_approx   = squeeze(capa_b_new_approx);


capa1_e_ergodic      = capa1_e_avg;
capa1_e_jensen       = squeeze(capa1_e_jensen);
capa1_e_new_approx   = squeeze(capa1_e_new_approx);

% ANALYTIC ERRORS
error_1st_order_b = capa_b_ergodic - capa_b_jensen;
error_2nd_order_b = capa_b_ergodic - capa_b_new_approx;

error_1st_order_1e = capa1_e_ergodic - capa1_e_jensen;
error_2nd_order_1e = capa1_e_ergodic - capa1_e_new_approx;


% SIMULATION ERRORS
% error_1st_order_1e_simu = capa1_e_ergodic - capa1_e_jensen_simu;          % IDEM QUE MODEL
% error_2nd_order_1e_simu = capa1_e_ergodic - capa1_e_new_approx_simu;

error_1st_order_2e_simu = capa2_e_ergodic - capa2_e_jensen_simu;
error_2nd_order_2e_simu = capa2_e_ergodic - capa2_e_new_approx_simu;

error_1st_order_5e_simu = capa5_e_ergodic - capa5_e_jensen_simu;
error_2nd_order_5e_simu = capa5_e_ergodic - capa5_e_new_approx_simu;

% SR

sr1_jensen = capa_b_jensen - capa1_e_jensen;
sr1_new_approx = capa_b_new_approx - capa1_e_new_approx;
sr1_ergodic = capa_b_ergodic - capa1_e_ergodic;


sr2_jensen = capa_b_jensen - capa2_e_jensen_simu;
sr2_new_approx = capa_b_new_approx - capa2_e_new_approx_simu;
sr2_ergodic = sr2_avg;

sr5_jensen = capa_b_jensen - capa5_e_jensen_simu;
sr5_new_approx = capa_b_new_approx - capa5_e_new_approx_simu;
sr5_ergodic = sr5_avg;


%% PLOT SECTION 

% Error matrices depend on : correlation , alpha , BOR , SNR --> 4 types of
% possible plots. For each type, 4 possible plots (error as a function of a
% particular parameter, by fixing the 3 others) --> 16 possible plots. 

% When alpha is fixed: alpha = 0.5
% When U is fixed: U = 4
% When correlation is fixed: correl = minimal
% When SNRb is fixed, SNRb = 10 dB

% Default parameters
set(0,'defaultAxesFontSize',15)
set(0,'DefaultLineMarkerSize',8);
set(0,'defaultLineLineWidth',1.4);
aa = find(alpha == 0.5);
bb = find(U == 4);
dd = find(min(b_subcar(1,:)));
% nn = find(snr_b == 10);


%% 1 ERROR PLOTS
for nn = 1:length(snr_b)
figure; 
subplot(4,3,1) 
plot(100*(1-alpha),squeeze(error_1st_order_b(:,bb,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(100*(1-alpha),squeeze(error_2nd_order_b(:,bb,dd,nn)),'Marker','square','color','b'); hold on; 
plot(100*(1-alpha),squeeze(error_1st_order_1e(:,bb,dd,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(100*(1-alpha),squeeze(error_2nd_order_1e(:,bb,dd,nn)),'Marker','v','color','r');
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
l1 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l1, 'location', 'bestoutside')
box on; grid on; 
xlabel('Percentage of AN energy injected : \%')
ylabel('Error : bit/channel use')

   
subplot(4,3,4) 
plot(x_axis(1,:),squeeze(error_1st_order_b(aa,bb,:,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(x_axis(1,:),squeeze(error_2nd_order_b(aa,bb,:,nn)),'Marker','square','color','b'); hold on; 
plot(x_axis(1,:),squeeze(error_1st_order_1e(aa,bb,:,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(x_axis(1,:),squeeze(error_2nd_order_1e(aa,bb,:,nn)),'Marker','v','color','r'); 
title(['BOR = ', num2str(U(bb)), ', $\alpha$ = ', num2str(alpha(aa)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
l2 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l2, 'location', 'bestoutside')
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('Error : bit/channel use')

subplot(4,3,7) 
plot(snr_b,squeeze(error_1st_order_b(aa,bb,dd,:)),'Marker','o','color','b','linestyle',':'); hold on;
plot(snr_b,squeeze(error_2nd_order_b(aa,bb,dd,:)),'Marker','square','color','b'); hold on; 
plot(snr_b,squeeze(error_1st_order_1e(aa,bb,dd,:)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(snr_b,squeeze(error_2nd_order_1e(aa,bb,dd,:)),'Marker','v','color','r'); 
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', $\alpha$ = ', num2str(alpha(aa)) ])
l3 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l3, 'location', 'bestoutside')
box on; grid on; 
xlabel('Bob SNR : dB')
ylabel('Error : bit/channel use')


subplot(4,3,10) 
plot(U,squeeze(error_1st_order_b(aa,:,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(U,squeeze(error_2nd_order_b(aa,:,dd,nn)),'Marker','square','color','b'); hold on; 
plot(U,squeeze(error_1st_order_1e(aa,:,dd,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(U,squeeze(error_2nd_order_1e(aa,:,dd,nn)),'Marker','v','color','r'); 
title(['$\alpha$ = ', num2str(alpha(aa)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])

l4 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l4, 'location', 'bestoutside')
box on; grid on; 
xlabel('BOR')
ylabel('Error : bit/channel use')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,3,2) 
plot(100*(1-alpha),squeeze(error_1st_order_b(:,bb,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(100*(1-alpha),squeeze(error_2nd_order_b(:,bb,dd,nn)),'Marker','square','color','b'); hold on; 
plot(100*(1-alpha),squeeze(error_1st_order_2e_simu(:,bb,dd,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(100*(1-alpha),squeeze(error_2nd_order_2e_simu(:,bb,dd,nn)),'Marker','v','color','r');
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
l1 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l1, 'location', 'bestoutside')
box on; grid on; 
xlabel('Percentage of AN energy injected : \%')
ylabel('Error : bit/channel use')

   
subplot(4,3,5) 
plot(x_axis(1,:),squeeze(error_1st_order_b(aa,bb,:,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(x_axis(1,:),squeeze(error_2nd_order_b(aa,bb,:,nn)),'Marker','square','color','b'); hold on; 
plot(x_axis(1,:),squeeze(error_1st_order_2e_simu(aa,bb,:,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(x_axis(1,:),squeeze(error_2nd_order_2e_simu(aa,bb,:,nn)),'Marker','v','color','r'); 
title(['BOR = ', num2str(U(bb)), ', $\alpha$ = ', num2str(alpha(aa)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
l2 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l2, 'location', 'bestoutside')
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('Error : bit/channel use')

subplot(4,3,8) 
plot(snr_b,squeeze(error_1st_order_b(aa,bb,dd,:)),'Marker','o','color','b','linestyle',':'); hold on;
plot(snr_b,squeeze(error_2nd_order_b(aa,bb,dd,:)),'Marker','square','color','b'); hold on; 
plot(snr_b,squeeze(error_1st_order_2e_simu(aa,bb,dd,:)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(snr_b,squeeze(error_2nd_order_2e_simu(aa,bb,dd,:)),'Marker','v','color','r'); 
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', $\alpha$ = ', num2str(alpha(aa)) ])
l3 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l3, 'location', 'bestoutside')
box on; grid on; 
xlabel('Bob SNR : dB')
ylabel('Error : bit/channel use')


subplot(4,3,11) 
plot(U,squeeze(error_1st_order_b(aa,:,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(U,squeeze(error_2nd_order_b(aa,:,dd,nn)),'Marker','square','color','b'); hold on; 
plot(U,squeeze(error_1st_order_2e_simu(aa,:,dd,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(U,squeeze(error_2nd_order_2e_simu(aa,:,dd,nn)),'Marker','v','color','r'); 
title(['$\alpha$ = ', num2str(alpha(aa)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])

l4 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l4, 'location', 'bestoutside')
box on; grid on; 
xlabel('BOR')
ylabel('Error : bit/channel use')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subplot(4,3,3) 
plot(100*(1-alpha),squeeze(error_1st_order_b(:,bb,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(100*(1-alpha),squeeze(error_2nd_order_b(:,bb,dd,nn)),'Marker','square','color','b'); hold on; 
plot(100*(1-alpha),squeeze(error_1st_order_5e_simu(:,bb,dd,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(100*(1-alpha),squeeze(error_2nd_order_5e_simu(:,bb,dd,nn)),'Marker','v','color','r');
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
l1 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l1, 'location', 'bestoutside')
box on; grid on; 
xlabel('Percentage of AN energy injected : \%')
ylabel('Error : bit/channel use')

   
subplot(4,3,6) 
plot(x_axis(1,:),squeeze(error_1st_order_b(aa,bb,:,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(x_axis(1,:),squeeze(error_2nd_order_b(aa,bb,:,nn)),'Marker','square','color','b'); hold on; 
plot(x_axis(1,:),squeeze(error_1st_order_5e_simu(aa,bb,:,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(x_axis(1,:),squeeze(error_2nd_order_5e_simu(aa,bb,:,nn)),'Marker','v','color','r'); 
title(['BOR = ', num2str(U(bb)), ', $\alpha$ = ', num2str(alpha(aa)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
l2 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l2, 'location', 'bestoutside')
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('Error : bit/channel use')

subplot(4,3,9) 
plot(snr_b,squeeze(error_1st_order_b(aa,bb,dd,:)),'Marker','o','color','b','linestyle',':'); hold on;
plot(snr_b,squeeze(error_2nd_order_b(aa,bb,dd,:)),'Marker','square','color','b'); hold on; 
plot(snr_b,squeeze(error_1st_order_5e_simu(aa,bb,dd,:)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(snr_b,squeeze(error_2nd_order_5e_simu(aa,bb,dd,:)),'Marker','v','color','r'); 
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', $\alpha$ = ', num2str(alpha(aa)) ])
l3 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l3, 'location', 'bestoutside')
box on; grid on; 
xlabel('Bob SNR : dB')
ylabel('Error : bit/channel use')


subplot(4,3,12) 
plot(U,squeeze(error_1st_order_b(aa,:,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(U,squeeze(error_2nd_order_b(aa,:,dd,nn)),'Marker','square','color','b'); hold on; 
plot(U,squeeze(error_1st_order_5e_simu(aa,:,dd,nn)),'Marker','diamond','color','r','linestyle',':'); hold on; 
plot(U,squeeze(error_2nd_order_5e_simu(aa,:,dd,nn)),'Marker','v','color','r'); 
title(['$\alpha$ = ', num2str(alpha(aa)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])

l4 =['Bob : 1st order'; 'Bob : 2nd order' ; ...
       'Eve : 1st order' ; 'Eve : 2nd order']; 
legend(l4, 'location', 'bestoutside')
box on; grid on; 
xlabel('BOR')
ylabel('Error : bit/channel use')


end

 

%% 2. SR PLOTS



for nn = 1:length(snr_b)
figure; 
subplot(4,3,1) 
plot(100*(1-alpha),squeeze(sr1_jensen(:,bb,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(100*(1-alpha),squeeze(sr1_new_approx(:,bb,dd,nn)),'Marker','square','color','g'); hold on; 
plot(100*(1-alpha),squeeze(sr1_ergodic(:,bb,dd,nn)),'Marker','diamond','color','r','linestyle','-.');  
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])

legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('Percentage of AN energy injected : \%')
ylim([-3 max(squeeze(sr1_jensen(:,bb,dd,nn)))+2])
ylabel('SR: bit/ch use')
   
subplot(4,3,4) 
plot(x_axis(1,:),squeeze(sr1_jensen(aa,bb,:,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(x_axis(1,:),squeeze(sr1_new_approx(aa,bb,:,nn)),'Marker','square','color','g'); hold on; 
plot(x_axis(1,:),squeeze(sr1_ergodic(aa,bb,:,nn)),'Marker','diamond','color','r','linestyle','-.'); 
title(['BOR = ', num2str(U(bb)), ', $\alpha$ = ', num2str(alpha(aa)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('SR: bit/ch use')

subplot(4,3,7) 
plot(snr_b,squeeze(sr1_jensen(aa,bb,dd,:)),'Marker','o','color','b','linestyle',':'); hold on;
plot(snr_b,squeeze(sr1_new_approx(aa,bb,dd,:)),'Marker','square','color','g'); hold on; 
plot(snr_b,squeeze(sr1_ergodic(aa,bb,dd,:)),'Marker','diamond','color','r','linestyle','-.'); hold on; 
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', $\alpha$ = ', num2str(alpha(aa)) ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('Bob SNR : dB')
ylabel('SR: bit/ch use')

subplot(4,3,10) 
plot(U,squeeze(sr1_jensen(aa,:,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(U,squeeze(sr1_new_approx(aa,:,dd,nn)),'Marker','square','color','g'); hold on; 
plot(U,squeeze(sr1_ergodic(aa,:,dd,nn)),'Marker','diamond','color','r','linestyle','-.'); 
title(['$\alpha$ = ', num2str(alpha(aa)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])

legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('BOR')
ylabel('SR: bit/ch use')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(4,3,2) 
plot(100*(1-alpha),squeeze(sr2_jensen(:,bb,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(100*(1-alpha),squeeze(sr2_new_approx(:,bb,dd,nn)),'Marker','square','color','g'); hold on; 
plot(100*(1-alpha),squeeze(sr2_ergodic(:,bb,dd,nn)),'Marker','diamond','color','r','linestyle','-.'); hold on; 
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
ylim([-3 max(squeeze(sr1_jensen(:,bb,dd,nn)))+2])
xlabel('Percentage of AN energy injected : \%')
ylabel('SR: bit/ch use')
   
subplot(4,3,5) 
plot(x_axis(1,:),squeeze(sr2_jensen(aa,bb,:,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(x_axis(1,:),squeeze(sr2_new_approx(aa,bb,:,nn)),'Marker','square','color','g'); hold on; 
plot(x_axis(1,:),squeeze(sr2_ergodic(aa,bb,:,nn)),'Marker','diamond','color','r','linestyle','-.');
title(['BOR = ', num2str(U(bb)), ', $\alpha$ = ', num2str(alpha(aa)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('SR: bit/ch use')

subplot(4,3,8) 
plot(snr_b,squeeze(sr2_jensen(aa,bb,dd,:)),'Marker','o','color','b','linestyle',':'); hold on;
plot(snr_b,squeeze(sr2_new_approx(aa,bb,dd,:)),'Marker','square','color','g'); hold on; 
plot(snr_b,squeeze(sr2_ergodic(aa,bb,dd,:)),'Marker','diamond','color','r','linestyle','-.');
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', $\alpha$ = ', num2str(alpha(aa)) ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('Bob SNR : dB')
ylabel('SR: bit/ch use')

subplot(4,3,11) 
plot(U,squeeze(sr2_jensen(aa,:,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(U,squeeze(sr2_new_approx(aa,:,dd,nn)),'Marker','square','color','g'); hold on; 
plot(U,squeeze(sr2_ergodic(aa,:,dd,nn)),'Marker','diamond','color','r','linestyle','-.'); 
title(['$\alpha$ = ', num2str(alpha(aa)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('BOR')
ylabel('SR: bit/ch use')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(4,3,3) 
plot(100*(1-alpha),squeeze(sr2_jensen(:,bb,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(100*(1-alpha),squeeze(sr2_new_approx(:,bb,dd,nn)),'Marker','square','color','g'); hold on; 
plot(100*(1-alpha),squeeze(sr2_ergodic(:,bb,dd,nn)),'Marker','diamond','color','r','linestyle','-.'); hold on; 
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
ylim([-3 max(squeeze(sr1_jensen(:,bb,dd,nn)))+2])
xlabel('Percentage of AN energy injected : \%')
ylabel('SR: bit/ch use')
   
subplot(4,3,6) 
plot(x_axis(1,:),squeeze(sr2_jensen(aa,bb,:,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(x_axis(1,:),squeeze(sr2_new_approx(aa,bb,:,nn)),'Marker','square','color','g'); hold on; 
plot(x_axis(1,:),squeeze(sr2_ergodic(aa,bb,:,nn)),'Marker','diamond','color','r','linestyle','-.');
title(['BOR = ', num2str(U(bb)), ', $\alpha$ = ', num2str(alpha(aa)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('$\Delta f_C / \Delta f_N$')
ylabel('SR: bit/ch use')

subplot(4,3,9) 
plot(snr_b,squeeze(sr2_jensen(aa,bb,dd,:)),'Marker','o','color','b','linestyle',':'); hold on;
plot(snr_b,squeeze(sr2_new_approx(aa,bb,dd,:)),'Marker','square','color','g'); hold on; 
plot(snr_b,squeeze(sr2_ergodic(aa,bb,dd,:)),'Marker','diamond','color','r','linestyle','-.');
title(['BOR = ', num2str(U(bb)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', $\alpha$ = ', num2str(alpha(aa)) ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('Bob SNR : dB')
ylabel('SR: bit/ch use')

subplot(4,3,12) 
plot(U,squeeze(sr2_jensen(aa,:,dd,nn)),'Marker','o','color','b','linestyle',':'); hold on;
plot(U,squeeze(sr2_new_approx(aa,:,dd,nn)),'Marker','square','color','g'); hold on; 
plot(U,squeeze(sr2_ergodic(aa,:,dd,nn)),'Marker','diamond','color','r','linestyle','-.'); 
title(['$\alpha$ = ', num2str(alpha(aa)), ', $\Delta f_C / \Delta f_N$ = ', num2str(x_axis(1,dd)), ...
    ', SNRb = ', num2str(snr_b(nn)), ' dB' ])
legend('1st order SR', '2nd order SR', 'Ergodic SR','location', 'bestoutside')
box on; grid on; 
xlabel('BOR')
ylabel('SR: bit/ch use')
end

