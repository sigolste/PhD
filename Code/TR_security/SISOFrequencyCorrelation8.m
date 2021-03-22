%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   Code that simulates the effect of correlation introduced at B and E on 
%   SR performances. Here, we fix the signal BW, the number of subcarriers,
%   therefore the subcarrier BW, the BOR, and we change the coherence BW
%   in order 
%
%
%   Code Started:   19.02.2021
%   Last Update:    22.02.2021
%
%  
%   Mathematical derivations: cfr FC35-36
%
%  © SIDNEY GOLSTEIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% close all;

warning('off','all')
rmpath('folderthatisnotonpath')

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
nb_run = 1500;                               % number of experiments
fc = 2e9 ;                                  % Carrier frequency
c = 3e8;                                    % Light speed

alpha_step = 5;                           % Percentage between subsequent alpha values
alpha = 0:alpha_step/100:1;         

% Communication parameters
Q = 16;
U = [8];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters

snr_b  = [15].';  % SNR @Bob
snr_e  = [15];  % SNR @Eve

% Noise energy
sigma_b = 1./U./10.^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U./10.^(snr_e/10);     % expected noise energy @Eve


%% Channel modeling

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance
sigma_tau = .5e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                        % Approximation of coherence bandwidth

% Correlation @ Bob
coef_freq_b = [1:.5:3].*N.'/6;                              % From min/6 ∆fc to max/6 ∆fc, depending on N
delta_f_n_b = coef_freq_b.*delta_f_c;   
b_subcar_b = delta_f_n_b./N.';                              % Bandwidth of each subcarrier

% Correlation @ Eve
coef_freq_e = coef_freq_b;                                  % From min/6 ∆fc to max/6 ∆fc, depending on N
delta_f_n_e = coef_freq_e.*delta_f_c;
b_subcar_e = delta_f_n_e./N.';                              % Bandwidth of each subcarrier
% 

x_axis_b  = b_subcar_b./delta_f_c;                              % number of coherence bandwidth where lies a subcarrier 
x_axis = x_axis_b(1,:);

%% Matrix instantiation
e_noise_e            = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_b            = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_an_TX              = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

e_sym_decod1_b       = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_decod1_b     = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));


e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

e_sym_decod3_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_decod3_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_an_decod3_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_denom_decod3_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

e_sym_decod4_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_decod4_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_an_decod4_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_denom_decod4_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

e_sym_decod5_e      = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_noise_decod5_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_an_decod5_e       = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
e_denom_decod5_e    = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

sym_b_exp4          = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
sym_decod1_e_exp4   = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
sym_decod2_e_exp4   = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
sym_decod5_e_exp4   = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

noise_b_exp4        = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
denom_decod1_e_exp4 = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
denom_decod2_e_exp4 = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));
denom_decod5_e_exp4 = zeros(nb_run,length(alpha),length(U),size(b_subcar_b,2),size(b_subcar_e,2),length(snr_b));

Hb       = zeros(nb_run,Q,length(U),size(b_subcar_b,2));
Hb1      = zeros(nb_run,Q,length(U),size(b_subcar_b,2));
He       = zeros(nb_run,Q,length(U),size(b_subcar_e,2));
He1      = zeros(nb_run,Q,length(U),size(b_subcar_e,2));
Tb       = zeros(Q,Q,length(U),size(b_subcar_b,2));
Te       = zeros(Q,Q,length(U),size(b_subcar_e,2));
abs_rho = zeros(Q,length(U),size(b_subcar_b,2));

% Channel Generation:
for bb = 1:length(U)
    for dd = 2:size(b_subcar_b,2) % Bob
        [Hb(:,:,bb,dd), Hb1(:,:,bb,dd), abs_rho(:,bb,dd), Tb(:,:,bb,dd)] = corr_frequency( Q , b_subcar_b(bb,dd) , sigma_tau , nb_run ) ;
    end
    for dd = 2:size(b_subcar_e,2) % Eve
        [He(:,:,bb,dd), He1(:,:,bb,dd), ~ , Te(:,:,bb,dd)] = corr_frequency( Q , b_subcar_e(bb,dd) , sigma_tau , nb_run ) ;
    end
end
% figure;
% subplot(1,2,1)


for iter = 1:nb_run

for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


for dd = 1:size(b_subcar_b,2)
for ee = 1:size(b_subcar_e,2)
%channel generation
% Hb_TX = diag(squeeze(Hb1(iter,:,bb,dd)).');          % Variable correlation at Bob
% Hwb_TX = diag(squeeze(Hb(iter,:,bb,dd)).');         
%Hb_TX = 1/sqrt(2)*(randn(1) + 1j*randn(1))*eye(Q);  % Bob fully correlated
if dd == 1
    Hb_TX = 1/sqrt(2)*(randn(1) + 1j*randn(1))*eye(Q);
else
    Hb_TX = diag(squeeze(Hb1(iter,:,bb,dd)).');          % Variable correlation at Eve
end


if ee == 1
    He_TX = 1/sqrt(2)*(randn(1) + 1j*randn(1))*eye(Q);
else
    He_TX = diag(squeeze(He1(iter,:,bb,ee)).');          % Variable correlation at Eve
end
%He_TX = 1/sqrt(2)*(randn(1) + 1j*randn(1))*eye(Q);  % Eve fully correlated

% Hwe_TX = diag(squeeze(He(iter,:,bb,ee)).');
% He_TX = diag(squeeze(He1(iter,:,bb,ee)).');          % Variable correlation at Eve
He_RX = ctranspose(He_TX);
% Hwe_RX = ctranspose(Hwe_TX);

% Hb_TX = channelRayleigh(Q, mu , sigma);             % Assumption: uncorrelated subcarrier for Eve channel
Hb_RX = ctranspose(Hb_TX);


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
[noise_b, e_noise_b(iter,aa,bb,dd,ee,nn) ] = addNoise(sym_b , snr_b(nn), 1/U(bb));

[noise_e, e_noise_e(iter,aa,bb,dd,ee,nn) ] = addNoise(sym_e , snr_e, 1/U(bb));

% AN symbol
an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob


%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_RX*He_TX;                 % matched filter
decod3 = matrix_despread*(Hb_RX/He_RX);                     % AN killer
gamma_E = (He_RX*Hb_TX)*matrix_spread;
gamma_EH = ctranspose(gamma_E);
decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an));   % LMMSE
decod5 = matrix_despread*He_TX;                             % Only He known by Eve
% 
sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_decod2_e = decod2*sym_e;
sym_decod3_e = decod3*sym_e;
sym_decod4_e = decod4*sym_e;
sym_decod5_e = decod5*sym_e;

sym_b_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(sym_decod1_b).^4); 
sym_decod1_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(sym_decod1_e).^4); 
sym_decod2_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(sym_decod2_e).^4); 
sym_decod5_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(sym_decod5_e).^4); 


% 
noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;
noise_b_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(noise_decod1_b).^4); 

noise_decod2_e = decod2*noise_e;
noise_decod3_e = decod3*noise_e;
noise_decod4_e = decod4*noise_e;
noise_decod5_e = decod5*noise_e;
% 
an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;
an_decod3_e = decod3*an_e;
an_decod4_e = decod4*an_e;
an_decod5_e = decod5*an_e;

denom_decod1_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(noise_decod1_e + an_decod1_e).^4);
denom_decod2_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(noise_decod2_e + an_decod2_e).^4);
% denom_decod3_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(noise_decod3_e + an_decod3_e).^4);
% denom_decod4_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(noise_decod4_e + an_decod4_e).^4);
denom_decod5_e_exp4(iter,aa,bb,dd,ee,nn) = mean(abs(noise_decod5_e + an_decod5_e).^4);



%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,aa,bb,dd,ee,nn)      = energy(sym_decod1_b); %  energy(matrix_despread*abs(Hb_RX).^2*matrix_spread);
e_noise_decod1_b(iter,aa,bb,dd,ee,nn)    = energy(noise_decod1_b);



% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,dd,ee,nn)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,dd,ee,nn)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,aa,bb,dd,ee,nn)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb,dd,ee,nn)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

% @ Eve : decod 3
e_sym_decod3_e(iter,aa,bb,dd,ee,nn)     = energy(sym_decod3_e);
e_noise_decod3_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod3_e);
e_an_decod3_e(iter,aa,bb,dd,ee,nn)      = energy(an_decod3_e);
e_denom_decod3_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod3_e + an_decod3_e);        % energy of the sinr denominator for decoder 3 @Eve

% @ Eve : decod 4
e_sym_decod4_e(iter,aa,bb,dd,ee,nn)     = energy(sym_decod4_e);
e_noise_decod4_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod4_e);
e_an_decod4_e(iter,aa,bb,dd,ee,nn)      = energy(an_decod4_e);
e_denom_decod4_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod4_e + an_decod4_e);        % energy of the sinr denominator for decoder 4 @Eve

% @ Eve : decod 5
e_sym_decod5_e(iter,aa,bb,dd,ee,nn)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod5_e);
e_an_decod5_e(iter,aa,bb,dd,ee,nn)      = energy(an_decod5_e);
e_denom_decod5_e(iter,aa,bb,dd,ee,nn)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve


end
end
end
waitbar(iter / nb_run)
end
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

e_avg_sym_decod4_e      = squeeze(mean(e_sym_decod4_e,1));
e_avg_noise_decod4_e    = squeeze(mean(e_noise_decod4_e,1));
e_avg_an_decod4_e       = squeeze(mean(e_an_decod4_e,1));
e_avg_denom_decod4_e    = squeeze(mean(e_denom_decod4_e,1));

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
sinr4_e = e_sym_decod4_e./e_denom_decod4_e;
sinr5_e = e_sym_decod5_e./e_denom_decod5_e;



% Instantaneous capacity
capa1_b = log2(1+sinr1_b);

capa1_e = log2(1+sinr1_e);
capa2_e = log2(1+sinr2_e);
capa3_e = log2(1+sinr3_e);
capa4_e = log2(1+sinr4_e);
capa5_e = log2(1+sinr5_e);

% Instantaneous SR
sr1 = secrecyCapacity(sinr1_b,sinr1_e);
sr2 = secrecyCapacity(sinr1_b,sinr2_e);
sr3 = secrecyCapacity(sinr1_b,sinr3_e);
sr4 = secrecyCapacity(sinr1_b,sinr4_e);
sr5 = secrecyCapacity(sinr1_b,sinr5_e);


% Ergodic SINR
sinr1_b_avg      = squeeze(mean(sinr1_b,1)); 

sinr1_e_avg = squeeze(mean(sinr1_e,1)); 
sinr2_e_avg = squeeze(mean(sinr2_e,1)); 
sinr3_e_avg = squeeze(mean(sinr3_e,1)); 
sinr4_e_avg = squeeze(mean(sinr4_e,1)); 
sinr5_e_avg = squeeze(mean(sinr5_e,1)); 


% Ergodic capacity
capa1_b_avg = squeeze(mean(capa1_b,1)); 

capa1_e_avg = squeeze(mean(capa1_e,1));
capa2_e_avg = squeeze(mean(capa2_e,1));
capa3_e_avg = squeeze(mean(capa3_e,1));
capa4_e_avg = squeeze(mean(capa4_e,1));
capa5_e_avg = squeeze(mean(capa5_e,1));

% Ergodic Secrecy Rate
sr1_avg = squeeze(mean(sr1));%capa1_b_correl_avg - capa1_e_avg; 
sr2_avg = squeeze(mean(sr2));%capa1_b_correl_avg - capa2_e_avg;
sr3_avg = squeeze(mean(sr3));%capa1_b_correl_avg - capa3_e_avg;
sr4_avg = squeeze(mean(sr4));%capa1_b_correl_avg - capa4_e_avg;
sr5_avg = squeeze(mean(sr5));%capa1_b_correl_avg - capa5_e_avg;




sr1_ergodic = sr1_avg;
sr2_ergodic = sr2_avg;
sr3_ergodic = sr3_avg;
sr4_ergodic = sr4_avg;
sr5_ergodic = sr5_avg;

%% Plot Section : max SR as a function of variable correlations at Bob and/or Eve
% a = b_subcar_e./delta_f_c;
a = b_subcar_b./delta_f_c;
figure;
plot(a ,max(sr1_avg),'Marker','o'); hold on;
plot(a ,max(sr2_avg),'Marker','square'); hold on;
plot(a ,max(sr3_avg),'Marker','v'); hold on;
plot(a ,max(sr4_avg),'Marker','<'); hold on;
plot(a ,max(sr5_avg),'Marker','diamond'); 
legend('SDS Decoder','MF Decoder','AN Killer Decoder', 'LMMSE Decoder','OC Decoder','Location','bestoutside')
box on; grid on;
xlim([min(a) max(a)])
xlabel('Bob correlation: $\Delta f_N / \Delta f_C$');
ylabel('Maximal Ergodic Secrecy Rate (bit/channel use)')
title('No correlation at Eve')



a = b_subcar_b./delta_f_c;
figure;
yyaxis left
plot(a ,capa1_b_avg,'Marker','o'); hold on;
ylabel('Capacity (bit/channel use)')
yyaxis right
plot(a,sinr1_b_avg,'Marker','v')
ylabel('SINR')
legend('Ergodic capacity','Ergodic SINR','location','bestoutside')
box on; grid on;
xlim([min(a) max(a)])
xlabel('Bob correlation: $\Delta f_N / \Delta f_C$');
%title('Maxima')




