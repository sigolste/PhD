%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   Implementation of new approx for Bob capa and Eve capa with SDS/MF/OC 
%   decoders when  frequency corrleation among Bob's subcarriers is  
%   introduced.
%
%
%
%   Code Started:   12.01.2021
%   Last Update:    05.02.2021
%
%   © SIDNEY GOLSTEIN
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%close all;

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
nb_run = 3000;                               % number of experiments
fc = 2e9 ;                                  % Carrier frequency
c = 3e8;                                    % Light speed

alpha_step = 5;                           % Percentage between subsequent alpha values
alpha = .5;%0:alpha_step/100:1;         

% Communication parameters
Q = 4;
U = 4;
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters

snr_b  = 10;  % SNR @Bob
snr_e  = 200;  % SNR @Eve

% Noise energy
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve


%% Channel modeling

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance
sigma_tau = .5e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                        % Approximation of coherence bandwidth
coef_freq = [1:5:30].*N/6;
delta_f_n = coef_freq.*delta_f_c;   
b_subcar = delta_f_n./N;                                    % Bandwidth of each subcarrier
x_axis  = delta_f_n./delta_f_c;
for dd = 1:length(b_subcar)
    [H(:,:,dd), H1(:,:,dd), abs_rho(:,dd), T(:,:,dd)] = corr_frequency( Q , b_subcar(dd) , sigma_tau , nb_run ) ;
end

%% Energy matrix instantiation
e_noise_e            = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_noise_b            = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_an_TX              = zeros(nb_run,length(alpha),length(U),length(b_subcar));

e_sym_decod1_b       = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_noise_decod1_b     = zeros(nb_run,length(alpha),length(U),length(b_subcar));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar));


e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),length(b_subcar));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar));

for iter = 1:nb_run

for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


for dd = 1:length(b_subcar)
%channel generation
Hb_TX = diag(squeeze(H1(iter,:,dd)).');
Hw_TX = diag(squeeze(H(iter,:,dd)).');

Hb_RX = ctranspose(Hb_TX);
Hw_RX = ctranspose(Hw_TX);

He_TX = channelRayleigh(Q, mu , sigma);             % Assumption: uncorrelated subcarrier for Eve channel
He_RX = ctranspose(He_TX);


%% Encoder

% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % Qx1


sym_precoded = Hb_TX*sym_spread;


% AN generation
[an,V1_correl(:,:,iter),S_correl(:,:,iter)] = generateAN_TEST(Hb_RX,Q,U(bb),matrix_despread,1/U(bb),"svd");% % Qx1, not weighted 


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




% Noise symbol
[noise_b, e_noise_b(iter,aa,bb,dd) ] = addNoise(sym_b , snr_b, 1/U(bb));

[noise_e, e_noise_e(iter,aa,bb,dd) ] = addNoise(sym_e , snr_e, 1/U(bb));

% AN symbol
an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob


%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_RX*He_TX;                 % matched filter
% decod3 = matrix_despread*(Hb_RX/He_RX);                     % AN killer
% gamma_E = (He_RX*Hb_TX)*matrix_spread;
% gamma_EH = ctranspose(gamma_E);
% decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an) + e_noise_e*eye(Q) );   % LMMSE
% decod5 = matrix_despread*He_TX;                             % Only He known by Eve
% 
sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_b_test(:,iter,aa,bb,dd) = sym_decod1_b; 
sym_decod1_e_test(:,iter,aa,bb,dd) = sym_decod1_e; 

% 
noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;
noise_b_test(:,iter,aa,bb,dd) = noise_decod1_b; 

noise_decod2_e = decod2*noise_e;
% noise_decod3_e = decod3*noise_e;
% noise_decod4_e = decod4*noise_e;
% noise_decod5_e = decod5*noise_e;
% 
an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;

denom_decod1_e(:,iter,aa,bb,dd) = noise_decod1_e + an_decod1_e;

% an_decod3_e = decod3*an_e;
% an_decod4_e = decod4*an_e;
% an_decod5_e = decod5*an_e;

%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,aa,bb,dd)      = energy(sym_decod1_b); %  energy(matrix_despread*abs(Hb_RX).^2*matrix_spread);
e_noise_decod1_b(iter,aa,bb,dd)    = energy(noise_decod1_b);



% @ Eve : decod 1
e_sym_decod1_e(iter,aa,bb,dd)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,aa,bb,dd)   = energy(noise_decod1_e);
e_an_decod1_e(iter,aa,bb,dd)      = energy(an_decod1_e);
e_denom_decod1_e(iter,aa,bb,dd)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
% e_sym_decod2_e(iter,aa,bb,dd)     = energy(sym_decod2_e);
% e_noise_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e);
% e_an_decod2_e(iter,aa,bb,dd)      = energy(an_decod2_e);
% e_denom_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

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
% % @ Eve : decod 5
% e_sym_decod5_e(iter,aa,bb,dd)     = energy(sym_decod5_e);
% e_noise_decod5_e(iter,aa,bb,dd)   = energy(noise_decod5_e);
% e_an_decod5_e(iter,aa,bb,dd)      = energy(an_decod5_e);
% e_denom_decod5_e(iter,aa,bb,dd)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve



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

% e_avg_sym_decod2_e      = squeeze(mean(e_sym_decod2_e,1));
% e_avg_noise_decod2_e    = squeeze(mean(e_noise_decod2_e,1));
% e_avg_an_decod2_e       = squeeze(mean(e_an_decod2_e,1));
% e_avg_denom_decod2_e    = squeeze(mean(e_denom_decod2_e,1));


% Instantaneous SINR
sinr1_b = e_sym_decod1_b./e_noise_decod1_b;

sinr1_e = e_sym_decod1_e./e_denom_decod1_e;
% sinr2_e = e_sym_decod2_e./e_denom_decod2_e;



% Instantaneous capacity
capa1_b = log2(1+sinr1_b);


capa1_e = log2(1+sinr1_e);
% capa2_e = log2(1+sinr2_e);


sr1 = secrecyCapacity(sinr1_b,sinr1_e);
% sr2 = secrecyCapacity(sinr1_b,sinr2_e);


% Ergodic SINR
sinr1_b_avg      = squeeze(mean(sinr1_b,1)); 

sinr1_e_avg = squeeze(mean(sinr1_e,1)); 
% sinr2_e_avg = squeeze(mean(sinr2_e,1)); 


% Ergodic capacity
capa1_b_avg = squeeze(mean(capa1_b,1)); 

capa1_e_avg = squeeze(mean(capa1_e,1));
% capa2_e_avg = squeeze(mean(capa2_e));

% Ergodic Secrecy Rate
sr1_avg = squeeze(mean(sr1));%capa1_b_correl_avg - capa1_e_avg; 
% sr2_avg = squeeze(mean(sr2));%capa1_b_correl_avg - capa2_e_avg;




H4H4_test = 2*U; % OK
% % 
% % % 2 ∑_i ∑_j!=i   |hb_{n+iN}|^2*|he_{n+iN}|^2*|hb_{n+jN}|^2*|he_{n+jiN}|^2
% H2H2 = mean(h2h2); % OK
for dd = 1:length(b_subcar)
H2H2_test(:,dd) = modelCorrelH2H2(N,U,T(:,:,dd)); % Ok
end

% 
for dd = 1:length(b_subcar)   
sinr_b_model(:,dd)    = sinrModelingFrequencyCorrelation(alpha,U,N,T(:,:,dd),snr_b,snr_e,"bob_correl").';
var_b_model(:,dd)     = modelVariance(alpha,N,U,T(:,:,dd),sinr_b_model(:,dd),snr_b,snr_e,"bob_correl");
capa_b_jensen(:,dd)       = log2(1+sinr_b_model(:,dd));
capa_b_new_approx(:,dd)   = log2(1+sinr_b_model(:,dd)) - 1/2*var_b_model(:,dd)./((1+sinr_b_model(:,dd)).^2);
end
capa_b_ergodic      = capa1_b_avg;


error_jensen_b = capa1_e_ergodic - capa_b_jensen;
error_new_approx_b = capa1_e_ergodic - capa_b_new_approx;


   
%% Eve Decoder 1
for dd = 1:length(b_subcar)   
sinr1_e_model(:,dd)    = sinrModelingFrequencyCorrelation(alpha,U,N,T(:,:,dd),snr_b,snr_e,"eve_decod1_correl").';
var1_e_model(:,dd)     = modelVariance(alpha,N,U,T(:,:,dd),sinr1_e_model(:,dd),snr_b,snr_e,"eve_decod1_correl");
capa1_e_jensen(:,dd)       = log2(1+sinr1_e_model(:,dd));
capa1_e_new_approx(:,dd)   = log2(1+sinr1_e_model(:,dd)) - 1/2*var1_e_model(:,dd)./((1+sinr1_e_model(:,dd)).^2);
end
capa1_e_ergodic      = capa1_e_avg;


%% SR 
sr1_jensen = capa_b_jensen - capa1_e_jensen;
sr1_new_approx = capa_b_new_approx - capa1_e_new_approx;
sr1_ergodic = capa_b_ergodic - capa1_e_ergodic;


%% PLOT SECTION
figure; 
plot(capa_b_ergodic, 'Marker', 'square'); 
hold on; plot(capa_b_jensen, 'Marker', 'o') ; 
plot(capa_b_new_approx, 'Marker', 'v');
if length(b_subcar) == 2
legend('Ergodic capacity: $E[\log_2(1+X)]$, correl','Ergodic capacity: $E[\log_2(1+X)]$, no correl', ...
        'Jensen inequality: $\log_2(1+E[X])$, correl', 'Jensen inequality: $\log_2(1+E[X])$, no correl',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$, correl', ...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$, no correl','location','best')
else
    legend('Ergodic capacity: $E[\log_2(1+X)]$', ...
        'Jensen inequality: $\log_2(1+E[X])$',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
end
title('BOB')
figure; 
plot(capa1_e_ergodic, 'Marker', 'square'); 
hold on; plot(capa1_e_jensen, 'Marker', 'o') ; 
plot(capa1_e_new_approx, 'Marker', 'v');
if length(b_subcar) == 2
legend('Ergodic capacity: $E[\log_2(1+X)]$, correl','Ergodic capacity: $E[\log_2(1+X)]$, no correl', ...
        'Jensen inequality: $\log_2(1+E[X])$, correl', 'Jensen inequality: $\log_2(1+E[X])$, no correl',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$, correl', ...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$, no correl','location','best')
else
    legend('Ergodic capacity: $E[\log_2(1+X)]$', ...
        'Jensen inequality: $\log_2(1+E[X])$',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
end   
title('EVE')

figure; 
plot(sr1_ergodic,'Marker', 'square'); 
hold on; plot(sr1_jensen, 'Marker', 'o') ; 
plot(sr1_new_approx,'Marker', 'v');
if length(b_subcar) == 2
legend('Ergodic SR, correl','Ergodic SR, no correl', ...
        'Jensen inequality SR correl', 'Jensen inequality SR, no correl',...
       'New approx SR, correl', ...
       'New approx SR, no correl','location','best')
else
    legend('Ergodic SR', ...
        'Jensen inequality SR',...
       'New approx SR','location','best')
end
   

