clear all;
close all; 
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',32)
set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'DefaultLineMarkerSize',15);
set(0, 'defaultFigurePosition',  [-1267  44   1256    872])
% figure('WindowStyle','Docked')
% set(0,'defaultMarkerIndices',[1:2:end])
%**************************************************************************
%
%   This script simulates an implementation of a secure Time Reversal
%   scheme in the frequency domain using OFDM. 
%
%   Frequency domain secure TR in a MISO configuration    
%   Complex Envelope (CE) simulations
%
%   OFDM to mitigate MPC effects
%
%   Artificial Noise added to secure the data transmission at Bob
%
%   Rice channel model implemented
%
%   Improved decoders:
%   1) Matched filtering at Eve before despreading --> Coherent combination
%   of the phases of RX symbols and then despreading --> TR gain at Eve
%   2) Eve estimates the AN
%
%   last update: 3 April 2019
%
%
%   by Sidney Golstein
%
%**************************************************************************


%% Initialisation of all parameters in function of the scenario. 

%scenario_name of type "yparam_xparam_variableparam_constantparam"
% 
% sc = ["ber_ebno_bor_alpha" "ber_ebno_alpha_bor" "ber_alpha_bor_ebno" "ber_alpha_ebno_bor" ...
%     "ber_bor_ebno_alpha" "ber_bor_alpha_ebno" "secrecy_bor_ebno_alpha" "secrecy_bor_alpha_ebno" ...
%     "secrecy_alpha_bor_ebno" "secrecy_alpha_ebno_bor" "secrecy_ebno_alpha_bor" "secrecy_ebno_bor_alpha"];   % All scenarios
% %     
% sc = ["ber_ebno_bor_alpha" "ber_ebno_alpha_bor" "ber_alpha_bor_ebno" "ber_alpha_ebno_bor" ...
%     "ber_bor_ebno_alpha" "ber_bor_alpha_ebno"]                                                             % BER scenarios 
sc = [ "secrecy_bor_ebno_alpha" "secrecy_bor_alpha_ebno" ...
    "secrecy_alpha_bor_ebno" "secrecy_alpha_ebno_bor" "secrecy_ebno_alpha_bor" "secrecy_ebno_bor_alpha"];   % Secrecy scenarios


% Individual scenarios
%sc = "ber_ebno_bor_alpha";
%sc = "ber_ebno_alpha_bor";
%sc = "ber_alpha_bor_ebno";
%sc = "ber_alpha_ebno_bor";
%sc = "ber_bor_ebno_alpha";
%sc = "ber_bor_alpha_ebno";
%sc = "secrecy_bor_ebno_alpha";
%sc = "secrecy_bor_alpha_ebno";
% sc = "secrecy_alpha_bor_ebno";
%sc = "secrecy_alpha_ebno_bor";
%sc = "secrecy_ebno_alpha_bor";
%sc = "secrecy_ebno_bor_alpha";
sc = "secrecy_alpha_bor_ebno_match_filt1";

for scenar=1:length(sc)
    close all;
    clearvars -except sc scenar
    scenario = sc{scenar};
    get_parameters


%% 1. transmitted symbol sequence. 

for cc = 1:nb_channels
for nb_bor = 1:length(BOR)
msg_TX      = randi( [0 1] , nb_bit(nb_bor) , 1 ) ;                                 % Random bit data stream
sym_TX      = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types

sym_TX      = reshape(sym_TX, [], nb_block) ;                                       % OFDM block per OFDM block


%% ICI: Creation du preambule qu'il faut aussi spread pour estimer le canal de Eve et Bob pour l'égalisation

%% 2. Spreading matrix of dimension [ nb_subcar x nb_symb_per_block ]. 
% Matrix composed of BOR diagonal matrices of dimensions [nb_symb_per_block x nb_symb_per_block]

spread_matrix   = zeros( nb_subcar , nb_symb_per_block(nb_bor) );
spread_code     = reshape( 2.*randi([0 1], 1, nb_subcar).' - 1 , [] , BOR(nb_bor) )/sqrt(BOR(nb_bor));   % Diagonal elements of the spreading matrix. BPSK spreading code


for i = 1:size(spread_code,2)
    spread_matrix(i*nb_symb_per_block(nb_bor)-nb_symb_per_block(nb_bor)+1:i*end/BOR(nb_bor),:) = diag(spread_code(:,i));
end
despread_matrix = ctranspose(spread_matrix); 

%% 3. Symbol stream after spreading
sym_spread = spread_matrix*sym_TX;


%% 4. Channel generation : 
% Extrended Pedestrian Model (epa), normalized or statistical

% 4.1) Legitimate receiver (bob)
H_bob = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"decorrelated");


% 4.2) Eavesdropper (Eve)

H_eve = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"decorrelated");



HBOB(:,cc) = H_bob;
HEVE(:,cc) = H_eve;


H_bob_TX    = ctranspose(H_bob).';
H_eve_TX    = ctranspose(H_eve).';
H_bob_RX    = H_bob;
H_eve_RX    = H_eve;




%% 5.  ARTIFICIAL NOISE INJECTION

W = TR_security_generate_AN(H_bob,nb_subcar,BOR(nb_bor),despread_matrix,"underdetermined"); % Resolution of the underdetermined system


w_signal = ifft(W);

% Normalization of W
tmp = sum(abs(W).^2)/nb_subcar;
energy_sym_useful_TX = mean(sum(abs(H_bob_TX.*sym_spread).^2))/nb_subcar;
W = W/sqrt(tmp)*sqrt(energy_sym_useful_TX);
energy_W = sum(abs(W).^2)/nb_subcar;%sum(abs(W).^2)/nb_subcar;


%% 6. Symbol stream generation
for aa = 1:length(alpha)

% Precoding
sym_useful =sqrt(alpha(aa))*H_bob_TX.*sym_spread;       % Energy --> sqrt needed
AN = sqrt(1 - alpha(aa))*W;
energy_AN = mean(abs(AN).^2);

sym_raditated_noiseless = sym_useful + AN;
energy_radiated(cc,nb_bor,aa) = mean(sum(abs(sym_raditated_noiseless).^2))/nb_subcar; %sum(mean(abs(sym_raditated_noiseless),2).^2)/nb_subcar;


for nn = 1:length(EbN0)



%% 7. Signal to transmit at the antennas
signal_useful_TX = ifft(sym_useful);
signal_AN_TX = ifft(AN);
signal_TX = ifft(sym_raditated_noiseless);



%% PART 3. RECEIVER SIDE

%% 1) OFDM 

% 1.1) Useful symbols
% sym_useful_RX_bob = TR_security_OFDM(signal_useful_TX, nb_subcar, nb_block, length_CP , h_bob );    % Noiseless useful symbols received @ bob
% sym_useful_RX_eve = TR_security_OFDM(signal_useful_TX, nb_subcar, nb_block, length_CP , h_eve );    % Noiseless useful symbols received @ eve

sym_useful_RX_bob = sym_useful.*H_bob_RX;
sym_useful_RX_eve = sym_useful.*H_eve_RX;


% 1.2) AN symbols (only @ Eve)
%sym_AN_RX_eve = TR_security_OFDM(signal_AN_TX, nb_subcar, 1, length_CP , h_eve );            % Noiseless AN symbols received @ eve
sym_AN_RX_eve = AN.*H_eve_RX;            % Noiseless AN symbols received @ eve



% 1.3) useful symbols + AN
% [sym_RX_bob] = TR_security_OFDM(signal_TX, nb_subcar, nb_block, length_CP , h_bob ); % Noiseless received symbols (useful + AN) @ Bob
% [sym_RX_eve] = TR_security_OFDM(signal_TX, nb_subcar, nb_block, length_CP , h_eve ); % Noiseless received symbols (useful + AN) @ Eve

sym_RX_bob = sym_raditated_noiseless.*H_bob_RX;
sym_RX_eve = sym_raditated_noiseless.*H_eve_RX;


%% 2) Noise addition on useful symbols + AN
[sym_RX_noisy_bob,energy_noise_bob] = TR_security_addNoise(sym_RX_bob, EbN0(nn), M, var(sym_raditated_noiseless));   % Noisy received symbols (useful + AN) @ Bob
[sym_RX_noisy_eve,energy_noise_eve] = TR_security_addNoise(sym_RX_eve, EbN0(nn), M, var(sym_raditated_noiseless));   % Noisy received symbols (useful + AN) @ Eve


%% Eve matched filter before despreading: 

% 1) Matched filtering (filt1)
filt1 = despread_matrix*diag(H_bob_RX.*H_eve_TX); 

% 2) AN killed (filt2)
filt2 = despread_matrix*(diag(H_bob_RX)/(diag(H_eve_RX)));


% 3) LMMSE (filt3)
gamma_E = diag(H_eve_RX.*H_bob_TX)*spread_matrix;
gamma_EH = ctranspose(gamma_E);
filt3 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*diag(abs(H_eve_RX).^2)*energy_W + mean(energy_noise_eve)*eye(nb_subcar) );

% Filt1
sym_RX_eve_filt1 =  filt1*sym_RX_eve;
sym_RX_noisy_eve_filt1 = filt1*sym_RX_noisy_eve;
sym_useful_RX_eve_filt1  = filt1*sym_useful_RX_eve;
sym_AN_RX_eve_filt1 = despread_matrix*(H_bob_RX.*abs(H_eve_RX).^2.*AN);

% Filt 2
sym_RX_eve_filt2      = filt2*sym_RX_eve;
sym_useful_RX_eve_filt2  = filt2*sym_useful_RX_eve;
sym_AN_RX_eve_filt2   = filt2*sym_AN_RX_eve;
sym_RX_noisy_eve_filt2 = filt2*sym_RX_noisy_eve;

% Filt 3
sym_RX_eve_filt3      = filt3*sym_RX_eve;
sym_useful_RX_eve_filt3  = filt3*sym_useful_RX_eve;
sym_AN_RX_eve_filt3   = filt3*sym_AN_RX_eve;
sym_RX_noisy_eve_filt3 = filt3*sym_RX_noisy_eve;

% filtertest(:,cc) = mean(abs(diag(filt1*ctranspose(filt1))).^2);


%% 3. Despreading

% 3.1) Useful symbols
sym_useful_RX_despread_bob = despread_matrix*sym_useful_RX_bob; % Despread useful symbols @ Bob
sym_useful_RX_despread_eve = despread_matrix*sym_useful_RX_eve; % Despread useful symbols @ Eve
%sym_useful_RX_despread_eve_filt1 = despread_matrix*sym_useful_eve_filt1;
% sym_useful_RX_despread_eve_filt2 = despread_matrix*sym_useful_eve_filt2;

% 3.2) AN symbols (only @ Eve)
sym_AN_RX_despread_eve = despread_matrix*sym_AN_RX_eve; % Despread AN symbols @ Eve
%sym_AN_RX_despread_eve_filt1 = despread_matrix*sym_AN_RX_eve_filt1;
% sym_AN_RX_despread_eve_filt2 = despread_matrix*sym_AN_RX_eve_filt2;


% 3.3) useful symbols + AN
sym_RX_despread_bob = despread_matrix*sym_RX_noisy_bob; % Depsread symbols (useful + AN) @ Bob
sym_RX_despread_eve = despread_matrix*sym_RX_noisy_eve; % Depsread symbols (useful + AN) @ Eve
%sym_RX_despread_eve_filt1 = despread_matrix*sym_RX_noisy_eve_filt1;
% sym_RX_despread_eve_filt2 = despread_matrix*sym_RX_noisy_eve_filt2;





%% 4. Equalization


% 4.1) Equalization channels (ZF)

% 4.1.1) At bob
H_bob_eq = sum(reshape(abs(H_bob).^2,nb_subcar/BOR(nb_bor),[]),2)/BOR(nb_bor); % Mean of subcarrier contributions for each symbol
%tmp(cc,:) = H_bob_eq;
% 4.1.2) At Eve
H_eve_eq = sum(reshape(H_bob_TX.*H_eve, nb_subcar/BOR(nb_bor),[]),2)/BOR(nb_bor).';
%tmp2(cc,:) = H_eve_eq;

% 4.2) Equalized symbol stream

% 4.2.1 Useful symbols
sym_useful_RX_despread_bob_eq = sym_useful_RX_despread_bob./H_bob_eq;
sym_useful_RX_despread_eve_eq = sym_useful_RX_despread_eve./H_eve_eq;

% 4.1.2 AN symbols (only @ Eve)
sym_AN_RX_despread_eve_eq = sym_AN_RX_despread_eve./H_eve_eq;

% 4.1.3 Useful symbols + AN
sym_RX_despread_bob_eq = sym_RX_despread_bob ./ H_bob_eq;
sym_RX_despread_eve_eq   = sym_RX_despread_eve ./ H_eve_eq ;


%% 5. Metrics

% 5.1) SINR Computation

% 5.1.1) Before equalization

% noise symbols (energy should be identical --> it is)
noise_despread_bob =  sym_RX_despread_bob - sym_useful_RX_despread_bob;
noise_despread_eve =  sym_RX_despread_eve - sym_AN_RX_despread_eve - sym_useful_RX_despread_eve;
noise_eve_filt1 = sym_RX_noisy_eve_filt1 - sym_AN_RX_eve_filt1 - sym_useful_RX_eve_filt1;
noise_eve_filt2 = sym_RX_noisy_eve_filt2 - sym_AN_RX_eve_filt2 - sym_useful_RX_eve_filt2;
noise_eve_filt3 = sym_RX_noisy_eve_filt3 - sym_AN_RX_eve_filt3 - sym_useful_RX_eve_filt3;


energy_noise_despread_bob = mean(mean(abs(noise_despread_bob).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_bob).^2));
energy_noise_despread_eve = mean(mean(abs(noise_despread_eve).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_eve).^2));
energy_noise_eve_filt1 = mean(mean(abs(noise_eve_filt1).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_eve).^2));
energy_noise_eve_filt2 = mean(mean(abs(noise_eve_filt2).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_eve).^2));
energy_noise_eve_filt3 = mean(mean(abs(noise_eve_filt3).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_eve).^2));



% useful symbol
energy_sym_useful_RX_despread_bob = mean(mean(abs(sym_useful_RX_despread_bob).^2));                         %mean(var(sym_useful_RX_despread_bob)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2)); %mean(var(sym_useful_RX_despread_bob)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2));
energy_sym_useful_RX_despread_eve = mean(mean(abs(sym_useful_RX_despread_eve).^2));                         %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));
energy_sym_useful_RX_eve_filt1 = mean(mean(abs(sym_useful_RX_eve_filt1).^2));        %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));
energy_sym_useful_RX_eve_filt2 = mean(mean(abs(sym_useful_RX_eve_filt2).^2));             %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));
energy_sym_useful_RX_eve_filt3 = mean(mean(abs(sym_useful_RX_eve_filt3).^2));             %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));


% AN symbols (only @ Eve)
energy_sym_AN_RX_despread_eve =  mean(abs(sym_AN_RX_despread_eve).^2); %var(sym_AN_RX_despread_eve) ; % 1/nb_subcar*sum(abs(sym_AN_RX_despread_eve).^2); % %
energy_sym_AN_RX_eve_filt1 =  mean(abs(sym_AN_RX_eve_filt1).^2);% mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve_filt1).^2));
energy_sym_AN_RX_eve_filt3 =  mean(abs(sym_AN_RX_eve_filt3).^2);% mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve_filt1).^2));

energy_an(cc,nn,nb_bor,aa) = energy_W;



energy_denom_eve_filt0 =  mean(mean(abs(sym_AN_RX_despread_eve + noise_despread_eve).^2));
energy_denom_eve_filt1 = mean(mean(abs(sym_AN_RX_eve_filt1 + noise_eve_filt1).^2));
energy_denom_eve_filt2 = mean(mean(abs(sym_AN_RX_eve_filt2 + noise_eve_filt2).^2));
energy_denom_eve_filt3 = mean(mean(abs(sym_AN_RX_eve_filt3 + noise_eve_filt3).^2));


% SINR after despread
SINR_bob_despread(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_bob/energy_noise_despread_bob;
SINR_eve_despread(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_eve/energy_denom_eve_filt0;
SINR_eve_filt1(cc,nn,nb_bor,aa) = energy_sym_useful_RX_eve_filt1/energy_denom_eve_filt1;
SINR_eve_filt2(cc,nn,nb_bor,aa) = energy_sym_useful_RX_eve_filt2/energy_denom_eve_filt2;
SINR_eve_filt3(cc,nn,nb_bor,aa) = energy_sym_useful_RX_eve_filt3/energy_denom_eve_filt3;
SINR_eve_filt3(isnan(SINR_eve_filt3)) = 0;

% Numerator (data)
numbob(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_bob;
numfilt0(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_eve; 
numfilt1(cc,nn,nb_bor,aa) = energy_sym_useful_RX_eve_filt1; 
numfilt2(cc,nn,nb_bor,aa) = energy_sym_useful_RX_eve_filt2;
numfilt3(cc,nn,nb_bor,aa) = energy_sym_useful_RX_eve_filt3;


% Noise
noisefilt1(cc,nn,nb_bor,aa) = energy_noise_eve_filt1 ;
noisefilt2(cc,nn,nb_bor,aa) = energy_noise_eve_filt2;
noisefilt3(cc,nn,nb_bor,aa) = energy_noise_eve_filt3;

% AN
anfilt1(cc,nn,nb_bor,aa) = energy_sym_AN_RX_eve_filt1;
anfilt3(cc,nn,nb_bor,aa) = energy_sym_AN_RX_eve_filt3;

% Denominateur (AN + noise)
denomfilt1(cc,nn,nb_bor,aa) = energy_denom_eve_filt1;

% 

%% TEST
% B = abs(despread_matrix*diag(H_bob_RX)).^2;
% C = abs(despread_matrix*diag(abs(H_eve_RX).^2)).^2;
% A = abs(despread_matrix*diag(W)).^2 ;
% 
% 
% for ii=1:nb_symb
%     cov_termA(:,ii) = nonzeros(A(ii,:));
%     cov_termB(:,ii) = nonzeros(B(ii,:));
%     cov_termC(:,ii) = nonzeros(C(ii,:));
%     covarianceAB = cov(cov_termA(:,ii),cov_termB(:,ii));
%     cov_termAB(ii) = covarianceAB(2);
%     covarianceAB_C = cov(cov_termA(:,ii).*cov_termB(:,ii),cov_termC(:,ii));
%     cov_termAB_C(ii) = covarianceAB_C(2);
% end

% tmp1 = despread_matrix*diag(abs(H_bob_RX).^2);
% tmp2 = despread_matrix*diag(abs(W).^2) ;
% % 
% % 
% % 
% for ii=1:nb_symb
%     cov_term1(:,ii) = nonzeros(tmp1(ii,:));
%     cov_term2(:,ii) = nonzeros(tmp2(ii,:));
%     covariance1 = cov(cov_term1(:,ii),cov_term2(:,ii));
%     cov_term(ii) = covariance1(2);
% end
% 
% covfilt1(cc,nn,nb_bor,aa) = mean(cov_term);

% % covABfilt1(cc,nn,nb_bor,aa) = mean(cov_termAB);% tmp(2);
% % covAB_Cfilt1(cc,nn,nb_bor,aa) = mean(cov_termAB_C);% tmp(2);
% % 
% 
% 
tmp =cov(abs(H_bob_RX(1:BOR(nb_bor):end)).^2,abs(W(1:BOR(nb_bor):end)).^2); % cov(abs(filt1*diag(H_eve_RX.*W)).^2);%;
covfilt1(cc,nn,nb_bor,aa) = tmp(2);

%%
% SINR at one subcarr

SINR_bob_despread_one_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_bob(1,:)).^2*BOR(nb_bor))/ energy_noise_despread_bob; 

SINR_eve_despread_one_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_eve(1,:)).^2*BOR(nb_bor))/ (energy_sym_AN_RX_despread_eve + energy_noise_despread_eve);



energy_signal_bob_one_subcar(cc,nn,nb_bor,aa) =  mean(abs(sym_useful_RX_despread_bob(1,:)).^2*BOR(nb_bor));
energy_noise_bob_one_subcar(cc,nn,nb_bor,aa) =  mean(energy_noise_bob);
energy_noise_eve_one_subcar(cc,nn,nb_bor,aa) =  mean(energy_noise_eve);



% 5.1.2) After equalization

% noise symbols (noise not identical since not the same equalization)
noise_despread_bob_eq =  sym_RX_despread_bob_eq - sym_useful_RX_despread_bob_eq;
noise_despread_eve_eq =  sym_RX_despread_eve_eq - sym_AN_RX_despread_eve_eq - sym_useful_RX_despread_eve_eq;

% mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_bob_eq).^2));
energy_noise_despread_bob_eq = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_bob_eq).^2));
energy_noise_despread_eve_eq = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_eve_eq).^2));

%useful symbol
energy_sym_useful_RX_despread_bob_eq = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_bob_eq).^2));
energy_sym_useful_RX_despread_eve_eq = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_eve_eq).^2));

% energy_sym_useful_RX_despread_bob_eq = mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob_eq).^2));
% energy_sym_useful_RX_despread_eve_eq = mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve_eq).^2));

% AN symbols (only @ Eve)
energy_sym_AN_RX_despread_eve_eq = 1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve_eq).^2); %1/nb_subcar*sum(abs(sym_AN_RX_despread_eve_eq).^2);


% SINR after eq
SINR_bob_despread_eq(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_bob_eq/energy_noise_despread_bob_eq;
SINR_eve_despread_eq(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_eve_eq/(energy_noise_despread_eve_eq + energy_sym_AN_RX_despread_eve_eq);


% Instantaneous capacity
capa_bob_despread(cc,nn,nb_bor,aa) = log2(1+SINR_bob_despread(cc,nn,nb_bor,aa));
capa_eve_despread(cc,nn,nb_bor,aa) = log2(1+SINR_eve_despread(cc,nn,nb_bor,aa));
capa_eve_filt1(cc,nn,nb_bor,aa) = log2(1+SINR_eve_filt1(cc,nn,nb_bor,aa));
capa_eve_filt2(cc,nn,nb_bor,aa) = log2(1+SINR_eve_filt2(cc,nn,nb_bor,aa));
capa_eve_filt3(cc,nn,nb_bor,aa) = log2(1+SINR_eve_filt3(cc,nn,nb_bor,aa));



% Instantaneous secrecy rate:
secrecy_capa_despread(cc,nn,nb_bor,aa) = capa_bob_despread(cc,nn,nb_bor,aa) - capa_eve_despread(cc,nn,nb_bor,aa);
secrecy_capa_filt1(cc,nn,nb_bor,aa) = capa_bob_despread(cc,nn,nb_bor,aa) - capa_eve_filt1(cc,nn,nb_bor,aa);
secrecy_capa_filt2(cc,nn,nb_bor,aa) = capa_bob_despread(cc,nn,nb_bor,aa) - capa_eve_filt2(cc,nn,nb_bor,aa);
secrecy_capa_filt3(cc,nn,nb_bor,aa) = capa_bob_despread(cc,nn,nb_bor,aa) - capa_eve_filt3(cc,nn,nb_bor,aa);






%%  6. De-Parralelize 
sym_RX_despread_bob     = reshape(sym_RX_despread_bob,[],1);
sym_RX_despread_bob_eq  = reshape(sym_RX_despread_bob_eq,[],1)/sqrt(alpha(aa));

sym_RX_despread_eve       = reshape(sym_RX_despread_eve,[],1);
sym_RX_despread_eve_eq    = reshape(sym_RX_despread_eve_eq,[],1)/sqrt(alpha(aa));


%% 7. Demodulation and bit recovery
msg_RX_bob = qamdemod(sym_RX_despread_bob_eq,M,'gray','UnitAveragePower',true,'OutputType','bit');
msg_RX_eve = qamdemod(sym_RX_despread_eve_eq,M,'gray','UnitAveragePower',true,'OutputType','bit');




%% 8. BER

[~,BER_bob(cc,nn,nb_bor,aa)] = biterr(msg_RX_bob,msg_TX); %/nb_symb(nb_bor);
[~,BER_eve(cc,nn,nb_bor,aa)] = biterr(msg_RX_eve,msg_TX); %/nb_symb(nb_bor);

% clc;
% fprintf(progress_bar(nb_bor,length(BOR),'BOR'))
% % fprintf(progress_bar(cc,nb_channels,'\nChannel'))
% fprintf(progress_bar(nn,length(EbN0),'\nNoise'))




end
end
clc;
fprintf(progress_bar(nn,length(EbN0),'Noise level'))
fprintf(progress_bar(aa,length(alpha),'AN'))
fprintf(progress_bar(nb_bor,length(BOR),'Back Of Rate'))
fprintf(progress_bar(cc,nb_channels,'Channel Status'))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%              RESULTS                 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BER_bob = squeeze(mean(BER_bob,1));             % Mean over all channels
BER_eve = squeeze(mean(BER_eve,1));             % Mean over all channels
energy_radiated = squeeze(mean(energy_radiated,1));
 
 


%% 1  Numerator

term_an = anfilt1;

numbob = squeeze(mean(numbob,1));
numfilt0 = squeeze(mean(numfilt0,1)); 
numfilt1 = squeeze(mean(numfilt1,1)); 
numfilt2 = squeeze(mean(numfilt2,1)); 
numfilt3 = squeeze(mean(numfilt3,1)); 

%% 2 Noise
noisefilt1 = squeeze(mean(noisefilt1,1)); 
noisefilt2 = squeeze(mean(noisefilt2,1));
noisefilt3 = squeeze(mean(noisefilt3,1));

%% 3 AN
anfilt1 = squeeze(mean(anfilt1,1));
anfilt3 = squeeze(mean(anfilt3,1));

covfilt1 = squeeze(mean(covfilt1,1));
% covABfilt1 = squeeze(mean(covABfilt1,1));
% covAB_Cfilt1 = squeeze(mean(covAB_Cfilt1,1));
% Denom
denomfilt1 = squeeze(mean(denomfilt1,1));


%% 4. THEORETICAL COMPUTATIONS
energy_an = squeeze(mean(energy_an,1));
energy_signal_bob_one_subcar = squeeze(mean(energy_signal_bob_one_subcar,1));
energy_noise_bob_one_subcar = squeeze(mean(energy_noise_bob_one_subcar,1));

%energy_signal_eve_one_subcar = squeeze(mean(energy_signal_eve_one_subcar,1));
energy_noise_eve_one_subcar = squeeze(mean(energy_noise_eve_one_subcar,1));




% Filt 0

% @ BOB
% figure; % numerator Bob
% plot(numbob','--'); hold on;
% plot((alpha.*(BOR'+1)./BOR.').')
% title('Bob SINR numerator')
% 
% 
% @ EVE
% figure;                 % Num of SINR e classic
% plot(numfilt0','--');
% hold on;
% plot((alpha./BOR')')
% title('Eve SINR numerator term, stupid decoder')
% 

% Filt 1 (@Eve only) 
% figure                 % Num of SINR e macthed filtering
% plot(numfilt1','--');
% hold on;
% plot((alpha.*((BOR.'+3)./BOR.')).')
% title('Eve SINR numerator, filt1')


figure  % AN matched filtering
plot(anfilt1','--'); hold on;
plot((2.*(1-alpha).*(energy_an+covfilt1)).');
title('Eve SINR AN term, matched filter')


% figure
% plot(noisefilt1','--'); hold on;
% plot((energy_noise_eve_one_subcar).')
% title('Eve SINR noise term, matched filter')


%  Filt 2 (@Eve only) 
% figure;
% plot(numfilt2','--'); hold on;
% plot((alpha.*(BOR.'+1)./BOR.').');
% title('Eve SINR numerator term, kill AN')

% figure;
% plot(noisefilt2','--'); hold on;
% plot((energy_noise_eve_one_subcar.*0.6).'); hold on; %90%
% plot((energy_noise_eve_one_subcar).');hold on; %75%
% plot((energy_noise_eve_one_subcar.*2).');hold on; %50%
% plot((energy_noise_eve_one_subcar.*4.8).');hold on; %25%
% plot((energy_noise_eve_one_subcar.*13.15).');hold on; %10%
% plot((energy_noise_eve_one_subcar.*27).');hold on; %5%
% title('Eve SINR noise term, kill AN')


% Filt 3 (@Eve only) 
% figure                 % Num of SINR e macthed filtering
% plot(numfilt3','--');
% title('Eve SINR numerator, LMMSE')
% 
% figure  % AN matched filtering
% plot(anfilt3','--'); 
% title('Eve SINR AN term, LMMSE')
% 
% figure
% plot(noisefilt3','--'); 
% title('Eve SINR noise term, LMMSE')




%% 5 THEORETICAL SINR 
% 5.1 @ Bob:
SINR_bob_theoretical = alpha.*(BOR'+1)./(energy_noise_bob_one_subcar.*BOR');

% 5.2 à Eve
SINR_eve_theoretical_filt0 = mean(alpha./BOR' ./ ((1-alpha).*energy_an  + energy_noise_eve_one_subcar));
SINR_eve_theoretical_filt1 = (alpha.*((BOR.'+3)./BOR.')) ./  (  (BOR.'+1)./BOR.'.*(1-alpha).*(energy_an+covfilt1) +  energy_noise_eve_one_subcar  ) ; 



%% 6 THEORETICAL SECRECY CAPACITY
secrecy_capa_theoretical_filt0 = log2(1+SINR_bob_theoretical)-log2(1+SINR_eve_theoretical_filt0);
secrecy_capa_theoretical_filt1 = log2(1+SINR_bob_theoretical)-log2(1+SINR_eve_theoretical_filt1);





%% 7. ERGODIC SINR

SINR_bob_ergo = squeeze(mean(SINR_bob_despread,1));
SINR_eve_ergo_filt0 = squeeze(mean(SINR_eve_despread,1));
SINR_eve_ergo_filt1 = squeeze(mean(SINR_eve_filt1,1));
SINR_eve_ergo_filt2 = squeeze(mean(SINR_eve_filt2,1));
SINR_eve_ergo_filt3 = squeeze(mean(SINR_eve_filt3,1));
% 
% 
% %Not important
% SINR_bob_ergo_one_subcar = squeeze(mean(SINR_bob_despread_one_subcar,1));
% SINR_eve_ergo_one_subcar = squeeze(mean(SINR_eve_despread_one_subcar,1));




%% 8 ERGODIC CAPACITIES
capa_bob_ergo = squeeze(mean(capa_bob_despread,1));
capa_eve_ergo_filt0 = squeeze(mean(capa_eve_despread,1));
capa_eve_ergo_filt1 = squeeze(mean(capa_eve_filt1,1));
capa_eve_ergo_filt2 = squeeze(mean(capa_eve_filt2,1));
capa_eve_ergo_filt3 = squeeze(mean(capa_eve_filt3,1));


%% 9 ERGODIC SECRECY CAPACITY
secrecy_capa_ergo_filt0 = squeeze(mean(secrecy_capa_despread,1));
secrecy_capa_ergo_filt1 = squeeze(mean(secrecy_capa_filt1,1));
secrecy_capa_ergo_filt2 = squeeze(mean(secrecy_capa_filt2,1));
secrecy_capa_ergo_filt3 = squeeze(mean(secrecy_capa_filt3,1));



secrecy_ergo_max = max(secrecy_capa_ergo_filt0');
secrecy_ergo_filt1_max = max(secrecy_capa_ergo_filt1');
secrecy_ergo_filt2_max = max(secrecy_capa_ergo_filt2');
secrecy_ergo_filt3_max = max(secrecy_capa_ergo_filt3');



%% Plot main results
plot_result;
% save('secrecy_capa.mat', 'secrecy_ergo_filt2_max',  'secrecy_ergo_filt1_max' , 'secrecy_ergo_max')

% figure;
% plot(100*alpha,SINR_eve_ergo_filt1,'--'); hold on ;
% plot(100*alpha,SINR_eve_theoretical_filt1)

% 
% 
% for ii = 1:nb_channels
%     for jj = 1:nb_bor
%         mean_an(ii,jj) = mean(term_an(1:ii,jj));
%     end
% end
% 
% figure;
% plot(mean_an)
% xlabel('Number of iterations')
% ylabel('Moving mean of AN term')
% legend('BOR = 2', 'BOR = 4')
% title([num2str(nb_subcar), 'subcarriers']);

end 