clear all;
close all; 
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',32)
set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'DefaultLineMarkerSize',15);
set(0, 'defaultFigurePosition',  [-1267  44   1256    872])

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
for nb_bor = 1:length(BOR)
msg_TX      = randi( [0 1] , nb_bit(nb_bor) , 1 ) ;                                 % Random bit data stream
sym_TX      = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types
sym_TX      = reshape(sym_TX, [], nb_block) ;                                       % OFDM block per OFDM block


%% ICI: Creation du preambule qu'il faut aussi spread pour estimer le canal de Eve et Bob pour l'�galisation

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


for cc = 1:nb_channels
%% 4. Channel generation : 
% Extrended Pedestrian Model (epa), normalized or statistical

% 4.1) Legitimate receiver (bob)

%h_bob  = TR_security_generate_rayleigh_channel(delta_t,42,nb_subcar,"epa");
%h_bob = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"normalized");
%h_bob = TR_security_generate_rayleigh_channel(nb_taps,nb_taps,nb_subcar,"statistical");
H_bob = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"decorrelated");
%H_bob = corr_rayleigh(nb_subcar,64,.5,1);
% 4.2) Eavesdropper (Eve)

%h_eve  = TR_security_generate_rayleigh_channel(delta_t,42,nb_subcar,"epa");
%h_eve = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"normalized");
%h_eve = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"statistical");
H_eve = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"decorrelated");
%H_eve = corr_rayleigh(nb_subcar,64,.5,1);

% Frequency domain:

% H_bob           = fft(h_bob,nb_subcar); %.*sqrt(1/nb_taps);
norm_H_bob_tmp  = abs(H_bob).^2;
% energy_H_bob    = sum(abs(H_bob).^2);
% H_bob           = H_bob/sqrt(energy_H_bob);


% H_eve           = fft(h_eve,nb_subcar); %.*sqrt(1/nb_taps);
norm_H_eve_tmp  = H_bob .* H_eve;
% energy_H_eve    = sum(abs(H_eve).^2);
% H_eve           = H_eve/sqrt(energy_H_eve);


norm_H_bob = sum(norm_H_bob_tmp,1);                           % Contributrion of all antennas (not use for now)
norm_H_eve = sum(norm_H_eve_tmp,1);
H_bob_TX    = ctranspose(H_bob).';
H_bob_RX    = H_bob;
H_eve_RX    = H_eve;




%% 5.  ARTIFICIAL NOISE INJECTION

W = TR_security_generate_AN(H_bob,nb_subcar,BOR(nb_bor),despread_matrix,"underdetermined"); % Resolution of the underdetermined system
% W = TR_security_generate_AN(H_bob,nb_subcar,BOR(nb_bor),despread_matrix,"determined");     % Resolution of the determined system
w_signal = ifft(W);
PAPR_AN(cc,nb_bor) = PAPR(w_signal);                


% Normalization of W
tmp = sum(abs(W).^2)/nb_subcar;
energy_sym_useful_TX = mean(sum(abs(H_bob_TX.*sym_spread).^2))/nb_subcar;
W = W/sqrt(tmp)*sqrt(energy_sym_useful_TX);
energy_W = sum(abs(W).^2)/nb_subcar;%sum(abs(W).^2)/nb_subcar;

%% 6. Symbol stream generation
for aa = 1:length(alpha)
    
sym_useful =sqrt(alpha(aa))*H_bob_TX.*sym_spread;       % Energy --> sqrt needed
AN = sqrt(1 - alpha(aa))*W;
energy_AN = sum(abs(AN).^2)/nb_subcar;
% energy_qtt = mean(sum(abs(sym_useful).^2)/nb_subcar);

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
matched_filt = diag(H_bob_RX)*diag(ctranspose(H_eve_RX)); 
sym_RX_eve_matched      = matched_filt*sym_RX_eve;
sym_useful_eve_matched  = matched_filt*sym_useful_RX_eve;
sym_AN_RX_eve_matched   = matched_filt*sym_AN_RX_eve;
sym_RX_noisy_eve_matched = matched_filt*sym_RX_noisy_eve;

%% 3. Despreading

% 3.1) Useful symbols
sym_useful_RX_despread_bob = despread_matrix*sym_useful_RX_bob; % Despread useful symbols @ Bob
sym_useful_RX_despread_eve = despread_matrix*sym_useful_RX_eve; % Despread useful symbols @ Eve
sym_useful_RX_despread_eve_matched = despread_matrix*sym_useful_eve_matched;

% 3.2) AN symbols (only @ Eve)
sym_AN_RX_despread_eve = despread_matrix*sym_AN_RX_eve; % Despread AN symbols @ Eve
sym_AN_RX_despread_eve_matched = despread_matrix*sym_AN_RX_eve_matched;

% 3.3) useful symbols + AN
sym_RX_despread_bob = despread_matrix*sym_RX_noisy_bob; % Depsread symbols (useful + AN) @ Bob
sym_RX_despread_eve = despread_matrix*sym_RX_noisy_eve; % Depsread symbols (useful + AN) @ Eve
sym_RX_despread_eve_matched = despread_matrix*sym_RX_noisy_eve_matched;




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
noise_despread_eve_matched = sym_RX_despread_eve_matched - sym_AN_RX_despread_eve_matched - sym_useful_RX_despread_eve_matched;

energy_noise_despread_bob = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_bob).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_bob).^2));
energy_noise_despread_eve = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_eve).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_eve).^2));
energy_noise_despread_eve_matched = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_eve_matched).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_eve).^2));


% useful symbol
energy_sym_useful_RX_despread_bob = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_bob).^2)); %mean(var(sym_useful_RX_despread_bob)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2)); %mean(var(sym_useful_RX_despread_bob)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2));
energy_sym_useful_RX_despread_eve = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_eve).^2)); %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));
energy_sym_useful_RX_despread_eve_matched = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_eve_matched).^2)); %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));

% AN symbols (only @ Eve)
energy_sym_AN_RX_despread_eve =  mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve).^2)); %var(sym_AN_RX_despread_eve) ; % 1/nb_subcar*sum(abs(sym_AN_RX_despread_eve).^2); % %
energy_an(cc,nn,nb_bor,aa) = energy_W;%energy_sym_AN_RX_despread_eve;


energy_denom_eve =  mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve + noise_despread_eve).^2));
energy_denom_eve_matched = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve_matched + noise_despread_eve_matched).^2));

% SINR after despread
SINR_bob_despread(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_bob/energy_noise_despread_bob;
SINR_eve_despread(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_eve/energy_denom_eve;
SINR_eve_despread_matched(cc,nn,nb_bor,aa) = energy_sym_useful_RX_despread_eve_matched/energy_denom_eve_matched;


% SINR at one subcarr

SINR_bob_despread_one_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_bob(1,:)).^2*BOR(nb_bor))/ energy_noise_despread_bob; 

SINR_eve_despread_one_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_eve(1,:)).^2*BOR(nb_bor))/ (energy_sym_AN_RX_despread_eve + energy_noise_despread_eve);


% SINR_bob_despread_one_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_bob(1,:)).^2)/ mean(noise_power_bob); %mean(abs(noise_despread_bob(1,:)).^2);
% SINR_bob_despread_second_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_bob(2,:)).^2)/ mean(noise_power_bob); %mean(abs(noise_despread_bob(2,:)).^2);
% 
% SINR_eve_despread_one_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_eve(1,:)).^2)/ (mean(noise_power_eve)  + energy_sym_AN_RX_despread_eve); %mean(abs(noise_despread_bob(1,:)).^2);
% SINR_eve_despread_second_subcar(cc,nn,nb_bor,aa) = mean(abs(sym_useful_RX_despread_eve(2,:)).^2)/ (mean(noise_power_eve)  + energy_sym_AN_RX_despread_eve); %mean(abs(noise_despread_bob(2,:)).^2);



energy_signal_bob_one_subcar(cc,nn,nb_bor,aa) =  mean(abs(sym_useful_RX_despread_bob(1,:)).^2*BOR(nb_bor));
energy_noise_bob_one_subcar(cc,nn,nb_bor,aa) =  mean(energy_noise_bob);
energy_noise_eve_one_subcar(cc,nn,nb_bor,aa) =  mean(energy_noise_eve);

% signal_power_one_subcar_bob(cc,nn,nb_bor,aa) =  mean(abs(sym_useful_RX_despread_bob(1,:)).^2);
% noise_power_one_subcar_bob(cc,nn,nb_bor,aa) =  mean(noise_power_bob);
% noise_power_one_subcar_eve(cc,nn,nb_bor,aa) =  mean(noise_power_eve);

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
capa_eve_despread_matched(cc,nn,nb_bor,aa) = log2(1+SINR_eve_despread_matched(cc,nn,nb_bor,aa));

% Instantaneous secrecy rate:
secrecy_capa_despread(cc,nn,nb_bor,aa) = capa_bob_despread(cc,nn,nb_bor,aa) - capa_eve_despread(cc,nn,nb_bor,aa);
secrecy_capa_despread_matched(cc,nn,nb_bor,aa) = capa_bob_despread(cc,nn,nb_bor,aa) - capa_eve_despread_matched(cc,nn,nb_bor,aa);





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
% fprintf(progress_bar(cc,nb_channels,'\nChannel'))
% fprintf(progress_bar(nn,length(EbN0),'\nNoise'))




end
end
clc;
fprintf(syyydlib.progress_bar(nn,length(EbN0),'Noise level'))
fprintf(syyydlib.progress_bar(aa,length(alpha),'AN'))
fprintf(syyydlib.progress_bar(nb_bor,length(BOR),'Back Of Rate'))
fprintf(syyydlib.progress_bar(cc,nb_channels,'Channel Status'))
end
end


%% PART 4. RESULTS

%% 1. BIT ERROR RATE
BER_bob = squeeze(mean(BER_bob,1));             % Mean over all channels
BER_eve = squeeze(mean(BER_eve,1));             % Mean over all channels


%% 2. RADIATED ENERGY
energy_radiated = squeeze(mean(energy_radiated,1));
 
 
%% 3. ERGODIC SINR

% 3.1 Ergodic SINR After despread
SINR_bob_ergo = squeeze(mean(SINR_bob_despread,1));
SINR_eve_ergo = squeeze(mean(SINR_eve_despread,1));
SINR_eve_ergo_matched = squeeze(mean(SINR_eve_despread_matched,1));

%Not important
SINR_bob_ergo_one_subcar = squeeze(mean(SINR_bob_despread_one_subcar,1));
SINR_eve_ergo_one_subcar = squeeze(mean(SINR_eve_despread_one_subcar,1));



%% 4 ERGODIC CAPACITIES
capa_bob_ergo = squeeze(mean(capa_bob_despread,1));
capa_eve_ergo = squeeze(mean(capa_eve_despread,1));
capa_eve_ergo_matched = squeeze(mean(capa_eve_despread_matched,1));

%% 5 ERGODIC SECRECY CAPACITY
secrecy_capa_ergo = squeeze(mean(secrecy_capa_despread,1));
secrecy_capa_ergo_matched = squeeze(mean(secrecy_capa_despread_matched,1));

secrecy_ergo_max = max(secrecy_capa_ergo');
secrecy_ergo_matched_max = max(secrecy_capa_ergo_matched');

plot_result;
save('secrecy_capa.mat', 'secrecy_ergo_matched_max' , 'secrecy_ergo_max')

end 