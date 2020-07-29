clear all;
%close all; 
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
%   SR computed after despreading
%
%   Frequency domain secure TR in a MISO configuration    
%   Complex Envelope (CE) simulations
%
%   OFDM to mitigate MPC effects
%
%   Artificial Noise added to secure the data transmission at Bob
%
%   Simulation of the effect of spatial correlation between Eve and Bob on
%   the SR. 
%   
%	Example code to call the function corr_spatial to generates 
%   nb_subca independent Rayleigh channel exhibiting spatial correlation 
%   assuming a three-dimensional Rayleigh channel
%
%   The spatial dependence of the correlation is given by eq. 24
%   from Ph. De Doncker, "Spatial Correlation Functions for Fields in Three
%   -Dimensional Rayleigh Channels," PIER, Vol. 40, 55-69, 2003
%
%   Inputs:
%   nb_subca: number of sub-carrier (or frequency samples within the
%   frequency bandwidth)
%   b_subca: bandwidth of each subcarrier (of frequency separation between 
%   frequency samples)
%   r: vector of the considered distances away from Bob
%   fc: carrier frequency
%   nb_realizations: number of samples of each sub-channel 
%
%   Outputs:
%   H: channel matrix whose size is (nb_realizations x nb_subca x nb_distance)

%
%
%   Code started: 23 September 2019
%   Last update: 23 September 2019
%
%
%   by Sidney Golstein & Julien Sarrazin
%
%**************************************************************************





tic

sc = "secrecy_alpha_bor_ebno";
%sc = "secrecy_alpha_ebno_bor";


close all;
clearvars -except sc scenar
load('alpha.mat', 'alpha_global')
load('secrecy_max.mat','max_value','max_value_theoretical')
scenario = sc;
get_parameters

nb_channels = 500;
% Number of realizations
nb_pos = 80;                                           % Number of spatial points between Bob and Eve
b_subcar = .1e6 ;                                           % Bandwidth of each subcarrier
r = linspace( 0 , 2*c/fc , nb_pos ) ;                    % Eve's distance from Bob . Distance from 0 to 1 with 40 points

dist = r./(c/fc);
%% Random Variable Calculation
[H , RHO_SPATIAL] = corr_spatial( nb_subcar , b_subcar , r , fc , nb_channels ) ;







%% TRANSMISSION PART
for nb_bor = 1:length(BOR)
for cc = 1:nb_channels
msg_TX      = randi( [0 1] , nb_bit(nb_bor) , 1 ) ;                                 % Random bit data stream
sym_TX      = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types
sym_TX      = reshape(sym_TX, [], nb_block) ;                                       % OFDM block per OFDM block
                                      


                                      
spread_matrix   = zeros( nb_subcar , nb_symb_per_block(nb_bor) );
spread_code     = reshape( 2.*randi([0 1], 1, nb_subcar).' - 1 , [] , BOR(nb_bor) )/sqrt(BOR(nb_bor));   % Diagonal elements of the spreading matrix. BPSK spreading code


for i = 1:size(spread_code,2)
    spread_matrix(i*nb_symb_per_block(nb_bor)-nb_symb_per_block(nb_bor)+1:i*end/BOR(nb_bor),:) = diag(spread_code(:,i));
end
despread_matrix = ctranspose(spread_matrix); 

%% 3. Symbol stream after spreading
sym_spread = spread_matrix*sym_TX;





for dd = 1:nb_pos
H_bob = squeeze(H(cc,:,1)).';
H_eve = squeeze(H(cc,:,dd)).';
H_bob_TX    = ctranspose(H_bob).';
H_bob_RX    = H_bob;
H_eve_RX    = H_eve;



sym_useful_TX = H_bob_TX.*sym_spread;      % ATTENTION sqrt(alpha) needed  
W = TR_security_generate_AN(H_bob,nb_subcar,BOR(nb_bor),despread_matrix, "underdetermined"); % Resolution of the underdetermined system


tmp = sum(abs(W).^2)/nb_subcar;
energy_sym_useful_TX = mean(sum(abs(H_bob_TX.*sym_spread).^2))/nb_subcar;
W = W/sqrt(tmp)*sqrt(energy_sym_useful_TX);




    
    
alpha_to_opt = alpha_global(nb_bor)*ones(nb_subcar,1);                      % Coef alpha to optimize

sym_raditated_noiseless = sqrt(alpha_to_opt).*sym_useful_TX + sqrt(ones(nb_subcar,1) - alpha_to_opt).*W;

%% RECEPTION BEFORE DESPREADING

% Total received symbol at B and E
sym_RX_bob = H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX + sqrt(ones(nb_subcar,1) - alpha_to_opt).*W);
sym_RX_eve = H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX + sqrt(ones(nb_subcar,1) - alpha_to_opt).*W);


% Noise addition on useful symbols + AN
[sym_RX_noisy_bob,energy_noise_bob] = TR_security_addNoise(sym_RX_bob, EbN0, M, var(sym_raditated_noiseless));   % Noisy received symbols (useful + AN) @ Bob
[sym_RX_noisy_eve,energy_noise_eve] = TR_security_addNoise(sym_RX_eve, EbN0, M, var(sym_raditated_noiseless));   % Noisy received symbols (useful + AN) @ Eve



% useful part of received symbol at B and E
sym_useful_RX_bob = H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX);
sym_useful_RX_eve = H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX);

% Artificial noise received symbol at E
sym_AN_RX_eve = H_eve_RX.*(sqrt(ones(nb_subcar,1) - alpha_to_opt).*W);  





%% RECEPTION AFTER DESPREADING

% Useful symbols
sym_useful_RX_despread_bob = despread_matrix*(H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)); % Despread useful symbols @ Bob
sym_useful_RX_despread_eve = despread_matrix*(H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)); % Despread useful symbols @ Eve

% AN symbols (only @ Eve)
sym_AN_RX_despread_eve = despread_matrix*(H_eve_RX.*(sqrt(ones(nb_subcar,1) - alpha_to_opt).*W)); % Despread AN symbols @ Eve

% useful symbols + AN
sym_RX_despread_bob = despread_matrix*sym_RX_noisy_bob; % Depsread symbols (useful + AN) @ Bob
sym_RX_despread_eve = despread_matrix*sym_RX_noisy_eve; % Depsread symbols (useful + AN) @ Eve

% noisy symbols
noise_despread_bob =  sym_RX_despread_bob - sym_useful_RX_despread_bob;
noise_despread_eve =  sym_RX_despread_eve - sym_AN_RX_despread_eve - sym_useful_RX_despread_eve;

%% Energy computation
energy_sym_useful_RX_despread_bob = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_bob).^2)); %mean(var(sym_useful_RX_despread_bob)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2)); %mean(var(sym_useful_RX_despread_bob)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2));
energy_sym_useful_RX_despread_eve = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_useful_RX_despread_eve).^2)); %mean(var(sym_useful_RX_despread_eve)) ; %mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2));
energy_noise_despread_bob = mean(1/nb_subcar*BOR(nb_bor)*sum(abs(noise_despread_bob).^2)); %mean(1/nb_subcar*sum(abs(noise_despread_bob).^2));
energy_denom_eve =  mean(1/nb_subcar*BOR(nb_bor)*sum(abs(sym_AN_RX_despread_eve + noise_despread_eve).^2));


%% SINR COMPUTATION
SINR_bob(cc,nb_bor,dd) = energy_sym_useful_RX_despread_bob/energy_noise_despread_bob;%mean(noise_power_bob); %energy_noise_despread_bob;
SINR_eve(cc,nb_bor,dd) = energy_sym_useful_RX_despread_eve/energy_denom_eve;

%% CAPA COMPUTATION

capa_bob(cc,nb_bor,dd) = log2(1+SINR_bob(cc,nb_bor,dd));
capa_eve(cc,nb_bor,dd) = log2(1+SINR_eve(cc,nb_bor,dd));


secrecy_capa(cc,nb_bor,dd) = capa_bob(cc,nb_bor,dd) - capa_eve(cc,nb_bor,dd) ; 
% capa_bob_per_subcar(:,cc,nb_bor) = log2(1 +  1/nb_subcar*BOR(nb_bor)*sum(abs(despread_matrix*(H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)))).^2 ./ sum(abs(noise_despread_bob).^2));    % Mean of capacity over all data blocks to have the mean capa per subcar
% capa_eve_per_subcar(:,cc,nb_bor) = log2(1 +  1/nb_subcar*BOR(nb_bor)*sum(abs(despread_matrix*(H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)))).^2 ./ sum(abs(noise_despread_eve + despread_matrix*(H_eve_RX.*(sqrt(ones(nb_subcar,1) - alpha_to_opt).*W))).^2));
% 
% secrecy_capa_per_subcar(:,cc,nb_bor) = capa_bob_per_subcar(:,cc,nb_bor) - capa_eve_per_subcar(:,cc,nb_bor);

%
end
clc;
fprintf(syyydlib.progress_bar(cc,nb_channels,'Channel Status'))
fprintf(syyydlib.progress_bar(nb_bor,length(BOR),'Back Of Rate'))
%fprintf(syyydlib.progress_bar(dd,nb_pos,'Eve Position'))


end
end


%% Mean over all blocks and all channels:

capa_bob_ergo = squeeze(mean(capa_bob,1));  % Capa averaged over all channel realization in function of BOR value and 
capa_eve_ergo = squeeze(mean(capa_eve,1));
secrecy_capa_ergo = capa_bob_ergo - capa_eve_ergo;


figure;
yyaxis left
plot(dist,secrecy_capa_ergo')
ylabel('Secrecy Capacity (bit/s/Hz)')
legendCell_secrecy = cellstr(num2str(BOR', 'BOR = %-d')); hold on;

yyaxis right
plot( dist , RHO_SPATIAL( 1 , : , 1))
ylabel('Spatial Correlation')

xlabel('$\lambda$')


legendCell = [legendCell_secrecy];
legend(legendCell)
