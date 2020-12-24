clear all;
% close all;

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

alpha_step = 5;                           % Percentage between subsequent alpha values
alpha = 1;%0:alpha_step/100:1;         

% Communication parameters
Q = 32;
U = 32;
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 10; % energy per bit over noise psd @Bob - dB
EbN0_e = [10]; % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance



sigma_tau = .5e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                        % Approximation of coherence bandwidth

coef_freq = 2*N/6;

delta_f_n = coef_freq.*delta_f_c;   

b_subcar = delta_f_n./N;
x_axis  = delta_f_n./delta_f_c;
for dd = 1:length(b_subcar)
    [H1(:,:,dd), abs_rho(:,dd), T(:,:,dd)] = corr_frequency( Q , b_subcar(dd) , sigma_tau , nb_run ) ;
end


e_sym_decod1_correl      = zeros(nb_run,length(alpha),length(U));
e_sym_decod1_decorrel      = zeros(nb_run,length(alpha),length(U));

sum_double = 0;                                                             % ∑_i ∑_j h_i h_j
for ii = 0: U-2
    for jj = ii+1 : U-1
        for kk = 1 : 1+ii*N
            sum_double = sum_double + 2*real(T(1+ii*N,kk) * ctranspose(T(1+jj*N,kk)) );
        end
    end
end
test_sum_double = 1/U^2*(U+sum_double);


                                                    
term_exponent_4 = 0;
term_exponent_2 = 0;
for ii = 0: U-2
    for jj = ii+1 : U-1
        for kk = 1 : 1+ii*N
            term_exponent_4 = term_exponent_4 + 2*abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,kk))^2;
        end
    end
end

for ii = 0: U-2
    for jj = ii+1 : U-1
        for kk = 1 : 1+ii*N
            for ll = 1 : 1+jj*N
                if ll ~= kk
                    term_exponent_2 = term_exponent_2 + abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,ll))^2;
                end   
            end
        end
    end
end

sum_double_square = 2*(term_exponent_4 + term_exponent_2);                  % ∑_i ∑_j |h_i|^2 |h_j|^2

sum_fourth = 2*U;

energy_data = alpha./(U^2)*(sum_fourth + sum_double_square);                  % E[ |S^h |Hb|^2 S|^2 ]


tmp3 = 0;

for ii = 0:U-1
    for kk = 1:1+ii*N
        tmp3 = tmp3 + 2*abs(T(1+ii*N,kk))^4;
        for ll = 1:1+ii*N
            if ll ~= kk
                tmp3 = tmp3 + 2*abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^2;
            end
        end
    end
end

tmp4 = 0;
for ii = 0:U-1
    for kk = 1:1+ii*N
        kk;
        %tmp4 = tmp4 + abs(T(1+ii*N,kk))^2;
        for ll = 1:1+ii*N
            ll;
           tmp4 = tmp4 + abs(T(1+ii*N,ll))^2*abs(T(1+ii*N,kk))^2;
        end
    end
end
tmp4 = 2*tmp4;

for iter = 1:nb_run




for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


for dd = 1:length(b_subcar)
%channel generation
H_correl_TX = diag(squeeze(H1(iter,:,dd)).');
H_correl_RX = ctranspose(H_correl_TX);
H_decorrel_TX = channelRayleigh(Q, mu , sigma); %diag(squeeze(H2(iter,:,dd)).');
H_decorrel_RX = ctranspose(H_decorrel_TX);


% check si TEST2 = sum_double_square

tmp = 0;
for ii = 0:U-1
    for jj = 0:U-1
        if jj ~= ii
            tmp = tmp + abs(H_correl_TX(1+ii*N,1+ii*N))^2*abs(H_correl_RX(1+jj*N,1+jj*N))^2;
        end
    end
end
TEST2(iter,dd) = tmp;


% check the value of ∑_i |h_n+iN|^4 = 2*U

tmp2 = 0;
for ii = 0:U-1
    tmp2 = tmp2 + abs(H_correl_TX(1+ii*N,1+ii*N))^4;
end
TEST3(iter,dd) = tmp2;


TEST4(iter,dd) = 1/U^2*(TEST2(iter,dd)+TEST3(iter,dd));
%abs(H_correl_TX(1,1)).^2*abs(H_correl_RX(3,3)).^2 + abs(H_correl_RX(3,3)).^2*abs(H_correl_TX(1,1)).^2;
H2(iter,:,dd) = diag(H_decorrel_TX); 
%% Encoder

% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % Qx1
sym_correl = sym_spread; %H_correl_TX*sym_spread; % Qx1, not weighted 
sym_decorrel = sym_spread;%H_decorrel_TX*sym_spread; % Qx1, not weighted 

sym_precoded_correl = H_correl_TX*sym_correl;
sym_precoded_decorrel = H_decorrel_TX*sym_decorrel;

% AN generation
% an = generateAN(Hb_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd"); % Qx1, not weighted


for aa = 1:length(alpha)
% Weighting
sym_correl_TX = sqrt(alpha(aa))*sym_precoded_correl;     % weighted
sym_decorrel_TX = sqrt(alpha(aa))*sym_precoded_decorrel;     % weighted

% an_TX = sqrt(1-alpha(aa))*an;                       % weighted
% e_an_TX(iter,aa) = energy(an_TX);

%% Receiver
% Whole stream
% sym_RX_b = Hb_RX*(sym_precoded_TX + an_TX);
% sym_RX_e = He_RX*(sym_precoded_TX + an_TX);

% Useful symbol
sym_correl_RX = H_correl_RX*sym_correl_TX; % Qx1
sym_decorrel_RX = H_decorrel_RX*sym_decorrel_TX;




% Noise symbol
% [noise_b, ~ ] = addNoise(sym_b , snr_b, energy(sym_precoded_TX+an_TX));     %addNoise(sym_b , snr_b, energy(sym_RX_b)); 
% [noise_e, e_noise_e ] = addNoise(sym_e , snr_e, energy(sym_precoded_TX+an_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   

% AN symbol
% an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob


%% Decoder
decod1 = matrix_despread;                                   % despreading
% decod2 = matrix_despread*Hb_RX*He_TX;                       % matched filter
% decod3 = matrix_despread*(Hb_RX/He_RX);                     % AN killer
% gamma_E = (He_RX*Hb_TX)*matrix_spread;
% gamma_EH = ctranspose(gamma_E);
% decod4 = sqrt(alpha(aa))*gamma_EH/( alpha(aa)*gamma_E*gamma_EH + (1-alpha(aa))*abs(He_RX).^2*energy(an) + e_noise_e*eye(Q) );   % LMMSE
% decod5 = matrix_despread*He_TX;                             % Only He known by Eve
% 
sym_decod1_correl = decod1*sym_correl_RX;
sym_decod1_decorrel = decod1*sym_decorrel_RX;
% sym_decod2_e = decod2*sym_e;
% sym_decod3_e = decod3*sym_e;
% sym_decod4_e = decod4*sym_e;
% sym_decod5_e = decod5*sym_e;
% 
% noise_decod1_b = decod1*noise_b;
% noise_decod1_e = decod1*noise_e;
% noise_decod2_e = decod2*noise_e;
% noise_decod3_e = decod3*noise_e;
% noise_decod4_e = decod4*noise_e;
% noise_decod5_e = decod5*noise_e;
% 
% an_decod1_e = decod1*an_e;
% an_decod2_e = decod2*an_e;
% an_decod3_e = decod3*an_e;
% an_decod4_e = decod4*an_e;
% an_decod5_e = decod5*an_e;

%% Energy of the different RX components 
% @ Bob
e_sym_decod1_correl(iter,aa,bb,dd) = energy(sym_decod1_correl);%  energy(matrix_despread*abs(H_correl_RX).^2*matrix_spread);
e_sym_decod1_decorrel(iter,aa,bb,dd) = energy(sym_decod1_decorrel);
% 
% % @ Eve : decod 1
% e_sym_decod1_e(iter,aa,bb,dd)     = energy(sym_decod1_e);
% e_noise_decod1_e(iter,aa,bb,dd)   = energy(noise_decod1_e);
% e_an_decod1_e(iter,aa,bb,dd)      = energy(an_decod1_e);
% e_denom_decod1_e(iter,aa,bb,dd)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve
% 
% % @ Eve : decod 2
% e_sym_decod2_e(iter,aa,bb,dd)     = energy(sym_decod2_e);
% e_noise_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e);
% e_an_decod2_e(iter,aa,bb,dd)      = energy(an_decod2_e);
% e_denom_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve
% 
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


e_avg_sym_decod1_correl      = squeeze(mean(e_sym_decod1_correl,1)); 
e_avg_sym_decod1_decorrel      = squeeze(mean(e_sym_decod1_decorrel,1));  















