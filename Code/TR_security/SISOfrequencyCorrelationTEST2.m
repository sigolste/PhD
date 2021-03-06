%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   Implementation of new approx for Bob capa and Eve capa with SDS/MF/OC 
%   decoders when  frequency corrleation among Bob's subcarriers is  
%   introduced. Modelization of SR thanks to first and second approximation
%   of the ergodic SR
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
alpha = 1;%0:alpha_step/100:1;         

% Communication parameters
Q = 16;
U = 2;
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters

snr_b  = 15;  % SNR @Bob
snr_e  = 15;  % SNR @Eve

% Noise energy
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve


%% Channel modeling

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance
sigma_tau_b = .5e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
sigma_tau_e = .2e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)

delta_f_c_b = 1 / 2 / pi / sigma_tau_b ;                        % Approximation of coherence bandwidth
delta_f_c_e = 1 / 2 / pi / sigma_tau_e ;                        % Approximation of coherence bandwidth

coef_freq_b = [1].*N/6;
coef_freq_e = [1].*N/6;

delta_f_n_b = coef_freq_b.*delta_f_c_b; 
delta_f_n_e = coef_freq_e.*delta_f_c_e;   

b_subcar_b = delta_f_n_b./N;                                    % Bandwidth of each subcarrier
b_subcar_e = delta_f_n_e./N;                                    % Bandwidth of each subcarrier

x_axis_b  = delta_f_n_b./delta_f_c_b;
x_axis_e  = delta_f_n_e./delta_f_c_e;

for dd = 1:length(b_subcar_b)
    [H(:,:,dd), H1(:,:,dd), abs_rho(:,dd), Tb(:,:,dd)] = corr_frequency( Q , b_subcar_b(dd) , sigma_tau_b , nb_run ) ;
end
for dd = 1:length(b_subcar_e)
    [I(:,:,dd), I1(:,:,dd), ~, Te(:,:,dd)] = corr_frequency( Q , b_subcar_e(dd) , sigma_tau_e , nb_run ) ;
end

%% Energy matrix instantiation
e_noise_e            = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_noise_b            = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_an_TX              = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));

e_sym_decod1_b       = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_noise_decod1_b     = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));

e_sym_decod1_e      = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_noise_decod1_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_an_decod1_e       = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_denom_decod1_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));


e_sym_decod2_e      = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_noise_decod2_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_an_decod2_e       = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));
e_denom_decod2_e    = zeros(nb_run,length(alpha),length(U),length(b_subcar_b));

for iter = 1:nb_run

for bb =1:length(U)
    
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


for dd = 1:length(b_subcar_b)
%channel generation
Hb_TX = diag(squeeze(H1(iter,:,dd)).');
Hw_TX = diag(squeeze(H(iter,:,dd)).');

Hb_RX = ctranspose(Hb_TX);
Hw_RX = ctranspose(Hw_TX);


H_test_TX = diag(squeeze(I1(iter,:,dd)).');
Hw_test_TX = diag(squeeze(I(iter,:,dd)).');

H_test_RX = ctranspose(H_test_TX);
Hw_test_RX = ctranspose(Hw_test_TX);


He_TX = diag(squeeze(I1(iter,:,dd)).'); 
%He_TX = channelRayleigh(Q, mu , sigma);             % Assumption: uncorrelated subcarrier for Eve channel
He_RX = ctranspose(He_TX);



% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                                                           %
% %                                                                           %
% %                                                                           %
% %                                                                           %
% %                                                                           %    
% %%                         TEST SECTION                                     %   
% %                                                                           %
% %                                                                           %
% %                                                                           %
% %                                                                           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%              ∑_i |hb|^8
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmp = 0;
% for ii = 0:U-1
%     tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^8;
% end
% h8h8(iter,dd) = tmp;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%          ∑_i ∑_j~=i |hb_i|^2 |Hb_j|^6
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmp = 0;
% for ii = 0:U-1
%     for jj = 0:U-1
%         if ii~= jj
%             tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_TX(1+jj*N,1+jj*N))^6;
%         end
%     end
% end
% h2h6(iter,dd) = tmp;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%          ∑_i ∑_j~=i |hb_i|^4 |Hb_j|^4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = 0:U-1
%         if ii~= jj
%             tmp = tmp + Hb_TX(1+ii*N,1+ii*N)*Hb_RX(1+jj*N,1+jj*N)*H_test_RX(1+ii*N,1+ii*N)*H_test_TX(1+jj*N,1+jj*N);
%         end
%     end
% end
% hh(iter,dd) = tmp;

tmp = 0;
for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj
            tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N)).^2*abs(Hb_TX(1+jj*N,1+jj*N)).^2* abs(He_TX(1+ii*N,1+ii*N)).^2*abs(He_TX(1+jj*N,1+jj*N)).^2;
        end
    end
end
h2h2(iter,dd) = tmp;


tmp = 0;
for ii = 0:U-1
    for jj =0:U-1
        if ii ~= jj
            for ll = 1:1+ii*N
                tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N)).^2*abs(Hb_TX(1+jj*N,1+jj*N)).^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2*abs(Hw_test_TX(ll,ll))^4;
            end
        end 
    end
end
h4e(iter,dd) = tmp;


tmp = 0;
for ii = 0:U-1
    for jj =0:U-1
        if ii ~= jj
            for ll = 1:1+ii*N
                for mm = 1:1+jj*N
                    if mm ~= ll
                        tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N)).^2*abs(Hb_TX(1+jj*N,1+jj*N)).^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2*abs(Hw_test_TX(ll,ll))^2*abs(Hw_test_TX(mm,mm))^2;
                    end
                end
            end
        end 
    end
end
h22e(iter,dd) = tmp;




tmp = 0;
for ii = 0:U-1
    for jj =0:U-1
        if ii ~= jj
            for kk = 1:1+ii*N
                for ll = kk+1:1+ii*N
                   tmp = tmp + 2*abs(Hb_TX(1+ii*N,1+ii*N)).^2*abs(Hb_TX(1+jj*N,1+jj*N)).^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                end
            end
        end 
    end
end
h22e_bis(iter,dd) = tmp;



% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ∑_i ∑_j~=i ∑_k ~= i,j |=hb_i|^2 |Hb_j|^2 |Hb_k|^4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (ii~= jj) && (ii~= kk) && (kk ~= jj)
%                 tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_TX(1+jj*N,1+jj*N))^2*abs(Hb_TX(1+kk*N,1+kk*N))^4;
%             end
%         end
%         
%     end
% end
% h4h2h2(iter,dd) = 2*tmp;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ∑_i ∑_j~=i ∑_k ~= i,j |=hb_i|^2 |Hb_j|^2 |Hb_k|^4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_TX(1+jj*N,1+jj*N))^2*abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%             end
%         end
%         
%     end
% end
% h2h2h2h2(iter,dd) = 24*tmp;
% 
% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ∑_i |=hb_i|^4 |He_i|^4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%
% tmp = 0;
% for ii = 0:U-1
%     tmp = tmp + abs(He_RX(1+ii*N,1+ii*N))^4*abs(Hb_RX(1+ii*N,1+ii*N))^4;
% end
% h4h4(iter,dd) = tmp;
% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = 0:U-1
%         if ii ~= jj
%             tmp = tmp + abs(He_RX(1+ii*N,1+ii*N))^2*abs(He_RX(1+jj*N,1+jj*N))^2* ... 
%                         abs(Hb_RX(1+ii*N,1+ii*N))^2*abs(Hb_RX(1+jj*N,1+jj*N))^2;
%         end
%     end
% end
% h2h2(iter,dd) = tmp;
% 
% %%%
% tmp = 0;
% tmp1 = 0;
% tmp2 = 0;
% for ii = 0:U-1
%     tmp1 = tmp + abs(He_RX(1+ii*N,1+ii*N))^4*abs(Hb_RX(1+ii*N,1+ii*N))^4;
%     for jj = 0:U-1
%         if jj ~= ii 
%                 tmp2 = tmp2 + 2*abs(He_RX(1+ii*N,1+ii*N))^2*abs(He_RX(1+jj*N,1+jj*N))^2* ... 
%                         abs(Hb_RX(1+ii*N,1+ii*N))^2*abs(Hb_RX(1+jj*N,1+jj*N))^2;
%         end
%     end
% end
% abs4_test(iter,dd ) = 1./U.^4 * (tmp1+tmp2);
% 
% 
% for ii = 0:U-1
%     tmp = tmp + He_RX(1+ii*N,1+ii*N).*Hb_TX(1+ii*N,1+ii*N);
% end
% abs4(iter,dd)  = abs(tmp)^4./U.^4;
% 
% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
sym_decod2_e = decod2*sym_e;

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
e_sym_decod2_e(iter,aa,bb,dd)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e);
e_an_decod2_e(iter,aa,bb,dd)      = energy(an_decod2_e);
e_denom_decod2_e(iter,aa,bb,dd)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

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

e_avg_sym_decod2_e      = squeeze(mean(e_sym_decod2_e,1));
e_avg_noise_decod2_e    = squeeze(mean(e_noise_decod2_e,1));
e_avg_an_decod2_e       = squeeze(mean(e_an_decod2_e,1));
e_avg_denom_decod2_e    = squeeze(mean(e_denom_decod2_e,1));


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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %    
%%                         MODELIZATION SECTION                             %   
%                           E[X^2] at Bob                                   %
%                                                                           %
%                                                                           %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % %% Correlation introduced 
% % 
% % % 1. ∑_i |hb_{n+iN}|^8 
% % % H8H8_simu = mean(h8h8)
% H8_model = modelCorrelH8(U);
% % 
% % % 2. ∑_i ∑_j~=i |hb_{n+iN}|^2 * |hb_{n+jN}|^6 
% % % H2H6_simu = mean(h2h6)
% H2H6_model = modelCorrelH2H6(N,U,Tb);
% % 
% % % 3. ∑_i ∑_j~=i |hb_{n+iN}|^4 * |hb_{n+jN}|^4 
% % % H4H4_simu = mean(h4h4)
% H4H4_model  = modelCorrelH4H4(N,U,Tb);
% % 
% % % 4. ∑_i ∑_j~=i ∑_k~=i,j |hb_{n+iN}|^4 * |hb_{n+jN}|^2 * |hb_{n+kN}|^2
% % % H4H2H2_simu = mean(h4h2h2)
% H4H2H2_model = modelCorrelH4H2H2(N,U,Tb);
% % 
% % % 5. ∑_i ∑_j~=i ∑_k~=i,j ∑_l~=i,j,k |hb_{n+iN}|^2 * |hb_{n+jN}|^2 * |hb_{n+kN}|^2 * |hb_{n+lN}|^2
% % % H2H2H2H2_simu = mean(h2h2h2h2)
% H2H2H2H2_model = modelCorrelH2H2H2H2(N,U,Tb);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %    
%%                         MODELIZATION SECTION                             %   
%                           E[X^2] at Eve                                   %
%                                                                           %
%                                                                           %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HH = mean(hh);
% tmp1 = 0;
% tmp2 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         if ii~= jj
%             for kk = 1:1+ii*N
%                 tmp1 = tmp1 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Te(1+ii*N,kk))*Te(1+jj*N,kk));
%                 for ll = 1:1+ii*N
%                     if ll ~= kk
%                         tmp2 = tmp2 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk)))*real(conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
%                     end
%                 end
%                
%             end
%         end
%     end
% end
% HH_test = tmp1+tmp2;

H2H2 = mean(h2h2)

tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
tmp8 = 0;


 for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj
        for kk = 1:1+ii*N
            tmp1 = tmp1 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,kk))^2;
            for ll = 1:1+ii*N
                if ll ~= kk
                    tmp2 = tmp2 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2;
                    tmp3 = tmp3 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2;
                end
            end
            for ll = 1:1+jj*N
                if ll ~= kk
                    tmp4 = tmp4 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,kk))^2;
                    for mm = 1:1+ii*N
                        if (mm ~= ll) && (mm ~= kk)
                            tmp8 = tmp8 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,mm))^2;
                        end
                    end
                end
            end
            for ll = kk+1:1+ii*N
                tmp5 = tmp5 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,kk))^2;
                tmp6 = tmp6 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2;
                for mm = 1:1+ii*N
                    if (mm ~= ll) && (mm ~= kk)
                        tmp7 = tmp7 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,mm))^2;
                    end
                end
            end
        end
        end
    end
end 




H4E = mean(h4e)
H4E_test = tmp1 + tmp2 + tmp3 + tmp4 +  tmp5 + tmp6 + tmp7 + tmp8


%%%%%%


tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
tmp8 = 0;
tmp9 = 0;
tmp10 = 0;
tmp11 = 0;
tmp12 = 0;
tmp13 = 0;
tmp14 = 0;
tmp15 = 0;
tmp16 = 0;
tmp17 = 0;
tmp18 = 0;
tmp19 = 0;
 for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj
            for kk = 1:1+ii*N
                for ll = 1:1+jj*N
                    if ll ~= kk
                        tmp1 = tmp1 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                        tmp2 = tmp2 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                        for mm = 1:1+ii*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp3 = tmp3 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,kk))^2;
                                tmp4 = tmp4 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,ll))^2;
                                for nn = 1:1+jj*N
                                    if (nn ~= mm) && (nn ~= ll) && (nn ~= kk)
                                        tmp5 = tmp5 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,nn))^2;
                                    end
                                end
                            end
                        end
                        for mm = 1:1+jj*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp6 = tmp6 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,mm))^2;
                            end
                        end
                    end
                end
                for ll = 1:1+ii*N
                    if ll ~= kk
                        tmp7 = tmp7 + 2*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                        tmp8 = tmp8 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,kk))^2;
                        for mm = 1:1+jj*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp9 = tmp9 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2;
                                tmp10 = tmp10 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2;
                            end
                        end
                    end
                end
                for ll = kk+1:1+ii*N
                    tmp12 = tmp12 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                    tmp15 = tmp15 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,kk))^2;
                    for mm = 1:1+jj*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp13 = tmp13 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,mm))^2;
                                tmp14 = tmp14 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2;
                                for nn = 1:1+ii*N
                                    if (nn ~= mm) && (nn ~= ll) && (nn ~= kk)
                                        tmp19 = tmp19 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,nn))^2*abs(Te(1+jj*N,mm))^2;
                                    end
                                end
                            end
                    end
                    for mm = 1:1+ii*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp16 = tmp16 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,kk))^2;
                                tmp17 = tmp17 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,ll))^2;
                                for nn = 1:1+jj*N
                                    if (nn ~= mm) && (nn ~= ll) && (nn ~= kk)
                                        tmp18 = tmp18 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,nn))^2;
                                    end
                                end
                            end
                    end
                end
            end
        end
    end
 end

H22E = mean(h22e)
H22E_test = tmp1 + tmp2 + tmp3 + tmp4 +  tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 +  tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19




tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
tmp8 = 0;
tmp9 = 0;
tmp10 = 0;
tmp11 = 0;
tmp12 = 0;
tmp13 = 0;
tmp14 = 0;
 for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj 
            for kk = 1:1+ii*N
                for ll = kk+1:1+ii*N
                    tmp1 = tmp1 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp2 = tmp2 + 4*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp3 = tmp3 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp4 = tmp4 + 2*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,kk))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp5 = tmp5 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll))*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    for mm = 1:1+ii*N
                        if (mm ~= kk) && (mm ~= ll)
                            tmp6 = tmp6 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,mm)*conj(Te(1+jj*N,mm))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                            for nn = mm+1:1+ii*N
                                if (nn ~= kk) && (nn ~= ll)
                                    tmp7 = tmp7 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll))*real(Te(1+ii*N,mm)*conj(Te(1+jj*N,mm))*conj(Te(1+ii*N,nn))*Te(1+jj*N,nn));
                                end
                            end
                        end
                    end
                    for mm = ll+1:1+ii*N
                        if (mm ~= kk) 
                            tmp8 = tmp8 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll))*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                            tmp9 = tmp9 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,mm))*Tb(1+jj*N,mm))*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                        end
                    end
                end
                for ll = 1:1+ii*N
                    if (ll ~= kk)
                        for mm = ll+1:1+ii*N
                            if mm ~= kk
                                tmp10 = tmp10 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*real(Te(1+ii*N,ll)*conj(Te(1+jj*N,ll))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                                tmp11 = tmp11 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,ll)*conj(Te(1+jj*N,ll))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                                tmp12 = tmp12 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,mm))*Tb(1+jj*N,mm))*real(Te(1+ii*N,ll)*conj(Te(1+jj*N,ll))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                            end
                        end
                    end
                end
                for ll = 1:1+jj*N
                    if (ll ~= kk)
                        for mm = ll+1:1+ii*N
                            if mm ~= kk
                                tmp13 = tmp13 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                            end
                        end
                        for mm = 1:1+ii*N
                            if (mm ~= kk) && (mm ~= ll)
                                for nn = mm+1 : 1+ii*N
                                    if (nn ~= kk) && (nn ~= ll)
                                        tmp14 = tmp14 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,mm)*conj(Te(1+jj*N,mm))*conj(Te(1+ii*N,nn))*Te(1+jj*N,nn));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
 end
H22E_bis = mean(h22e_bis)
H22E_bis_test = tmp1 + tmp2 + tmp3 + tmp4 +  tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 +  tmp13 + tmp14


%% 1 SDS DECODER. 
% 1 ∑_i |hb_{n+iN}|^4 * |he_{n+iN}|^4
% 
% H4H4 = mean(h4h4); % OK

figure;
plot(1./U^2*(4*U+H4E_test+H22E_test+H22E_bis_test),'-o') ; hold on; plot(e_avg_sym_decod2_e,':')

H4H4_test = 2*U; % OK
% % 
% % % 2 ∑_i ∑_j!=i   |hb_{n+iN}|^2*|he_{n+iN}|^2*|hb_{n+jN}|^2*|he_{n+jiN}|^2
% H2H2 = mean(h2h2); % OK
for dd = 1:length(b_subcar_b)
H2H2_test(:,dd) = modelCorrelH2H2(N,U,Tb(:,:,dd)); % Ok
end


for dd = 1:length(b_subcar_b)   
sinr_b_model(:,dd)    = sinrModelingFrequencyCorrelation(alpha,U,N,Tb(:,:,dd),snr_b,snr_e,"bob_correl").';
var_b_model(:,dd)     = modelVariance(alpha,N,U,Tb(:,:,dd),sinr_b_model(:,dd),snr_b,snr_e,"bob_correl");
capa_b_jensen(:,dd)       = log2(1+sinr_b_model(:,dd));
capa_b_new_approx(:,dd)   = log2(1+sinr_b_model(:,dd)) - 1/2*var_b_model(:,dd)./((1+sinr_b_model(:,dd)).^2);
end
capa_b_ergodic      = capa1_b_avg;


   
%% Eve Decoder 1
for dd = 1:length(b_subcar_b)   
sinr1_e_model(:,dd)    = sinrModelingFrequencyCorrelation(alpha,U,N,Tb(:,:,dd),snr_b,snr_e,"eve_decod1_correl").';
var1_e_model(:,dd)     = modelVariance(alpha,N,U,Tb(:,:,dd),sinr1_e_model(:,dd),snr_b,snr_e,"eve_decod1_correl");
capa1_e_jensen(:,dd)       = log2(1+sinr1_e_model(:,dd));
capa1_e_new_approx(:,dd)   = log2(1+sinr1_e_model(:,dd)) - 1/2*var1_e_model(:,dd)./((1+sinr1_e_model(:,dd)).^2);
end
capa1_e_ergodic      = capa1_e_avg;




%% SR
sr1_jensen = capa_b_jensen - capa1_e_jensen;
sr1_new_approx = capa_b_new_approx - capa1_e_new_approx;
sr1_ergodic = capa_b_ergodic - capa1_e_ergodic;



%% PLOT SECTION
set(0,'defaultAxesFontSize',19)
set(0,'DefaultLineMarkerSize',15);
set(0,'defaultLineLineWidth',2);

% 1. Approximation of the different approx

% 1.2. HIGH Correl as a fonction of alpha
figure; 
subplot(3,1,1)
plot(100*(1-alpha),capa_b_ergodic(:,1), 'Marker', 'square'); hold on; 
plot(100*(1-alpha),capa_b_jensen(:,1), 'Marker', 'o') ; hold on;
plot(100*(1-alpha),capa_b_new_approx(:,1), 'Marker', 'Te');
box on; grid on; 
xlabel('Percentage of AN injected ($\%$)')
ylabel('Capacity (bit/channel use)')
if length(b_subcar_b) == 2
    legend('Ergodic: $\Delta f_N / \Delta f_C = 1/6$', ...
       'First order : $\Delta f_N / \Delta f_C = 1/6$',...
       'Second order: $\Delta f_N / \Delta f_C = 1/6$','location','bestoutside')
else
    legend('Ergodic capacity: $E[\log_2(1+X)]$', ...
        'Jensen inequality: $\log_2(1+E[X])$',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
end
title('At Bob')


subplot(3,1,2) 
plot(100*(1-alpha),capa1_e_ergodic(:,1), 'Marker', 'square'); hold on; 
plot(100*(1-alpha),capa1_e_jensen(:,1), 'Marker', 'o') ; hold on; 
plot(100*(1-alpha),capa1_e_new_approx(:,1), 'Marker', 'Te');
box on; grid on; 
xlabel('Percentage of AN injected ($\%$)')
ylabel('Capacity (bit/channel use)')
if length(b_subcar_b) == 2
legend('Ergodic: $\Delta f_N / \Delta f_C = 1/6$', ...
       'First order : $\Delta f_N / \Delta f_C = 1/6$',...
       'Second order: $\Delta f_N / \Delta f_C = 1/6$','location','bestoutside')

else
    legend('Ergodic capacity: $E[\log_2(1+X)]$', ...
        'Jensen inequality: $\log_2(1+E[X])$',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
end   
title('At Eve: SDS Decoder')




subplot(3,1,3)
plot(100*(1-alpha),sr1_ergodic(:,1),'Marker', 'square'); hold on; 
plot(100*(1-alpha),sr1_jensen(:,1), 'Marker', 'o') ; hold on; 
plot(100*(1-alpha),sr1_new_approx(:,1),'Marker', 'Te');
box on; grid on; 
xlabel('Percentage of AN injected ($\%$)')
ylabel('Secrecy Rate (bit/channel use)')
if length(b_subcar_b) == 2
legend('Ergodic: $\Delta f_N / \Delta f_C = 1/6$', ...
        'First order: $\Delta f_N / \Delta f_C = 1/6$', ...
       'Second order: $\Delta f_N / \Delta f_C = 1/6$', ...
       'location','bestoutside')
else
    legend('Ergodic SR', ...
        'Jensen inequality SR',...
       'New approx SR','location','best')
end
title('SDS Decoder')



% % 1.3.No Correl as a fonction of alpha

figure; 
subplot(3,1,1)
plot(100*(1-alpha),capa_b_ergodic(:,2), 'Marker', 'square'); hold on; 
plot(100*(1-alpha),capa_b_jensen(:,2), 'Marker', 'o') ; hold on;
plot(100*(1-alpha),capa_b_new_approx(:,2), 'Marker', 'Te');
box on; grid on; 
xlabel('Percentage of AN injected ($\%$)')
ylabel('Capacity (bit/channel use)')
if length(b_subcar_b) == 2
    legend('Ergodic: $\Delta f_N / \Delta f_C = 100$', ...
       'First order : $\Delta f_N / \Delta f_C = 100$',...
       'Second order: $\Delta f_N / \Delta f_C = 100$','location','bestoutside')
else
    legend('Ergodic capacity: $E[\log_2(1+X)]$', ...
        'Jensen inequality: $\log_2(1+E[X])$',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
end
title('At Bob')


subplot(3,1,2) 
plot(100*(1-alpha),capa1_e_ergodic(:,2), 'Marker', 'square'); hold on; 
plot(100*(1-alpha),capa1_e_jensen(:,2), 'Marker', 'o') ; hold on; 
plot(100*(1-alpha),capa1_e_new_approx(:,2), 'Marker', 'Te');
box on; grid on; 
xlabel('Percentage of AN injected ($\%$)')
ylabel('Capacity (bit/channel use)')
if length(b_subcar_b) == 2
legend('Ergodic: $\Delta f_N / \Delta f_C = 100$', ...
       'First order : $\Delta f_N / \Delta f_C = 100$',...
       'Second order: $\Delta f_N / \Delta f_C = 100$','location','bestoutside')

else
    legend('Ergodic capacity: $E[\log_2(1+X)]$', ...
        'Jensen inequality: $\log_2(1+E[X])$',...
       'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
end   
title('At Eve: SDS Decoder')




subplot(3,1,3)
plot(100*(1-alpha),sr1_ergodic(:,2),'Marker', 'square'); hold on; 
plot(100*(1-alpha),sr1_jensen(:,2), 'Marker', 'o') ; hold on; 
plot(100*(1-alpha),sr1_new_approx(:,2),'Marker', 'Te');
box on; grid on; 
xlabel('Percentage of AN injected ($\%$)')
ylabel('Secrecy Rate (bit/channel use)')
if length(b_subcar_b) == 2
legend('Ergodic: $\Delta f_N / \Delta f_C = 100$', ...
        'First order: $\Delta f_N / \Delta f_C = 100$', ...
       'Second order: $\Delta f_N / \Delta f_C = 100$', ...
       'location','bestoutside')
else
    legend('Ergodic SR', ...
        'Jensen inequality SR',...
       'New approx SR','location','best')
end
title('SDS Decoder')



% Error plot as a function of the correlation at Bob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %
%%                         MODELIZATION SECTION                             %   
%                            Variance at Bob                                %
%                                                                           %
%                                                                           %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%

% % FOR H4H2H2
% % First term of H^2 H^2 with H^4 --> term4
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     tmp1 = tmp1 + 6*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*abs(Tb(1+kk*N,ll))^4;
%                     for mm = 1:1+kk*N
%                         if mm ~= ll
%                             tmp2 = tmp2 + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*abs(Tb(1+kk*N,mm))^4;
%                             tmp3 = tmp3 + 6*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*abs(Tb(1+kk*N,ll))^2*abs(Tb(1+kk*N,mm))^2;
%                             for nn = 1:1+kk*N
%                                 if (nn ~= mm) && (nn ~= ll)
%                                     tmp4 = tmp4 + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*abs(Tb(1+kk*N,mm))^2*abs(Tb(1+kk*N,nn))^2;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% TERM4_test = 8*(tmp1 + tmp2 + tmp3 + tmp4);
% TERM4 = mean(term4);            
% 
% 
% % Second term of H^2 H^2 with H^4 --> term22
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% tmp5 = 0;
% tmp6 = 0;
% tmp7 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 
%                 for ll = 1:1+ii*N
%                     for mm = 1:1+jj*N
%                         if mm ~= ll
%                             tmp1 = tmp1 + 3*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+kk*N,ll))^4*abs(Tb(1+jj*N,mm))^2;
%                             tmp2 = tmp2 + 3*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,mm))^4;
%                             tmp3 = tmp3 + 4*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+kk*N,ll))^2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,mm))^2;
%                             for nn = 1:1+kk*N
%                                 if (nn ~= mm) && (nn ~= ll)
%                                     tmp4 = tmp4 + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,nn))^4;
%                                     tmp5 = tmp5 + 4*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+kk*N,ll))^2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,nn))^2;
%                                     tmp6 = tmp6 + 4*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,mm))^2*abs(Tb(1+kk*N,nn))^2;
%                                     for oo = 1:1+kk*N
%                                         if (oo ~= nn) && (oo ~= mm) && (oo ~= ll)
%                                             tmp7 = tmp7 + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,nn))^2*abs(Tb(1+kk*N,oo))^2;
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end          
% TERM22_test = 4*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7);
% TERM22 = mean(term22);
% 
% 
% % Third term of H^2 H^2 with H^4 --> term22_bis
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% tmp5 = 0;
% tmp6 = 0;
% tmp7 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 
%                 for ll = 1:1+ii*N
%                     for mm = ll+1:1+jj*N
%                         if mm ~= ll
%                             tmp1 = tmp1 + 3*real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,ll))^4;
%                             tmp2 = tmp2 + 3*real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,mm))^4;
%                             tmp3 = tmp3 + 4*real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,ll))^2*abs(Tb(1+kk*N,mm))^2;
%                             for nn = 1:1+kk*N
%                                 if (nn ~= mm) && (nn ~= ll)
%                                     tmp4 = tmp4 + real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,nn))^4;
%                                     tmp5 = tmp5 + 4*real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,ll))^2*abs(Tb(1+kk*N,nn))^2;
%                                     tmp6 = tmp6 + 4*real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,mm))^2*abs(Tb(1+kk*N,nn))^2;
%                                     for oo = 1:1+kk*N
%                                         if (oo ~= nn) && (oo ~= mm) && (oo ~= ll)
%                                             tmp7 = tmp7 + real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm)))*abs(Tb(1+kk*N,nn))^2*abs(Tb(1+kk*N,oo))^2;
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end          
% TERM22_bis_test = 8*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7);
% TERM22_bis = mean(term22_bis);
% 
% 
% 
% % Fourth term of H^2 H^2 with H^4 --> term1111
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     for mm = ll+1:1+ii*N
%                         for nn = 1:1+jj*N
%                             for oo = nn+1:1+jj*N
%                                 if (oo ~= mm) && (oo ~= ll) && (nn ~= mm) && (nn ~= ll)
%                                     tmp = tmp + real(conj(Tb(1+ii*N,ll))*Tb(1+ii*N,mm)*Tb(1+jj*N,nn)*conj(Tb(1+jj*N,oo)) ...
%                                          * Tb(1+kk*N,ll)*conj(Tb(1+kk*N,mm))*conj(Tb(1+kk*N,nn))*Tb(1+jj*N,oo));
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% TERM1111_test = 32*tmp; % Very small 
% %TERM1111 = mean(term1111)
% 
% % Fifth term of H^2 H^2 with H^4 --> term211
% tmp1 = 0;
% tmp2 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     for mm = 1:1+jj*N
%                         for nn = mm+1:1+jj*N
%                             if (nn ~= ll) && (mm ~= ll)    
%                                 tmp1 = tmp1 + 2*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+kk*N,ll))^2* ...
%                                         real(conj(Tb(1+jj*N,mm))*Tb(1+kk*N,mm)*Tb(1+jj*N,nn)*conj(Tb(1+kk*N,nn)));    
%                             end
%                         end
%                     end
%                 end
%                 for ll = 1:1+ii*N
%                     for mm = 1:1+kk*N
%                         for nn = 1:1+jj*N
%                             for oo = nn+1:1+jj*N
%                                 if (oo ~= mm) && (oo ~= ll) && (nn ~= ll) && (nn ~= ll) && (mm ~= ll)
%                                     tmp2 = tmp2 + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+kk*N,mm))^2* ...
%                                         real(conj(Tb(1+jj*N,nn))*Tb(1+kk*N,nn)*Tb(1+jj*N,oo)*conj(Tb(1+kk*N,oo)));  
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end      
% TERM211_test = 16*(tmp1 + tmp2);
% TERM211 = mean(term211);
% 
% % Sixth term of H^2 H^2 with H^4 --> term211_bis
% tmp1 = 0;
% tmp2 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+jj*N
%                     for mm = 1:1+ii*N
%                         for nn = mm+1:1+ii*N
%                             if (nn ~= ll) && (mm ~= ll)    
%                                 tmp1 = tmp1 + 2*abs(Tb(1+jj*N,ll))^2*abs(Tb(1+kk*N,ll))^2* ...
%                                         real(conj(Tb(1+ii*N,mm))*Tb(1+kk*N,mm)*Tb(1+ii*N,nn)*conj(Tb(1+kk*N,nn)));    
%                             end
%                         end
%                     end
%                 end
%                 for ll = 1:1+jj*N
%                     for mm = 1:1+kk*N
%                         for nn = 1:1+ii*N
%                             for oo = nn+1:1+ii*N
%                                 if (oo ~= mm) && (oo ~= ll) && (nn ~= ll) && (nn ~= ll) && (mm ~= ll)
%                                     tmp2 = tmp2 + abs(Tb(1+jj*N,ll))^2*abs(Tb(1+kk*N,mm))^2* ...
%                                         real(conj(Tb(1+ii*N,nn))*Tb(1+kk*N,nn)*Tb(1+ii*N,oo)*conj(Tb(1+kk*N,oo)));  
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end      
% TERM211_bis_test = 16*(tmp1 + tmp2);
% TERM211_bis = mean(term211_bis);




% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     tmp = tmp + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2 ...
%                         *abs(Hw_TX(ll,ll))^4*abs(Hb_TX(1+kk*N,1+kk*N))^4;
%                 end
%             end
%         end
%         
%     end
% end
% term4(iter,dd) = 2*tmp;
% 
% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     for mm = 1:1+jj*N
%                         if ll ~= mm
%                             tmp = tmp + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,mm))^2 ...
%                                   *abs(Hw_TX(ll,ll))^2*abs(Hw_TX(mm,mm))^2*abs(Hb_TX(1+kk*N,1+kk*N))^4;
%                         end
%                     end
%                 end
%             end
%         end
%         
%     end
% end
% term22(iter,dd) = 2*tmp;
% 
% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     for mm = ll+1:1+jj*N
%                         if ll ~= mm
%                             tmp = tmp + real(conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll)*Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))) ...
%                                   *abs(Hw_TX(ll,ll))^2*abs(Hw_TX(mm,mm))^2*abs(Hb_TX(1+kk*N,1+kk*N))^4;
%                         end
%                     end
%                 end
%             end
%         end
%         
%     end
% end
% term22_bis(iter,dd) = 4*tmp;



% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     for mm = 1:1+ii*N
%                         for nn = 1:1+jj*N
%                             for oo = 1:1+jj*N
%                                 if  (oo ~= nn) && (oo ~= mm) && (oo ~= ll) && (nn ~= mm) && (nn ~= ll) && (mm ~= ll)
%                                     tmp = tmp + conj(Tb(1+ii*N,ll))*Tb(1+ii*N,mm)*Tb(1+jj*N,nn)*conj(Tb(1+jj*N,oo))* ...
%                                         Hw_TX(ll,ll)*Hw_RX(mm,mm)*Hw_TX(nn,nn)*Hw_RX(oo,oo)*abs(Hb_TX(1+kk*N,1+kk*N))^4;    
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end              
% term1111(iter,dd) = 2*tmp;


% 
% 
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+ii*N
%                     for mm = 1:1+jj*N
%                         for nn = 1:1+jj*N
%                             if (nn ~= mm) && (nn ~= ll) && (mm ~= ll)    
%                                 tmp = tmp + abs(Tb(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*conj(Tb(1+jj*N,mm))*Tb(1+jj*N,nn)* ...
%                                             Hw_TX(mm,mm)*Hw_RX(nn,nn)*abs(Hb_TX(1+kk*N,1+kk*N))^4;    
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end              
% term211(iter,dd) = 2*tmp;



% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 0:U-1
%             if (kk ~= ii) && (kk ~= jj)
%                 for ll = 1:1+jj*N
%                     for mm = 1:1+ii*N
%                         for nn = 1:1+ii*N
%                             if (nn ~= mm) && (nn ~= ll) && (mm ~= ll)    
%                                 tmp = tmp + abs(Tb(1+jj*N,ll))^2*abs(Hw_TX(ll,ll))^2*conj(Tb(1+ii*N,mm))*Tb(1+ii*N,nn)* ...
%                                             Hw_TX(mm,mm)*Hw_RX(nn,nn)*abs(Hb_TX(1+kk*N,1+kk*N))^4;    
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end              
% term211_bis(iter,dd) = 2*tmp;







% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for zz = 0:U-1
%             for kk = 1:1+ii*N
%                 tmp1 = tmp1 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Hw_TX(kk,kk))^4*abs(Hb_TX(1+zz*N,1+zz*N))^4;
%                 for ll = 1:1+jj*N
%                     if ll ~= kk
%                         tmp2 = tmp2 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Hw_TX(kk,kk))^2*abs(Hw_TX(ll,ll))^2*abs(Hb_TX(1+zz*N,1+zz*N))^4;
%                     end
%                 end
%             end
%             for kk = 1:1+ii*N
%                 for ll = kk+1:1+ii*N
%                     tmp3 = tmp3 + 2*real(conj(Tb(1+ii*N,kk))*Tb(1+jj*N,kk)*Tb(1+ii*N,ll)*conj(Tb(1+jj*N,ll)))...
%                            *abs(Hw_TX(kk,kk))^2*abs(Hw_TX(ll,ll))^2*abs(Hb_TX(1+zz*N,1+zz*N))^4;
%                 end
%             end
%         end
%     end
% end
% h4h2h2_bis(iter,dd) = 2*(tmp1+tmp2+tmp3);
% 



% % FOR H2H2H2H2
% % First term of hb_i^2 hb_j^2 with hb_k^2 hb_l^2 --> term4
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 for mm = 1:1+ii*N
%                     tmp = tmp + abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2 *abs(Hw_TX(mm,mm))^4* ...
%                           abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%                 end
%             end
%         end
%         
%     end
% end
% term4(iter,dd) = 24*tmp;
% 
% 
% % Second term of hb_i^2 hb_j^2 with hb_k^2 hb_l^2 --> term22
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 for mm = 1:1+ii*N
%                     for nn = 1:1+jj*N
%                         if nn ~= mm
%                             tmp = tmp + abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*abs(Hw_TX(mm,mm))^2*abs(Hw_TX(nn,nn))^2* ...
%                             abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%                         end
%                     end
%                 end
%             end
%         end
%         
%     end
% end
% term22(iter,dd) = 24*tmp;
% 
% 
% % Third term of hb_i^2 hb_j^2 with hb_k^2 hb_l^2 --> term22_bis
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 for mm = 1:1+ii*N
%                     for nn = mm+1:1+ii*N
%                             tmp = tmp + real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                 abs(Hw_TX(mm,mm))^2*abs(Hw_TX(nn,nn))^2* ...
%                                 abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%                     end
%                 end
%             end
%         end
%         
%     end
% end
% term22_bis(iter,dd) = 48*tmp;
% 
% 
% 
% 
% % Fourth term of hb_i^2 hb_j^2 with hb_k^2 hb_l^2 --> term1111
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 for mm = 1:1+ii*N
%                     for nn = 1:1+ii*N
%                         for oo= 1:1+jj*N
%                             for pp = 1:1+jj*N
%                                 if (pp ~= oo) && (pp ~= nn) && (pp ~= mm) && (oo ~= nn) && (oo ~= mm) && (mm ~= nn)
%                                     tmp = tmp + conj(Tb(1+ii*N,mm))*Tb(1+ii*N,nn)*conj(Tb(1+jj*N,oo))*Tb(1+jj*N,pp)*...
%                                         Hw_RX(mm,mm)*Hw_TX(nn,nn)*Hw_RX(oo,oo)*Hw_TX(pp,pp)* ...
%                                         abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% term1111(iter,dd) = 24*tmp; 
% 
% 
% % Fith term of hb_i^2 hb_j^2 with hb_k^2 hb_l^2 --> term211
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 
%                 for mm = 1:1+ii*N
%                     for nn = 1:1+jj*N
%                         for oo= 1:1+jj*N
%                             if (oo ~= nn) && (oo ~= mm) && (mm ~= nn)
%                                 tmp = tmp + abs(Tb(1+ii*N,mm))^2*abs(Hw_TX(mm,mm))^2* ...
%                                     Tb(1+jj*N,nn)*conj(Tb(1+jj*N,oo))*Hw_RX(nn,nn)*Hw_TX(oo,oo)* ...
%                                     abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end
% term211(iter,dd) = 24*tmp; 
% 
% 
% % Fith term of hb_i^2 hb_j^2 with hb_k^2 hb_l^2 --> term211_bis
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1
%                 
%                 for mm = 1:1+jj*N
%                     for nn = 1:1+ii*N
%                         for oo= 1:1+ii*N
%                             if (oo ~= nn) && (oo ~= mm) && (mm ~= nn)
%                                 tmp = tmp + abs(Tb(1+jj*N,mm))^2*abs(Hw_TX(mm,mm))^2* ...
%                                     Tb(1+ii*N,nn)*conj(Tb(1+ii*N,oo))*Hw_RX(nn,nn)*Hw_TX(oo,oo)* ...
%                                     abs(Hb_TX(1+kk*N,1+kk*N))^2*abs(Hb_TX(1+ll*N,1+ll*N))^2;
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end
% term211_bis(iter,dd) = 24*tmp; 
%
% % First term of hb_i^2 hb_j^2 (with h^4) multiplied by hb_k^2 hb_l^2  
% % -> term4
% 
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% tmp5 = 0;
% tmp6 = 0;
% tmp7 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1   %%%%%%%%
%                 
%                 for mm = 1:1+ii*N
%                     tmp1 = tmp1 + 12*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                      abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,mm))^2;
%                     for nn = 1:1+kk*N
%                         if nn ~= mm
%                             tmp2 = tmp2 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                             abs(Tb(1+kk*N,nn))^2*abs(Tb(1+ll*N,nn))^2;
%                             tmp3 = tmp3 + 3*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                             abs(Tb(1+ll*N,mm))^2*abs(Tb(1+kk*N,nn))^2;
%                             for oo = 1:1+ll*N
%                                 if (oo ~= nn) && (oo ~= mm)
%                                     tmp4 = tmp4 + abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                                   abs(Tb(1+kk*N,nn))^2*abs(Tb(1+ll*N,oo))^2;
%                                 end
%                             end
%                             for oo = nn+1:1+kk*N
%                                 if oo ~= mm
%                                     tmp7 = tmp7 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                                     real(Tb(1+kk*N,nn)*conj(Tb(1+ll*N,nn))*conj(Tb(1+kk*N,oo))*Tb(1+ll*N,oo));
%                                 end
%                             end
%                         end
%                     end
%                     for nn = 1:1+ll*N
%                         if nn ~= mm
%                             tmp5 = tmp5 + 3*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                             abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,nn))^2;
%                         end
%                     end
%                     for nn = mm+1:1+kk*N
%                         tmp6 = tmp6 + 12*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,mm))^2*...
%                                          real(Tb(1+kk*N,mm)*conj(Tb(1+ll*N,mm))*conj(Tb(1+kk*N,nn))*Tb(1+ll*N,nn));
%                     end
%                   
%                 end             %%%%%%%%
%              end
%         end
%         
%     end
% end
% 
% TERM4_test = 48*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7);  % 24x2 --> looks ok
% TERM4 = mean(term4);
% 
% % Second term of hb_i^2 hb_j^2 (with h^2 h^2) multiplied by hb_k^2 hb_l^2  
% % -> term22
% 
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% tmp5 = 0;
% tmp6 = 0;
% tmp7 = 0;
% tmp8 = 0;
% tmp9 = 0;
% tmp10 = 0;
% tmp11 = 0;
% tmp12 = 0;
% tmp13 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1   %%%%%%%%
%                 
%                 for mm = 1:1+ii*N
%                     for nn = 1:1+jj*N
%                         if nn ~= mm
%                             tmp1 = tmp1 + 6*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                             abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,mm))^2;
%                             tmp2 = tmp2 + 6*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                             abs(Tb(1+kk*N,nn))^2*abs(Tb(1+ll*N,nn))^2;
%                             tmp3 = tmp3 + 4*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                             abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,nn))^2;
%                            
%                             for oo = 1:1+kk*N
%                                if (oo ~= nn) && (oo ~= mm)
%                                    tmp4 = tmp4 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                    abs(Tb(1+kk*N,oo))^2*abs(Tb(1+ll*N,oo))^2;
%                                    tmp5 = tmp5 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                    abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,oo))^2;
%                                    tmp6 = tmp4;
%                                    for pp = 1:1+ll*N
%                                        if (pp ~= oo) && (pp ~= nn) && (pp ~= mm)
%                                            tmp9 = tmp9 + abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                          abs(Tb(1+kk*N,oo))^2*abs(Tb(1+ll*N,pp))^2;
%                                        end
%                                    end
%                                    for pp = oo+1:1+kk*N
%                                        tmp12 = tmp12 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                         real(Tb(1+kk*N,oo)*conj(Tb(1+ll*N,oo))*conj(Tb(1+kk*N,pp))*Tb(1+ll*N,pp));
%                                    end
%                                end
%                             end
%                             for oo = 1:1+ll*N
%                                 if (oo ~= nn) && (oo ~= mm)
%                                     tmp7 = tmp7 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                     abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,oo))^2;
%                                     tmp8 = tmp8 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                     abs(Tb(1+kk*N,nn))^2*abs(Tb(1+ll*N,oo))^2;
%                                 end
%                             end
%                         end
%                         if nn ~= mm
%                             for oo = mm+1:1+kk*N
%                                 if oo ~= nn
%                                     tmp10 = tmp10 + 8*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                 real(conj(Tb(1+kk*N,mm))*Tb(1+ll*N,mm)*Tb(1+kk*N,oo)*conj(Tb(1+ll*N,oo)));
%                                 end
%                             end
%                             for oo = nn+1:1+kk*N
%                                 if oo ~= mm
%                                     tmp11 = tmp11 + 8*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                                 real(conj(Tb(1+kk*N,nn))*Tb(1+ll*N,nn)*Tb(1+kk*N,oo)*conj(Tb(1+ll*N,oo)));
%                                 end
%                             end
%                         end
%                     end
%                     for nn = mm+1:1+jj*N
%                         tmp13 = tmp13 + 8*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+jj*N,nn))^2*...
%                                         real(Tb(1+kk*N,mm)*conj(Tb(1+ll*N,mm))*conj(Tb(1+kk*N,nn))*Tb(1+ll*N,nn));
%                     end
%               
%                 end             %%%%%%%%
%              end
%         end
%         
%     end
% end
% 
% TERM22_test = 24*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 + tmp13);
% 
% TERM22 = mean(term22);
% 
% 
% % Third term of hb_i^2 hb_j^2 (with conj(Tb Tb^* Tb Tb^*) h^2 h^2) multiplied
% % by hb_k^2 hb_l^2 -> term22_bis
% 
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% tmp5 = 0;
% tmp6 = 0;
% tmp7 = 0;
% tmp8 = 0;
% tmp9 = 0;
% tmp10 = 0;
% tmp11 = 0;
% tmp12 = 0;
% tmp13 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1   %%%%%%%%
%                 
%                 for mm = 1:1+ii*N
%                     for nn = mm+1:1+ii*N
%                             tmp1 = tmp1 + 6*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                             abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,mm))^2;
%                             tmp2 = tmp2 + 6*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                             abs(Tb(1+kk*N,nn))^2*abs(Tb(1+ll*N,nn))^2;
%                             tmp3 = tmp3 + 4*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                             abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,nn))^2;
%                            
%                             for oo = 1:1+kk*N
%                                if (oo ~= nn) && (oo ~= mm)
%                                    tmp4 = tmp4 + 2*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                    abs(Tb(1+kk*N,oo))^2*abs(Tb(1+ll*N,oo))^2;
%                                    tmp5 = tmp5 + 2*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                    abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,oo))^2;
%                                    tmp6 = tmp4;
%                                    for pp = 1:1+ll*N
%                                        if (pp ~= oo) && (pp ~= nn) && (pp ~= mm)
%                                            tmp9 = tmp9 + real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                          abs(Tb(1+kk*N,oo))^2*abs(Tb(1+ll*N,pp))^2;
%                                        end
%                                    end
%                                    for pp = oo+1:1+kk*N
%                                        tmp12 = tmp12 + 2*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                         real(Tb(1+kk*N,oo)*conj(Tb(1+ll*N,oo))*conj(Tb(1+kk*N,pp))*Tb(1+ll*N,pp));
%                                    end
%                                end
%                             end
%                             for oo = 1:1+ll*N
%                                 if (oo ~= nn) && (oo ~= mm)
%                                     tmp7 = tmp7 + 2*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                     abs(Tb(1+kk*N,mm))^2*abs(Tb(1+ll*N,oo))^2;
%                                     tmp8 = tmp8 + 2*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                     abs(Tb(1+kk*N,nn))^2*abs(Tb(1+ll*N,oo))^2;
%                                 end
%                             end
%                             for oo = mm+1:1+kk*N
%                                 if oo ~= nn
%                                     tmp10 = tmp10 + 8*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                 real(conj(Tb(1+kk*N,mm))*Tb(1+ll*N,mm)*Tb(1+kk*N,oo)*conj(Tb(1+ll*N,oo)));
%                                 end
%                             end
%                             for oo = nn+1:1+kk*N
%                                 if oo ~= mm
%                                     tmp11 = tmp11 + 8*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                                 real(conj(Tb(1+kk*N,nn))*Tb(1+ll*N,nn)*Tb(1+kk*N,oo)*conj(Tb(1+ll*N,oo)));
%                                 end
%                             end
%                     end
%                     for nn = mm+1:1+ii*N
%                         tmp13 = tmp13 + 8*real(Tb(1+ii*N,mm)*conj(Tb(1+jj*N,mm))*conj(Tb(1+ii*N,nn))*Tb(1+jj*N,nn))*...
%                                         real(Tb(1+kk*N,mm)*conj(Tb(1+ll*N,mm))*conj(Tb(1+kk*N,nn))*Tb(1+ll*N,nn));
%                     end
%               
%                 end             %%%%%%%%
%              end
%         end
%         
%     end
% end
% 
% TERM22_bis_test = 48*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 + tmp13);
% TERM22_bis = mean(term22_bis);
% 
% 
% % 4 Fourth term of hb_i^2 hb_j^2 (with h h^* h h^*) multiplied 
% % by hb_k^2 hb_l^2  -> term1111
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1   %%%%%%%%
%                 for mm = 1:1+ii*N
%                     for nn = mm+1:1+ii*N
%                         for oo = 1:1+jj*N
%                             if (oo ~= nn) && (oo ~= mm)
%                                 for pp = oo+1:1+jj*N
%                                     if (pp ~= nn) && (pp ~= mm)
%                                         tmp = tmp + real(conj(Tb(1+ii*N,mm))*Tb(1+kk*N,mm)*Tb(1+ii*N,nn)*conj(Tb(1+kk*N,nn))* ... 
%                                                          conj(Tb(1+jj*N,oo))*Tb(1+ll*N,oo)*Tb(1+jj*N,pp)*conj(Tb(1+ll*N,pp)));
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% TERM1111_test = 16*24*tmp;
% TERM1111 = mean(term1111);
% 
% 
% % 5 Fifth term of hb_i^2 hb_j^2 (with h^2 h^* h) multiplied 
% % by hb_k^2 hb_l^2  -> term211
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1   %%%%%%%%
%                 
%                 for mm = 1:1+ii*N
%                     for nn = 1:1+jj*N
%                         for oo = nn+1:1+jj*N
%                             if (oo ~= mm) && (nn ~= mm)
%                                 tmp1 = tmp1 + 4*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+kk*N,mm))^2* ... 
%                                        real(conj(Tb(1+jj*N,nn)) * Tb(1+ll*N,nn) * Tb(1+jj*N,oo) * conj(Tb(1+ll*N,oo)));
%                             end
%                         end
%                     end
%                     for nn = 1:1+kk*N
%                         for oo = 1:1+jj*N
%                             for pp = oo+1:1+jj*N
%                                if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
%                                    tmp2 = tmp2 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+kk*N,nn))^2* ... 
%                                        real(conj(Tb(1+jj*N,oo)) * Tb(1+ll*N,oo) * Tb(1+jj*N,pp) * conj(Tb(1+ll*N,pp)));
%                                end
%                             end
%                         end
%                     end
%                     for nn = 1:1+jj*N
%                         for oo = nn+1:1+jj*N
%                             if (oo ~= mm) && (nn ~= mm)
%                                 tmp3 = tmp3 + 4*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+ll*N,mm))^2* ... 
%                                        real(conj(Tb(1+jj*N,nn)) * Tb(1+kk*N,nn) * Tb(1+jj*N,oo) * conj(Tb(1+kk*N,oo)));
%                             end
%                         end
%                     end
%                     for nn = 1:1+ll*N
%                         for oo = 1:1+jj*N
%                             for pp = oo+1:1+jj*N
%                                if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
%                                    tmp4 = tmp4 + 2*abs(Tb(1+ii*N,mm))^2*abs(Tb(1+ll*N,nn))^2* ... 
%                                        real(conj(Tb(1+jj*N,oo)) * Tb(1+kk*N,oo) * Tb(1+jj*N,pp) * conj(Tb(1+kk*N,pp)));
%                                end
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end
% TERM211_test = 48*(tmp1 + tmp2 + tmp3 + tmp4);
% TERM211 = mean(term211);
% 
% % 6 Sixth term of hb_i^2 hb_j^2 (with h^2 h^* h) multiplied 
% % by hb_k^2 hb_l^2  -> term211_bis
% tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
% tmp4 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = jj+1:U-1
%             for ll = kk+1:U-1   %%%%%%%%Tb(1+ii*N,nn)*conj(Tb(1+ii*N,oo))
%                 
%                 for mm = 1:1+jj*N
%                     for nn = 1:1+ii*N
%                         for oo = nn+1:1+ii*N
%                             if (oo ~= mm) && (nn ~= mm)
%                                 tmp1 = tmp1 + 4*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,mm))^2* ... 
%                                        real(conj(Tb(1+ii*N,nn)) * Tb(1+ll*N,nn) * Tb(1+ii*N,oo) * conj(Tb(1+ll*N,oo)));
%                             end
%                         end
%                     end
%                     for nn = 1:1+kk*N
%                         for oo = 1:1+ii*N
%                             for pp = oo+1:1+ii*N
%                                if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
%                                    tmp2 = tmp2 + 2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+kk*N,nn))^2* ... 
%                                        real(conj(Tb(1+ii*N,oo)) * Tb(1+ll*N,oo) * Tb(1+ii*N,pp) * conj(Tb(1+ll*N,pp)));
%                                end
%                             end
%                         end
%                     end
%                     for nn = 1:1+ii*N
%                         for oo = nn+1:1+ii*N
%                             if (oo ~= mm) && (nn ~= mm)
%                                 tmp3 = tmp3 + 4*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+ll*N,mm))^2* ... 
%                                        real(conj(Tb(1+ii*N,nn)) * Tb(1+kk*N,nn) * Tb(1+ii*N,oo) * conj(Tb(1+kk*N,oo)));
%                             end
%                         end
%                     end
%                     for nn = 1:1+ll*N
%                         for oo = 1:1+ii*N
%                             for pp = oo+1:1+ii*N
%                                if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
%                                    tmp4 = tmp4 + 2*abs(Tb(1+jj*N,mm))^2*abs(Tb(1+ll*N,nn))^2* ... 
%                                        real(conj(Tb(1+ii*N,oo)) * Tb(1+kk*N,oo) * Tb(1+ii*N,pp) * conj(Tb(1+kk*N,pp)));
%                                end
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end
% TERM211_bis_test = 48*(tmp1 + tmp2 + tmp3 + tmp4);
% TERM211_bis = mean(term211_bis);






%
%

% First term of |hb_i|^4 with |hb_j|^4 --> term4
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 1:1+ii*N
%             tmp = tmp + abs(Tb(1+ii*N,kk))^4*abs(Hw_TX(kk,kk))^4*abs(Hb_TX(1+jj*N,1+jj*N))^4;
%         end
%     end
% end
% term4(iter,dd) = tmp;
% 
% 
% % Third term of |hb_i|^4 with |hb_j|^4 --> term22
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 1:1+ii*N
%             for ll = 1:1+ii*N
%                 if ll ~= kk
%                     tmp = tmp + 2*abs(Tb(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
%                         abs(Tb(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*abs(Hb_TX(1+jj*N,1+jj*N))^4;
%                 end
%             end
%         end
%     end
% end
% term22(iter,dd) = tmp;
% 
% % tmp = 0;
% % for ii = 0:U-1
% %     for jj = ii+1:U-1
% %         for kk = 1:1+ii*N
% %             for ll = 1:1+ii*N
% %                 if ll ~= kk
% %                     for mm = 1:1+jj*N
% %                         tmp = tmp + 2*abs(Tb(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
% %                         abs(Tb(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*abs(Tb(1+jj*N,mm))^4*abs(Hw_TX(mm,mm))^4;
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% % term22_tmp1_tmp3(iter,dd) = tmp;
% % 
% % tmp = 0;
% % for ii = 0:U-1
% %     for jj = ii+1:U-1
% %         for kk = 1:1+ii*N
% %             for ll = 1:1+ii*N
% %                 if ll ~= kk
% %                     for mm = 1:1+jj*N
% %                         for nn = 1:1+jj*N
% %                             if mm ~= nn
% %                                 tmp = tmp + 2*abs(Tb(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
% %                                 abs(Tb(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*...
% %                                 2*abs(Tb(1+jj*N,mm))^2*abs(Hw_TX(mm,mm))^2*abs(Tb(1+jj*N,nn))^2*abs(Hw_TX(nn,nn))^2;
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% % term22_tmp2_tmp4_tmp5(iter,dd) = tmp;
% 
% 
% % Fourth term of  |hb_i|^4 with |hb_j|^4 --> term211
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 1:1+ii*N
%             for ll = 1:1+ii*N
%                 for mm = 1:1+ii*N
%                     if (mm ~= ll) && (mm ~= kk) && (kk ~= ll)
%                         tmp = tmp + 2*abs(Tb(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
%                         conj(Tb(1+ii*N,ll))*Hw_RX(ll,ll)*Tb(1+ii*N,mm)*Hw_TX(mm,mm)*abs(Hb_TX(1+jj*N,1+jj*N))^4;
%                     end
%                 end
%             end
%         end
%     end
% end
% term211(iter,dd) = tmp;
% 
% 
% 
% 
% 
% 
% 
% 

%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % No correlation: 
% % Approximation of E[ 1+ log(X) ] = log(1 + E[X]) - var(X)/( 2( 1+E[X] )^2 )
% % where X = SINR_b ou SINR_e
% % 1 At Bob
% % 1.1: var(X) = E[X^2] - (E[X])^2
% % E[X^2] = E[ |S^H |Hb|^2 S|^4 / |S^H v_b|^4 ]
% % 1.1.1 E[ | S^H |Hb|^2 S |^4  ]
% 
% sym_b_TEST = mean(squeeze(mean(abs(sym_b_test).^4,2)),1);
% SYM_B_TEST = alpha.^2./U^4*(24*U + 36*U*(U-1) + 12*U*(U-1)*(U-2) + U*(U-1)*(U-2)*(U-3));  %-> OK si no correl
% % figure;
% % plot(sym_b_TEST); hold on; plot(SYM_B_TEST,'o');
% % legend('Bob square energy symbol simu','Bob square energy symbol analytic')
% 
% % 1.1.2: E[ |S^H v_b|^4 ] 
% noise_b_TEST = mean(squeeze(mean(abs(noise_b_test).^4,2)),1);
% NOISE_B_TEST = 2*sigma_b^2;                                                 %-> OK si no correl
% % figure;
% % plot(noise_b_TEST); hold on; plot(ones(length(alpha)).*NOISE_B_TEST,'o');   % -> OK si no correl
% % legend('Bob square energy noise simu','Bob square energy noise analytic')
% 
% % --> (E[X])^2 
% SINRb = alpha.*(U+1)./(U.*sigma_b);
% SINRsquareb = SINRb.^2;
% 
% % --> variance @Bob : var(X) = E[X^2] - (E[X])^2
% varTESTb = SYM_B_TEST./NOISE_B_TEST - SINRsquareb;
% 
% % 1.2 New capa approximation New approx:
% TESTcapab = log2(1+SINRb) - varTESTb./(2*(1+SINRb).^2); % New approx of Bob capa 
% 
% 
% 
% 
% % 2 At Eve: SDS
% % 2.1: var(X) = E[X^2] - (E[X])^2
% % E[X^2] = E[ |S^H Hb^* He S|^4 / |S^H v_e + S^H He w|^4 ]
% % 1.1.1 E[ | S^H Hb^* He S |^4  ]
% sym_e_TEST1 = mean(squeeze(mean(abs(sym_e_test1).^4,2)),1);
% SYM_E_TEST1 = 2*alpha.^2.*(U+1)./U^3;
% 
% % 1.1.2 E[ |S^H v_e + S^H He w |^4 ]
% denom_e_TEST1 = mean(squeeze(mean(abs(denom_e_test1).^4,2)),1);
% DENOM_E_TEST1 = 1./U*(2*sigma_e^2 + 4./U^2.*(1-alpha).^2 + 2*(U-1)*sigma_e^2 ...
%                 + 2./U^2.*(1-alpha).^2.*(U-1)+4./U.*(1-alpha).*sigma_e ...
%                 + 4./U.*(1-alpha).*(U-1).*sigma_e);
% 
% % --> (E[X])^2 
% SINRe1 = alpha./(U.*sigma_e + (1-alpha));
% SINRsquare_e1 = SINRe1.^2;
% 
% % --> variance @Bob : var(X) = E[X^2] - (E[X])^2
% varTESTe1 = SYM_E_TEST1./DENOM_E_TEST1 - SINRsquare_e1;
% 
% % 1.2 New capa approximation New approx:
% TESTcapae1 = log2(1+SINRe1) - varTESTe1./(2*(1+SINRe1).^2); % New approx of Bob capa 
% 
% 
% 
% 
% % 3. Test different capacity approximation
% A = 1+sinr1_correl_b;
% B = sinr1_correl_b; 
% C = log2(1+sinr1_correl_b);
% 
% E = 1+sinr1_e;
% F = sinr1_e; 
% G = log2(1+sinr1_e);
% 
% 
% TESTA = log2(squeeze(mean(A)));   % log2(E[1+sinr]) --> not computed 
% TESTB = log2(1+squeeze(mean(B))); % log2(1+E[sinr]) --> Jensen's aproximation
% TESTC = squeeze(mean(C));         % E[log2(1+sinr)] --> exact ergodic capa
% 
% TESTE = log2(squeeze(mean(E)));   % log2(E[1+sinr]) --> not computed
% TESTF = log2(1+squeeze(mean(F))); % log2(1+E[sinr]) --> Jensen's aproximation
% TESTG = squeeze(mean(G));         % E[log2(1+sinr)] --> exact ergodic capa
% 
% 
% 
% figure; 
% plot(TESTC,'Marker','square'); hold on;  plot(TESTB,'-o'); hold on; plot(TESTcapab, 'Marker','diamond')
% title('Computation of capacity BOB')
% legend('Ergodic capacity: $E[\log_2(1+X)]$','Jensen inequality: $\log_2(1+E[X])$ ',...
%        'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
% box on; grid on;
% ylabel('Capacity (bit/channel use)')
% 
% figure; 
% plot(TESTG,'Marker','square'); hold on; plot(TESTF,'-o'); hold on; plot(TESTcapae1, 'Marker','diamond')
% title('Computation of capacity Eve')
% legend('Ergodic capacity: $E[\log_2(1+X)]$','Jensen inequality: $\log_2(1+E[X])$ ',...
%        'New approx: $\log_2(1+E[X]) - \frac{var(X)}{2(1+E[X])^2}$','location','best')
% box on; grid on;
% ylabel('Capacity (bit/channel use)')
% 
% 
% 
% %%%%%%%
% %%%%%%%
% %%%%%%%
% %%%%%%%
% %%%%%%%
% %%%%%%%
% %%%%%%%
% %%%%%%%
% %%%%%%%
% 
% 
% %% MODELING
% 
% % SINR modeling
% % sinr1_decorrel_model_b = sinrModeling(alpha,U,snr_b,snr_e,1,1,1,"bob_decod1");
% for dd = 1:length(b_subcar_b)
% sinr1_correl_model_b(:,dd) = sinrModelingFrequencyCorrelation(alpha,U,N,squeeze(Tb(:,:,dd)),snr_b,snr_e,"bob_correl");
% 
% sinr1_model_e(:,dd) = sinrModelingFrequencyCorrelation(alpha,U,N,Tb(:,:,dd),snr_b,snr_e,"eve_decod1_correl");
% sinr2_model_e(:,dd) = sinrModelingFrequencyCorrelation(alpha,U,N,Tb(:,:,dd),snr_b,snr_e,"eve_decod2_correl");
% end
% 
% % Capacity modeling
% 
% capa1_correl_model_b = log2(1+sinr1_correl_model_b);
% % capa1_decorrel_model_b = log2(1+sinr1_decorrel_model_b);
% 
% capa1_model_e = log2(1+sinr1_model_e);
% capa2_model_e = log2(1+sinr2_model_e);
% 
% % Secrecy Rate modeling (only for correlated Bob's channel)
% 
% sr1_model = capa1_correl_model_b - capa1_model_e;
% sr2_model = capa1_correl_model_b - capa2_model_e;
% 
% 
% %% PLOT
% 
% 
% % 1: Check each term of Bob SINR
% % %% TEST OF INDIVIDUAL TERM
% % 
% % % 1: Data energy at Bob  := e_avg_sym_decod1_correl_b
% % % PROBLEM QUAND U GRANDET GDE CORRELATION
% TT = Tb; 
% 
% 
% for dd = 1:length(b_subcar_b)
% Tb = TT(:,:,dd); 
% term_exponent_4 = 0;
% term_exponent_2 = 0;
% term_exponent_TEST = 0;
% for ii = 0: U-2
%     %term_exponent_TEST = term_exponent_TEST +  2*abs(Tb(1+ii*N,kk))^4;
%     for jj = ii+1 : U-1
%         %jj
%         %term_exponent_TEST = term_exponent_TEST +  2*abs(Tb(1+ii*N,kk))^4;
%         for kk = 1 : 1+ii*N
%              %[{1+ii*N,kk}  {1+jj*N,kk}]
%             term_exponent_4 = term_exponent_4 + 2*abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,kk))^2;
%         end
%     end
% end
% for ii = 0: U-2
%     for jj = ii+1 : U-1
%         for kk = 1 : 1+ii*N
%             for ll = 1 : 1+jj*N
%                 if ll ~= kk
%                     %[{1+ii*N,kk}  {1+jj*N,ll}]
%                     term_exponent_2 = term_exponent_2 + abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,ll))^2   ;
%                 end   
%             end
%         end
%     end
% end
% 
% term_indep = 0;
% for ii = 0: U-2 
%     for jj = ii+1 : U-1
%         for kk = 1 : 1+ii*N
%             for ll = kk+1 : 1+ii*N
%                 %if ll  kk
%                     %[{1+ii*N,kk}  {1+ii*N,ll} , {1+jj*N,kk}  {1+jj*N,ll}]
%                     term_indep = term_indep + 2*real(conj(Tb(1+ii*N,kk)) * Tb(1+ii*N,ll) * Tb(1+jj*N,kk) * conj(Tb(1+jj*N,ll)))   ;
%                 %end   
%             end
%         end
%     end
% end
% 
% 
% 
% 
% sum_double_square(dd) = 2*(term_exponent_4 + term_exponent_2 + term_indep);                  % ∑_i ∑_j |h_i|^2 |h_j|^2
% sum_fourth(dd) = 2*U;
% energy_data_correl(:,dd) = alpha./(U^2)*(sum_fourth(dd) + sum_double_square(dd));
% end


% figure;
% plot(e_avg_sym_decod1_correl_b,'-o'); hold on;
% plot(energy_data_correl);
% title('Data energy bob')
% legend('correl simu','no correl simu','correl model','no correl model')
% 
% 
% % 1: Check Bob SINR
% figure;
% plot(sinr1_b_correl_avg,'-o'); hold on;
% plot(energy_data_correl./sigma_b); hold on
% title('Bob SINR')
% legend('correl simu','no correl simu','correl model','no correl model')
% 
% 
% 
% 
% % 3: Check Eve SINR
% 
% % 4: Check each term of Eve SINR
% %
% figure;
% plot(sr1_avg,'-o'); hold on;
% plot(sr1_model) ; hold on;
% title('SR SDS')
% legend('correl simu','no correl simu','correl model','no correl model')

% figure;
% plot(sinr1_model_e); hold on; 
% plot(sinr1_e_avg)
% 
% 
% figure;
% plot(log2(sinr2_model_e)) ; hold on;
% plot(log2(sinr2_e_avg))
% 
% figure; 
% subplot(1,2,1)
% plot(sr1_avg); hold on; 
% plot(sr1_model)
% legend('SDS, simu', 'SDS model')
% 
% subplot(1,2,2)
% plot(sr2_avg); hold on; 
% plot(sr2_model)
% legend('MF, simu', 'MF model')


% figure;
% plot(e_avg_sym_decod1_correl_b); hold on; 
% plot(energy_data_correl)
% legend('Simu', 'Model')
% 

% CHECK QUEL DES 2 TERMES EST PAS BON (sum fourth ou sum double square)

% test1 = squeeze(mean(TEST1,1));
% test2 = squeeze(mean(TEST2,1));
% test3 = squeeze(mean(TEST3,1));
%test4 = squeeze(mean(TEST4));
% 
% figure;
% plot(test1); hold on; plot(sum_double_square)
% legend('sum double square simu', 'sum double square model')
% 
% figure;
% plot(test2); hold on; plot(sum_fourth)
% legend('sum fourth simu', 'sum fourth model')


%% SINR MODELING 
% sum_double = 0;                                                             % ∑_i ∑_j h_i h_j
% for ii = 0: U-2
%     for jj = ii+1 : U-1
%         for kk = 1 : 1+ii*N
%             sum_double = sum_double + 2*real(Tb(1+ii*N,kk) * ctranspose(Tb(1+jj*N,kk)) );
%         end
%     end
% end
% test_sum_double = 1/U*(U+sum_double);
% 
% 
%                                                     
%                  % E[ |S^h |Hb|^2 S|^2 ]
% sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob correl/decorrel


% 
% tmp3 = 0;
% 
% for ii = 0:U-1
%     for kk = 1:1+ii*N
%         tmp3 = tmp3 + 2*abs(Tb(1+ii*N,kk))^4;
%         for ll = 1:1+ii*N
%             if ll ~= kk
%                 tmp3 = tmp3 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+ii*N,ll))^2;
%             end
%         end
%     end
% end
% 
% tmp4 = 0;
% for ii = 0:U-1
%     for kk = 1:1+ii*N
%         kk;
%         %tmp4 = tmp4 + abs(Tb(1+ii*N,kk))^2;
%         for ll = 1:1+ii*N
%             ll;
%            tmp4 = tmp4 + abs(Tb(1+ii*N,ll))^2*abs(Tb(1+ii*N,kk))^2;
%         end
%     end
% end
% tmp4 = 2*tmp4;



% TEST energy AN for MF filter @E

% TERM AN:

% 
%         term_exponent_4 = 0;    % |t_ik|^2 |t_jk|^2 |h_k|^4
%         term_exponent_2 = 0;    % |t_ik|^2 |t_jl|^2 |h_k|^2 |h_l|^2
%         term_indep = 0;         % t_ik t_il^* t_jk t_jl^* |h_k|^2 |h_l|^2
% 
%         for ii = 0: U-1
%             for jj = ii+1 : U-1
%                 for kk = 1 : 1+ii*N
%                      %[{1+ii*N,kk}  {1+jj*N,kk}] % tuples
%                     term_exponent_4 = term_exponent_4 + 2*abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,kk))^2;
%                 end
%             end
%         end
%         for ii = 0: U-2
%             for jj = ii+1 : U-1
%                 for kk = 1 : 1+ii*N
%                     for ll = 1 : 1+jj*N
%                         if ll ~= kk
%                             %[{1+ii*N,kk}  {1+jj*N,ll}] % tuples
%                             term_exponent_2 = term_exponent_2 + abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,ll))^2   ;
%                         end   
%                     end
%                 end
%             end
%         end
% 
%         for ii = 0: U-2
%             for jj = ii+1 : U-1
%                 for kk = 1 : 1+ii*N
%                     for ll = kk+1 : 1+ii*N
%                         %[{1+ii*N,kk}  {1+ii*N,ll} , {1+jj*N,kk} {1+jj*N,ll}] % tuples
%                         term_indep = term_indep + 2*real(conj(Tb(1+ii*N,kk)) * Tb(1+ii*N,ll) * Tb(1+jj*N,kk) * conj(Tb(1+jj*N,ll)))   ;
%                     end
%                 end
%             end
%         end
% 
% 
%         sum_double_square_an = 2*(term_exponent_4 + term_exponent_2 + term_indep);                  % ∑_i ∑_{j!=i}  |h_i|^2 |h_j|^2
%         sum_fourth_an = 2*U;
%         
%         v1_2 = 1./U;
%         v1_4 = 2./U.^2 ./ (1./(U.^2).*(sum_fourth_an + sum_double_square_an));
%         
%          e_an = (1-alpha).*U./(U-1).*(v1_2 - v1_4 )%1./(U+1);%*
%          
%          
%          
%          figure; plot(e_an); hold on; plot(e_avg_an_decod2_e)
%          title('LAST TEST')
% 
%          
%          
% test4 = 1/N/(U-1)*mean(TEST4);
%% IS OK
% % check si TEST2 = sum_double_square
%
% tmp = 0;
% for ii = 0:U-1
%     for jj = 0:U-1
%         if jj ~= ii
%             tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_RX(1+jj*N,1+jj*N))^2;
%         end
%     end
% end
% TEST2(iter,dd) = tmp;
%
%
% % check the value of ∑_i |h_n+iN|^4 = 2*U
%
% tmp2 = 0;
% for ii = 0:U-1
%     tmp2 = tmp2 + abs(Hb_TX(1+ii*N,1+ii*N))^4;
% end
% TEST3(iter,dd) = tmp2;
%
%
% TEST4(iter,dd) = 1/U^2*(TEST2(iter,dd)+TEST3(iter,dd));
% %abs(Hb_TX(1,1)).^2*abs(Hb_RX(3,3)).^2 + abs(Hb_RX(3,3)).^2*abs(Hb_TX(1,1)).^2;

% H2(iter,:,dd) = diag(H_decorrel_TX);


        
        




% % ∑_i ∑_{j!=i}  |Hb_i|^2 |Hb_j|^2  = sum_double_square
% tmp1 = 0;
% for ii = 0:U-1
%     for jj = 0:U-1
%         if jj ~= ii
%             tmp1 = tmp1 + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_TX(1+jj*N,1+jj*N))^2;
%         end
%     end
% end
% TEST1(iter,dd) = tmp1;
% 
% 
% % ∑_i |Hb_i|^4  = 2*U
% tmp2 = 0;
% for ii = 0:U-1
%     tmp2 = tmp2 + abs(Hb_TX(1+ii*N,1+ii*N))^4;
% end
% TEST2(iter,dd) = tmp2; % SHOULD BE = 2U --> OK
% 
% 
% % (test1 + test2)/U^2 = e_avg_sym_decod1_correl_b
% 
% % ∑_i ∑_j > i  |Hb_i|^2 |Hb_j|^2  = TEST1
% tmp3 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         tmp3 = tmp3 + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_RX(1+jj*N,1+jj*N))^2;
%     end
% end
% TEST3(iter,dd) = 2*tmp3; 

% % remplacer H_correl par son expression H_correl,k = ∑_{i=1}^k t_{k,i} h^*_i dans test3 --> = TEST1 OK OK OK 
% tmp4 = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 1:1+ii*N
%             for ll = 1:1+ii*N
%                 for KK = 1:1+jj*N
%                     for LL = 1:1+jj*N
%                         tmp4 = tmp4 + conj(Tb(1+ii*N,kk))*H_decorrel_TX(kk,kk)*Tb(1+ii*N,ll)*H_decorrel_RX(ll,ll)*conj(Tb(1+jj*N,KK))*H_decorrel_TX(KK,KK)*Tb(1+jj*N,LL)*H_decorrel_RX(LL,LL);
%                     end
%                 end
%             end
%         end
%     end
% end
% TEST4(iter,dd) = 2*tmp4; % OK


