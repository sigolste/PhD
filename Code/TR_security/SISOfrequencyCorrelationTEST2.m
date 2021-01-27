%
%
%
%
%   Created 12.01.2021
%   Correct new approx for Bob capa and Eve capa with SDS decoder when  no
%   frequency correlation between Bob's subcarriers. 
%   See SISOfrequencyCorrelationTEST.m to  
%
%
%
%
%
%   © SIDNEY GOLSTEIN

clear all;
% close all;

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
alpha = 1;%0:alpha_step/100:1;         

% Communication parameters
Q = 16;
U = 4;
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters

snr_b  = 10;  % SNR @Bob
snr_e  = 10;  % SNR @Eve

% Noise energy
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve


%% Channel modeling

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance
sigma_tau = .5e-6 ;                                         % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                        % Approximation of coherence bandwidth
coef_freq = [1].*N/6;
delta_f_n = coef_freq.*delta_f_c;   
b_subcar = delta_f_n./N;                                    % Bandwidth of each subcarrier
x_axis  = delta_f_n./delta_f_c;
for dd = 1:length(b_subcar)
    [H(:,:,dd), H1(:,:,dd), abs_rho(:,dd), T(:,:,dd)] = corr_frequency( Q , b_subcar(dd) , sigma_tau , nb_run ) ;
end

%% Energy matrix instantiation
e_noise_e = zeros(nb_run,length(alpha),length(U),length(b_subcar));
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %    
%%                         TEST SECTION                                     %   
%                                                                           %
%                                                                           %
%                                                                           %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              ∑_i |hb|^8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = 0;
for ii = 0:U-1
    tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^8;
end
h8h8(iter,dd) = tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          ∑_i ∑_j~=i |hb_i|^2 |Hb_j|^6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = 0;
for ii = 0:U-1
    for jj = 0:U-1
        if ii~= jj
            tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_TX(1+jj*N,1+jj*N))^6;
        end
    end
end
h2h6(iter,dd) = tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          ∑_i ∑_j~=i |hb_i|^4 |Hb_j|^4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        if ii~= jj
            tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^4*abs(Hb_TX(1+jj*N,1+jj*N))^4;
        end
    end
end
h4h4(iter,dd) = 2*tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ∑_i ∑_j~=i ∑_k ~= i,j |=hb_i|^2 |Hb_j|^2 |Hb_k|^4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (ii~= jj) && (ii~= kk) && (kk ~= jj)
                tmp = tmp + abs(Hb_TX(1+ii*N,1+ii*N))^2*abs(Hb_TX(1+jj*N,1+jj*N))^2*abs(Hb_TX(1+kk*N,1+kk*N))^4;
            end
        end
        
    end
end
h4h2h2(iter,dd) = 2*tmp;


tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for zz = 0:U-1
            for kk = 1:1+ii*N
                tmp1 = tmp1 + abs(T(1+ii*N,kk))^2*abs(T(1+jj*N,kk))^2*abs(Hw_TX(kk,kk))^4*abs(Hb_TX(1+zz*N,1+zz*N))^4;
                for ll = 1:1+jj*N
                    if ll ~= kk
                        tmp2 = tmp2 + abs(T(1+ii*N,kk))^2*abs(T(1+jj*N,ll))^2*abs(Hw_TX(kk,kk))^2*abs(Hw_TX(ll,ll))^2*abs(Hb_TX(1+zz*N,1+zz*N))^4;
                    end
                end
            end
            for kk = 1:1+ii*N
                for ll = kk+1:1+ii*N
                    tmp3 = tmp3 + 2*real(conj(T(1+ii*N,kk))*T(1+jj*N,kk)*T(1+ii*N,ll)*conj(T(1+jj*N,ll)))...
                           *abs(Hw_TX(kk,kk))^2*abs(Hw_TX(ll,ll))^2*abs(Hb_TX(1+zz*N,1+zz*N))^4;
                end
            end
        end
    end
end
h4h2h2_bis(iter,dd) = 2*(tmp1+tmp2+tmp3);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Encoder

% Spreading + TR precoding
sym_spread = matrix_spread*sym_TX;  % Qx1


sym_precoded = Hb_TX*sym_spread;


% AN generation
[an,V1_correl(:,:,iter),S_correl(:,:,iter)] = generateAN_TEST(Hb_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd");% % Qx1, not weighted 


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
[noise_b, e_noise_b(iter,aa,bb,dd) ] = addNoise(sym_b , snr_b, energy(sym_precoded_TX+an_TX));

[noise_e, e_noise_e(iter,aa,bb,dd)  ] = addNoise(sym_e , snr_e, energy(sym_precoded_TX+an_TX));

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

% 
noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;

noise_decod2_e = decod2*noise_e;
% noise_decod3_e = decod3*noise_e;
% noise_decod4_e = decod4*noise_e;
% noise_decod5_e = decod5*noise_e;
% 
an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;

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

capa1_e_avg = squeeze(mean(capa1_e));
% capa2_e_avg = squeeze(mean(capa2_e));

% Ergodic Secrecy Rate
sr1_avg = squeeze(mean(sr1));%capa1_b_correl_avg - capa1_e_avg; 
% sr2_avg = squeeze(mean(sr2));%capa1_b_correl_avg - capa2_e_avg;



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

%% Correlation introduced

% 1. ∑_i |hb_{n+iN}|^8 
H8H8_simu = mean(h8h8);
H8H8_model = modelCorrelH8(U);


% 2. ∑_i ∑_j~=i |hb_{n+iN}|^2 * |hb_{n+jN}|^6 
H2H6_simu = mean(h2h6);
H2H6_model = modelCorrelH2H6(N,U,T);



% 3. ∑_i ∑_j~=i |hb_{n+iN}|^4 * |hb_{n+jN}|^4 


H4H4_simu = mean(h4h4);
H4H4_model  = modelCorrelH4H4(N,U,T);


% 4. ∑_i ∑_j~=i ∑_k~=i,j |hb_{n+iN}|^4 * |hb_{n+jN}|^2 * |hb_{n+kN}|^2
H4H2H2 = mean(h4h2h2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
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






% First term of |hb_i|^4 with |hb_j|^4 --> term4
% tmp = 0;
% for ii = 0:U-1
%     for jj = ii+1:U-1
%         for kk = 1:1+ii*N
%             tmp = tmp + abs(T(1+ii*N,kk))^4*abs(Hw_TX(kk,kk))^4*abs(Hb_TX(1+jj*N,1+jj*N))^4;
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
%                     tmp = tmp + 2*abs(T(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
%                         abs(T(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*abs(Hb_TX(1+jj*N,1+jj*N))^4;
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
% %                         tmp = tmp + 2*abs(T(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
% %                         abs(T(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*abs(T(1+jj*N,mm))^4*abs(Hw_TX(mm,mm))^4;
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
% %                                 tmp = tmp + 2*abs(T(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
% %                                 abs(T(1+ii*N,ll))^2*abs(Hw_TX(ll,ll))^2*...
% %                                 2*abs(T(1+jj*N,mm))^2*abs(Hw_TX(mm,mm))^2*abs(T(1+jj*N,nn))^2*abs(Hw_TX(nn,nn))^2;
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
%                         tmp = tmp + 2*abs(T(1+ii*N,kk))^2*abs(Hw_TX(kk,kk))^2*...
%                         conj(T(1+ii*N,ll))*Hw_RX(ll,ll)*T(1+ii*N,mm)*Hw_TX(mm,mm)*abs(Hb_TX(1+jj*N,1+jj*N))^4;
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
% for dd = 1:length(b_subcar)
% sinr1_correl_model_b(:,dd) = sinrModelingFrequencyCorrelation(alpha,U,N,squeeze(T(:,:,dd)),snr_b,snr_e,"bob_correl");
% 
% sinr1_model_e(:,dd) = sinrModelingFrequencyCorrelation(alpha,U,N,T(:,:,dd),snr_b,snr_e,"eve_decod1_correl");
% sinr2_model_e(:,dd) = sinrModelingFrequencyCorrelation(alpha,U,N,T(:,:,dd),snr_b,snr_e,"eve_decod2_correl");
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
% TT = T; 
% 
% 
% for dd = 1:length(b_subcar)
% T = TT(:,:,dd); 
% term_exponent_4 = 0;
% term_exponent_2 = 0;
% term_exponent_TEST = 0;
% for ii = 0: U-2
%     %term_exponent_TEST = term_exponent_TEST +  2*abs(T(1+ii*N,kk))^4;
%     for jj = ii+1 : U-1
%         %jj
%         %term_exponent_TEST = term_exponent_TEST +  2*abs(T(1+ii*N,kk))^4;
%         for kk = 1 : 1+ii*N
%              %[{1+ii*N,kk}  {1+jj*N,kk}]
%             term_exponent_4 = term_exponent_4 + 2*abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,kk))^2;
%         end
%     end
% end
% for ii = 0: U-2
%     for jj = ii+1 : U-1
%         for kk = 1 : 1+ii*N
%             for ll = 1 : 1+jj*N
%                 if ll ~= kk
%                     %[{1+ii*N,kk}  {1+jj*N,ll}]
%                     term_exponent_2 = term_exponent_2 + abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,ll))^2   ;
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
%                     term_indep = term_indep + 2*real(conj(T(1+ii*N,kk)) * T(1+ii*N,ll) * T(1+jj*N,kk) * conj(T(1+jj*N,ll)))   ;
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
%             sum_double = sum_double + 2*real(T(1+ii*N,kk) * ctranspose(T(1+jj*N,kk)) );
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
%         tmp3 = tmp3 + 2*abs(T(1+ii*N,kk))^4;
%         for ll = 1:1+ii*N
%             if ll ~= kk
%                 tmp3 = tmp3 + 2*abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^2;
%             end
%         end
%     end
% end
% 
% tmp4 = 0;
% for ii = 0:U-1
%     for kk = 1:1+ii*N
%         kk;
%         %tmp4 = tmp4 + abs(T(1+ii*N,kk))^2;
%         for ll = 1:1+ii*N
%             ll;
%            tmp4 = tmp4 + abs(T(1+ii*N,ll))^2*abs(T(1+ii*N,kk))^2;
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
%                     term_exponent_4 = term_exponent_4 + 2*abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,kk))^2;
%                 end
%             end
%         end
%         for ii = 0: U-2
%             for jj = ii+1 : U-1
%                 for kk = 1 : 1+ii*N
%                     for ll = 1 : 1+jj*N
%                         if ll ~= kk
%                             %[{1+ii*N,kk}  {1+jj*N,ll}] % tuples
%                             term_exponent_2 = term_exponent_2 + abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,ll))^2   ;
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
%                         term_indep = term_indep + 2*real(conj(T(1+ii*N,kk)) * T(1+ii*N,ll) * T(1+jj*N,kk) * conj(T(1+jj*N,ll)))   ;
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
%                         tmp4 = tmp4 + conj(T(1+ii*N,kk))*H_decorrel_TX(kk,kk)*T(1+ii*N,ll)*H_decorrel_RX(ll,ll)*conj(T(1+jj*N,KK))*H_decorrel_TX(KK,KK)*T(1+jj*N,LL)*H_decorrel_RX(LL,LL);
%                     end
%                 end
%             end
%         end
%     end
% end
% TEST4(iter,dd) = 2*tmp4; % OK


