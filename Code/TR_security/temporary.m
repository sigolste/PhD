
%% SDS Decoder
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve

% 1 Best alpha
T1 = U-U*sigma_tilde-sigma_tilde + 1;
T2 = (U*sigma_e + 1)*(U-U*sigma_tilde-sigma_tilde+1)-(U*sigma_b+sigma_tilde);
T3 = (U*sigma_b + sigma_tilde)*(U*sigma_e+1);
T4 = -(U*sigma_e*sigma_tilde+sigma_tilde);
alpha_opt = (sqrt(T1.^2.*T3.^2+T1.*T2.*T3.*T4-T1.*T3.*T4.^2)-T1.*T3)./(T1.*T4);

% 2 Max of CSI error still allowing SINR_b > SINR_e
U = 4;
snr_b = [-7.1600]; 
snr_e = [20]; 
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
alpha = [0:0.05:1];

sigma_tilde = 0.1;
sigma_max_decod1 = 1 +U*(sigma_e-sigma_b)./(U*(U*sigma_e+1-alpha)+(1-alpha)); % ICSI 15 rectp
sr1_ICSI_model = (-alpha.^2.*(U+1)*(1-sigma_tilde) + alpha.*((U*sigma_e+1)*(U+1)*(1-sigma_tilde)-(U*sigma_b+sigma_tilde)) + (U*sigma_b+sigma_tilde)*(U*sigma_e+sigma_tilde))./(-alpha.*sigma_tilde*(U*sigma_e+1)+(U*sigma_b+sigma_tilde)*(U*sigma_e+sigma_tilde));

% 3 SNR Formula
alpha = [0:0.05:1].';
snr_e = [20] ;
sigma_e = 1./U/10^(snr_e/10); 

delta = 16; % 2^SR
sigma_tilde = [0:0.01:0.9];

A = U*sigma_e + 1;

% SNR when sigma_e \neQ 0 ICSI P16 recto
snr_b_decod1 = 10*log10((alpha+A*(delta-1))./(-alpha.^2.*(U+1).*(1-sigma_tilde) + alpha.*(A.*((U+1).*(1-sigma_tilde) + delta.*sigma_tilde) - sigma_tilde) + sigma_tilde.*A.*(1-delta))); % OK

 
% % SNR when sigma_e = 0 ICSI P16 recto
% alpha = [0:0.05:1].';
% 
% delta = 16; % 2^SR
% sigma_tilde = [0:0.01:0.9];
% snr_b_noEveNoise_decod1 = 10*log10((alpha+(delta-1))./(-alpha.^2.*(U+1).*(1-sigma_tilde) + alpha.*((U+1).*(1-sigma_tilde) + delta.*sigma_tilde - sigma_tilde) + sigma_tilde.*(1-delta))); % OK


% 4 Guaranteeing SR:
% 4.1 SNR when sigma_e = 0
alpha = [0:0.05:1].';

delta = 2^1; % 2^SR
sigma_tilde = [0:0.01:0.9];
snr_b_noEveNoise_decod1 = 10*log10((alpha+(delta-1))./(-alpha.^2.*(U+1).*(1-sigma_tilde) + alpha.*((U+1).*(1-sigma_tilde) + delta.*sigma_tilde - sigma_tilde) + sigma_tilde.*(1-delta))); % OK ICSI P16 recto

% 4.2 Condition on sigma
sigma_tilde_lim_decod1 = 1 - (delta-1)./((U+1)+(delta-1));      % Ok: max error acceptable to guarantee SR= ∆ when sigma_e = 0 ICSI p17 recto

% 4.3 Best alpha
delta = 2; % 2^SR
sigma_tilde = [0:0.01:0.9];
A = (U+1).*(1-sigma_tilde);
B = sigma_tilde.*(delta-1);
C = A+B;
d_t = 4.*A.^2.*(delta-1).^2 - 4.*A.*(-(delta-1).*C-B);
alpha_opt_noEveNoise_decod1 = (-2.*A.*(delta-1)  + sqrt(d_t) ) ./ (2.*A); % OK : BEST ALPHA WHEN SIGMA_E = 0 




%% MF Decoder
% sigma_tilde = 0
snr_e = [20] ;
snr_b = [20] ;
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
U = 4; 
alpha = [0:0.05:1].';
sigma_tilde = [0:0.05:0.3];

A = (1-sigma_tilde)*(U+1)*U;
B = (1-sigma_tilde)*(U+1)*U*((U+1)*sigma_e+1)-U*(U*sigma_b+sigma_tilde);
C = (U*sigma_b+sigma_tilde)*U*((U+1)*sigma_e+1);
D = sigma_tilde*(U^2+3*U+3);
E = (U*sigma_b+sigma_tilde)*(U^2+3*U+3)-sigma_tilde*U*((U+1)*sigma_e+1);

% sigma_e = [0:.02:10];
% sigma_e =  (U*sigma_b*(U^2+3*U+3)+U^2+2*U+3)/(U*(U+1))*20;
sr2_ICSI_model = log2((-A.*alpha.^2 +B.*alpha + C)./(-D.*alpha.^2+E.*alpha+C));
% 
% sigmaT2 = 1-(U*sigma_b+1)/((U+1)^2*sigma_e+U+2);
% sigmaT5 = (U*sigma_b*(U^2+3*U+3))./(-U^2-2*U-3+(U+1)*U.*sigma_e);
% 
% sol = (2*C.*(A-D)-sqrt(4*C.^2.*(A-D).^2-4.*B.*D.*(B.*C-C.^2)))./(2.*B.*D);
% 
e_sym_decod2_model_e = alpha.*(U+1)./U;
e_noise_decod2_model_e = (U+1)./(U+3).*sigma_e;
e_an_decod2_model_e = (1-alpha)./(U+3);
sinr_e = e_sym_decod2_model_e./(e_noise_decod2_model_e + e_an_decod2_model_e);
capa_e = 1+sinr_e;

% Max acceptable sigma_tilde to have positive SR as a function of alpha
sigma_max_decod2 = (alpha.*(U+1)-U*sigma_b.*sinr_e)./((1-alpha).*sinr_e+alpha.*U) ;
  
% SNR bob when SNR eve = infty     
% Test sign of term 4
U = 4';
delta =0;
A = U.*((U+1).*sigma_e+1);
B = U.^2+3.*U+3;
sigma_term4 = (A.*(U+1))./(2.^delta.*B+(U+1).*A+U - 2.^delta.*A);           % term 4 = 0 if sigma = this expression ---> OK
% Test if SNR b corresponds to the formula verso p12


% Guaranteeing SR:
U = 64;
snr_e = 1000;                     % Infinite SNR @E
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
delta = 1;

alpha = [0:0.01:1]';
A = U.*((U+1).*sigma_e+1);
B = U.^2+3.*U+3;
sigma_tilde = [0:0.01:0.4];

num_snr_b_noEveNoise = (alpha.*(U+2^delta*B)+A.*(2^delta-1)); 


t1 = 2^delta*B.*sigma_tilde-U*(U+1)*(1-sigma_tilde);
t2 = 2^delta*A.*sigma_tilde-sigma_tilde.*2^delta*B+(U+1).*(1-sigma_tilde).*A-sigma_tilde*U;
t3 = sigma_tilde.*A*(1-2^delta);
denom_snr_b_noEveNoise = (alpha.^2*(t1) + alpha.*(t2) + t3);


snr_test_b = (alpha.*(U+2^delta*B)+A*(2^delta-1))./(alpha.^2*(2^delta*B.*sigma_tilde-U*(U+1)*(1-sigma_tilde)) + alpha.*(2^delta*A.*sigma_tilde-sigma_tilde.*2^delta*B+(U+1).*(1-sigma_tilde).*A-sigma_tilde*U) + sigma_tilde.*A*(1-2^delta)); 

snr_b_noEveNoise_decod2 = 10*log10(snr_test_b);

% Formula max sigma acceptable to ensure SR = ∆ when sigma_e = 0 (verso
% ICSI p13)
U = [2 4 8 16].';
B = U.^2+3.*U+3;
delta = [0:0.1:3.5];
sigma_tilde_lim_decod2 = (U.*(U+1))./(2.^delta.*B+U.*(U+1));

sigma_tilde = 0.04;
% Formula optimal alpha to ensure SR = ∆ : test (resto ICSI p14)
X = sigma_tilde.*(2.^delta.*B+U.*(U+1)) - U.*(U+1);
Y = sigma_tilde.*(2.^delta.*U-2.^delta.*B-U.*(U+2))+U.*(U+1);
Z = U+2.^delta.*B;
T = U*(1-2.^delta);

alpha_opt_noEveNoise_decod2 = (X.*T-sqrt(X.^2.*T.^2+Z.*X.*T.*(sigma_tilde.*Z+Y)))./(Z.*X);
A = U;
alpha = alpha_opt_noEveNoise_decod2;
snr_opt_NoEveNoise_decod2_b = 10*log10((alpha.*(U+2.^delta.*B)+A.*(2.^delta-1))./(alpha.^2.*(2.^delta.*B.*sigma_tilde-U.*(U+1).*(1-sigma_tilde)) + alpha.*(2.^delta.*A.*sigma_tilde-sigma_tilde.*2.^delta.*B+(U+1).*(1-sigma_tilde).*A-sigma_tilde.*U) + sigma_tilde.*A.*(1-2.^delta))).'; 
figure;
plot(delta,snr_opt_NoEveNoise_decod2_b,'-o');box on; grid on;
xlabel('Targetted SR  $\Delta$ (bit/channel use)')
ylabel('Required Bob SNR (dB)')
legend('U = 2', 'U = 4' , 'U = 8' , 'U = 16')

%% OC Decoder 

sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
sigma_tilde = 0.80;
% 1 Best alpha
T1 = U-U*sigma_tilde-sigma_tilde + 1;
T2 = (U*sigma_e/2 + 1)*(U-U*sigma_tilde-sigma_tilde+1)-(U*sigma_b+sigma_tilde);
T3 = (U*sigma_b + sigma_tilde)*(U*sigma_e/2+1);
T4 = -(U*sigma_e/2*sigma_tilde+sigma_tilde);
alpha_opt_decod5 = (sqrt(T1.^2.*T3.^2+T1.*T2.*T3.*T4-T1.*T3.*T4.^2)-T1.*T3)./(T1.*T4);

% Max sigma to have SNRb > SNRe
alpha = [0:0.05:1];
sigma_max_decod5 = 1 +U*(sigma_e/2-sigma_b)./(U*(U*sigma_e/2+1-alpha)+(1-alpha)); % ICSI 15 recto to have SNRb > SNRe as a function of alpha


% 3 SNR Formula
alpha = [0:0.05:1].';
snr_e = [20] ;
sigma_e = 1./U/10^(snr_e/10); 

delta = 16; % 2^SR
sigma_tilde = [0:0.01:0.9];

A = U*sigma_e/2 + 1;

% SNR when sigma_e \neQ 0 ICSI P16 recto
snr_b_decod5 = 10*log10((alpha+A*(delta-1))./(-alpha.^2.*(U+1).*(1-sigma_tilde) + alpha.*(A.*((U+1).*(1-sigma_tilde) + delta.*sigma_tilde) - sigma_tilde) + sigma_tilde.*A.*(1-delta))); % OK

 
% % SNR when sigma_e = 0 ICSI P16 recto
% alpha = [0:0.05:1].';
% 
% delta = 16; % 2^SR
% sigma_tilde = [0:0.01:0.9];
% snr_b_noEveNoise_decod5 = 10*log10((alpha+(delta-1))./(-alpha.^2.*(U+1).*(1-sigma_tilde) + alpha.*((U+1).*(1-sigma_tilde) + delta.*sigma_tilde - sigma_tilde) + sigma_tilde.*(1-delta))); % OK
% 

% 4 Guaranteeing SR:
% 4.1 SNR when sigma_e = 0
alpha = [0:0.05:1].';

delta = 2^4; % 2^SR
sigma_tilde = [0:0.01:0.9];
snr_b_noEveBoise_decod5 = 10*log10((alpha+(delta-1))./(-alpha.^2.*(U+1).*(1-sigma_tilde) + alpha.*((U+1).*(1-sigma_tilde) + delta.*sigma_tilde - sigma_tilde) + sigma_tilde.*(1-delta))); % OK idem que decoder 1

% 4.2 Condition on sigma
sigma_tilde_lim_decod5 = 1 - (delta-1)./((U+1)+(delta-1));      % Ok: max error acceptable to guarantee SR= ∆ when sigma_e = 0 ICSI p17 recto

% 4.3 Best alpha
delta = 16; % 2^SR
sigma_tilde = [0:0.01:0.9];
 A = (U+1).*(1-sigma_tilde);
 B = sigma_tilde.*(delta-1);
 C = A+B;
 d_t = 4.*A.^2.*(delta-1).^2 - 4.*A.*(-(delta-1).*C-B);
 alpha_opt_decod5 = (-2.*A.*(delta-1)  + sqrt(d_t) ) ./ (2.*A); % OK : BEST ALPHA WHEN SIGMA_E = 0 idem que decod1 
 
 
 


% AA = U^2*sigma_b*((U^2+3*U+3)*((U+1)*sigma_e+1));
% BB = 