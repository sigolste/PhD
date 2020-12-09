sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
T1 = U-U*sigma_tilde-sigma_tilde + 1;
T2 = (U*sigma_e + 1)*(U-U*sigma_tilde-sigma_tilde+1)-(U*sigma_b+sigma_tilde);
T3 = (U*sigma_b + sigma_tilde)*(U*sigma_e+1);
T4 = -(U*sigma_e*sigma_tilde+sigma_tilde);
alpha_opt = (sqrt(T1.^2.*T3.^2+T1.*T2.*T3.*T4-T1.*T3.*T4.^2)-T1.*T3)./(T1.*T4)


% sigma_tilde = 0
A = (1-sigma_tilde)*(U+1)*U;
B = (1-sigma_tilde)*(U+1)*U*((U+1)*sigma_e+1)-U*(U*sigma_b+sigma_tilde);
C = (U*sigma_b+sigma_tilde)*U*((U+1)*sigma_e+1);
D = sigma_tilde*(U^2+3*U+3);
E = (U*sigma_b+sigma_tilde)*(U^2+3*U+3)-sigma_tilde*U*((U+1)*sigma_e+1);

% sigma_e = [0:.02:10];
% sigma_e =  (U*sigma_b*(U^2+3*U+3)+U^2+2*U+3)/(U*(U+1))*20;
test = log2((-A.*alpha.^2.'+B.*alpha.' + C)./(-D.*alpha.^2.'+E.*alpha.'+C))
sigmaT2 = 1-(U*sigma_b+1)/((U+1)^2*sigma_e+U+2)
sigmaT5 = (U*sigma_b*(U^2+3*U+3))./(-U^2-2*U-3+(U+1)*U.*sigma_e)

sol = (2*C.*(A-D)-sqrt(4*C.^2.*(A-D).^2-4.*B.*D.*(B.*C-C.^2)))./(2.*B.*D)


e_sym_decod2_model_e = alpha.*(U+1)./U;
        e_noise_decod2_model_e = (U+1)./(U+3).*sigma_e;
        e_an_decod2_model_e = (1-alpha)./(U+3);
        sinr_e = e_sym_decod2_model_e./(e_noise_decod2_model_e + e_an_decod2_model_e);
        
        lea = 1+sinr_e
        (alpha.*(U+1)-U*sigma_b.*lea)./((1-alpha).*lea+alpha.*U)
  
        

AA = U^2*sigma_b*((U^2+3*U+3)*((U+1)*sigma_e+1));
BB = 