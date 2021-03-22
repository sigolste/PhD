function sinr = sinrModelingICSI(alpha,U,snr_b,snr_e,sigma,type)
sigma_b = 1./U.'./10.^(snr_b./10);    % expected noise energy @Bob
sigma_e = 1./U.'./10.^(snr_e./10);    % expected noise energy @Eve
% sigma_an = 1./U;
U = U.';
alpha = alpha.';
switch type 
    case "bob_decod1"
        % cfr ICSI4 verso
        e_num = alpha./U.*((U+1)*(1-sigma)+sigma);
        e_an = (1-alpha).*sigma./U;
        e_noise = sigma_b;
        sinr = e_num./(e_an + e_noise);
    case "eve_decod1"
        sinr = alpha./(U.*sigma_e + (1-alpha));
        
    case "eve_decod2"
        e_num = alpha.*(U+1)./U;
        e_noise = (U+1)./(U+3).*sigma_e;
        e_an = (1-alpha)./(U+3);
        sinr = e_num./(e_noise + e_an);
        
     case "eve_decod5"
        % Consideration: normalisation of decod 2 by sqrt((U+1)/(U+3)). 
        % Else, multiply all by (U+3)/(U+1). Also, we consider that energy 
        % of an signal is 1/U
        e_num = 2*alpha./U;
        e_noise = sigma_e;
        e_an = 2./U.*(1-alpha);
        sinr = e_num./(e_noise + e_an);
end