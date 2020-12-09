function sinr = sinrModelingICSI(alpha,U,snr_b,snr_e,sigma,type)
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Bob
% sigma_an = 1./U;

alpha = alpha.';
switch type 
    case "bob_decod1"
        num = alpha/U*((U+1)*(1-sigma)+sigma);
        an = (1-alpha)*sigma/U;
        sinr = num./(an + sigma_b);
    case "eve_decod1"
        sinr = alpha./(U.*sigma_e + (1-alpha));
        
    case "eve_decod2"
        % Consideration: normalisation of decod 2 by sqrt((U+1)/(U+3)). 
        % Else, multiply all by (U+3)/(U+1). Also, we consider that energy 
        % of an signal is 1/U
        e_sym_decod2_model_e = alpha.*(U+1)./U;
        e_noise_decod2_model_e = (U+1)./(U+3).*sigma_e;
        e_an_decod2_model_e = (1-alpha)./(U+3);
        sinr = e_sym_decod2_model_e./(e_noise_decod2_model_e + e_an_decod2_model_e);
        
     case "eve_decod5"
        % Consideration: normalisation of decod 2 by sqrt((U+1)/(U+3)). 
        % Else, multiply all by (U+3)/(U+1). Also, we consider that energy 
        % of an signal is 1/U
        e_sym_decod5_model_e = 2*alpha./U;
        e_noise_decod5_model_e =sigma_e;
        e_an_decod5_model_e = 2./U.*(1-alpha);
        sinr = e_sym_decod5_model_e./(e_noise_decod5_model_e + e_an_decod5_model_e);
end