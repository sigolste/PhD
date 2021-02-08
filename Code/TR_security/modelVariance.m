function var = modelVariance(alpha,N,U,T,sinr,snr_b,snr_e,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of variance of Bob/Eve SINR when frequency correlation is 
% introduced.
% This term is needed in order to better approximate the capacity at Bob/
% Eve when frequency correlation among subcarriers is introduced.
%
% var = E[X^2] - (E[X])^2 where X is Bob/Eve SINR
%
% INPUTS:
%   alpha : percentage of data energy injected in the communication
%   N : Number of symbol sent (N = Q/U)
%   U : Back-off rate
%   T : Choleski decomposition of correlation matrix RHO, s.t. T T^H = RHO
%   sinr : Bob or ve SINR depending on the model
%   snr_b: Bob's SNR
%   snr_e: Eve's SNR
%   type : Scenario that has to be modelized
%
% OUTPUT:
%   var : The modelized variance of Bob/Eve SINR
%
%
%
% Code started : 27.01.2021
% Last update  : 04.02.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preliminary:
% Approximation of E[1+ log(X)] = log(1 + E[X]) - var(X)/(2( 1+E[X] )^2)
% where X = SINR_b ou SINR_e
sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
alpha  = alpha.';
switch type
    case "bob_noCorrel"
    case "bob_correl"
        sinr_square     = sinr.^2;
        H8_model        = modelCorrelH8(U);
        H2H6_model      = modelCorrelH2H6(N,U,T);
        H4H4_model      = modelCorrelH4H4(N,U,T);
        H4H2H2_model    = modelCorrelH4H2H2(N,U,T);
        H2H2H2H2_model  = modelCorrelH2H2H2H2(N,U,T);
        e_sym_square    = alpha.^2./U.^4.*(H8_model + 4*H2H6_model + ...
                          3*H4H4_model + 6*H4H2H2_model + H2H2H2H2_model);
        e_noise_square  = 2*sigma_b^2;
        var = e_sym_square./e_noise_square - sinr_square;


    case "eve_decod1_noCorrel"
    case "eve_decod2_noCorrel"
    case "eve_decod5_noCorrel"
    case "eve_decod1_correl"
        sinr_square     = sinr.^2;
        e_denom_square  = 1./U*(2*sigma_e^2 + 4./U^2.*(1-alpha).^2 + 2*(U-1)*sigma_e^2 ...
                + 2./U^2.*(1-alpha).^2.*(U-1)+4./U.*(1-alpha).*sigma_e ...
                + 8./U.*(1-alpha).*(U-1).*sigma_e);     % e_denom_square computed verso FC33                                                                         
        H4_model        = modelCorrelH4(U);
        H2H2_model      = modelCorrelH2H2(N,U,T);
        e_sym_square    = 2.*alpha.^2./U.^4*(H4_model+H2H2_model);
        var = e_sym_square./e_denom_square - sinr_square;

    case "eve_decod2_correl"
    case "eve_decod5_correl"
end