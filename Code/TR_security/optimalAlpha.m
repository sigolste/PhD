function alpha_opt = optimalAlpha(U,snr_b,snr_e,type)
% Compute the percentage of data energy to inject in the communication to
% maximize the SR, depending on the scenario investigated.
%
% INPUTS:
%   U: back-off rate
%   snr_b: signal to noise ratio @Bob
%   snr_e: signal to noise ratio @Eve
%   type: investigatted scenario
%       model1: Bob and Eve with the same decoding structure
%       model2: Matched filtering @Eve
%       model5: Eve only knows her own channel He
%
% OUTPUT:
%   alpha_opt: optimal amount of data energy to radiate to maximize the SR
%
% By Sidney Golstein, May 2020.


sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Bob
switch type
    case "model1"
        a = -(U+1);
        b = (U+1).*U.*sigma_e + (U+1) - U.*sigma_b;
        alpha_opt = -b./(2.*a);
        
    case "model2"
        T1 = U+1;
        T2 = (U+1).^2.*sigma_e+(U+1)-U.*sigma_b;
        T3 = U.*(U+1).*sigma_b.*sigma_e+U.*sigma_b;
        T4 = (U+1).*(U+3).*sigma_b-U.*sigma_b;
        alpha_opt = (sqrt(T1.^2.*T3.^2+T1.*T2.*T3.*T4-T1.*T3.*T4.^2)-T1.*T3)./(T1.*T4);
        
    case "model5"
        a = -2*(U+1);
        b = (U+1).*(2+U.*sigma_e) - 2.*U.*sigma_b;
        alpha_opt = -b./(2.*a);
end
