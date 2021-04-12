function sigma_max = maxSigmaEveNoisy(U,snr_b,snr_e,alpha,type)
% Compute the maximal allowed ICSI error that can be made in order to
% ensure a positive secrecy rate, i.e. \gamma_b > \gamma_e. In this
% scenario, Eve is not noiseless, i.e., not the worst case scenario in
% terms of securty. 
% See maxSigmaNoiseless.m in order to find the condition on sigma that can
% be made at Alice to ensure a targetted SR when Eve is AWGN free, i.e.
% worst case scenario. 
% 
%
% INPUTS:
%   U: back-off rate
%   snr_b: signal to noise ratio @Bob
%   snr_e: signal to noise ratio @Eve
%   alpha : Percentage of data energy injected
%   type: investigatted scenario
%       decod1: Bob and Eve with the same decoding structure
%       decod2: MF decoder at Eve
%       decod5: Eve only knows her own channel He
%
% OUTPUT:
%   sigma_max : maximum allowed channel state information estimaton error 
%   made @Alice.
%   Remark: sigma_max \in [0,1]
%   
%
% © Ir. Sidney Golstein, Mar. 2021.


sigma_b = 1./U.'./10.^(snr_b./10);    % expected noise energy @Bob
sigma_e = 1./U.'./10.^(snr_e./10);    % expected noise energy @Eveswitch type
sigma_max = zeros(length(alpha),length(U),length(snr_b),length(snr_e));
switch type
    case "decod1"
        % Cfr derivation recto ICSI21
        for aa = 1:length(alpha)
        for bb = 1:length(U)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            sigma_max(aa,bb,nb,ne) = 1 - ( U(bb)*(sigma_b(bb,nb) - sigma_e(bb,ne)) ) ...
                                          /(U(bb)*(U(bb)*sigma_e(bb,ne) + 1- alpha(aa) ) + 1-alpha(aa) );
            if (sigma_max(aa,bb,nb,ne) > 1) % No condition on sigma.
                sigma_max(aa,bb,nb,ne) = 1;
            end
            if (sigma_max(aa,bb,nb,ne) < 0)
                sigma_max(aa,bb,nb,ne) = NaN;
            end
        end
        end
        end
        end
    case "decod2"
        % Cfr derivation recto ICSI21
        for aa = 1:length(alpha)
        for bb = 1:length(U)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            sigma_max(aa,bb,nb,ne) = ( U(bb)*(U(bb)+1)*( ( (U(bb)+1)*sigma_e(bb,ne) + 1-alpha(aa) ) - (U(bb)+3)*sigma_b(bb,nb) ) ) ...
                                          /( (1-alpha(aa))*(U(bb)+1)*(U(bb)+3) + U(bb)^2*((U(bb)+1)*sigma_e(bb,ne) + 1-alpha(aa) ) );
            if (sigma_max(aa,bb,nb,ne) > 1) 
                sigma_max(aa,bb,nb,ne) = 1;
            end
            if (sigma_max(aa,bb,nb,ne) < 0)
                sigma_max(aa,bb,nb,ne) = NaN;
            end
        end
        end
        end
        end

    case "decod5"
        % Cfr derivation recto ICSI21
        for aa = 1:length(alpha)
        for bb = 1:length(U)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            sigma_max(aa,bb,nb,ne) = 1 - ( U(bb)*(sigma_b(bb,nb) - sigma_e(bb,ne)/2) ) ...
                                          /(U(bb)*(U(bb)*sigma_e(bb,ne)/2 + 1- alpha(aa) ) + 1-alpha(aa) );
            if (sigma_max(aa,bb,nb,ne) > 1) 
                sigma_max(aa,bb,nb,ne) = 1;
            end
            if (sigma_max(aa,bb,nb,ne) < 0)
                sigma_max(aa,bb,nb,ne) = NaN;
            end
        end
        end
        end
        end

end
