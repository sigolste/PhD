function sigma_max = ICSImaxSigmaEveNoiseless(U,SR,type)
% Compute the maximal allowed ICSI error that can be made in order to
% ensure to target a desired SR = ∆ bit/channel use. In this
% scenario, Eve is noiseless, i.e., worst case scenario in terms of securty 
% 
%
% INPUTS:
%   U: back-off rate
%   SR : targetted Secrecy Rate (linear scale)
%   type: investigatted scenario
%       decod1: Bob and Eve with the same decoding structure
%       decod2: MF decoder at Eve
%       decod5: Eve only knows her own channel He
%
% OUTPUT:
%   sigma_max : maximum allowed channel state information estimaton error 
%   made @Alice when Eve is noiseless
%   Remark: sigma_max \in [0,1]
%   
%
% © Ir. Sidney Golstein, 26 Mar. 2021.


switch type
    case "decod1"
        % Cfr derivation recto ICSI17
        sigma_max = 1 - (SR - 1)./((U+1) + (SR-1));
    case "decod2"
        % Cfr derivation verso ICSI21
        B = U.^2 +3.*U+3;
        sigma_max = 1 - (SR.*B)./(SR.*B + U.*(U+1));


    case "decod5"
        % Cfr derivation recto ICSI17
        sigma_max = 1 - (SR - 1)./((U+1) + (SR-1));

end
