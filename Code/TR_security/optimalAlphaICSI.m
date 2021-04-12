function alpha_opt = optimalAlphaICSI(U,snr_b,snr_e,sigma_tilde,type)
% Compute the percentage of data energy to inject in the communication to
% maximize the SR, depending on the scenario investigated.
%
% INPUTS:
%   U: back-off rate
%   snr_b: signal to noise ratio @Bob
%   snr_e: signal to noise ratio @Eve
%   sigma_tilde : error of CSI estimation @A (between 0 and 1)
%   type: investigatted scenario
%       decod1: Bob and Eve with the same decoding structure
%       decod5: Eve only knows her own channel He
%
% OUTPUT:
%   alpha_opt: optimal amount of data energy to radiate to maximize the SR
%
% By Sidney Golstein, Mar. 2021.


sigma_b = 1./U.'./10.^(snr_b./10);    % expected noise energy @Bob
sigma_e = 1./U.'./10.^(snr_e./10);    % expected noise energy @Eveswitch type
switch type
    case "decod1"
        % Cfr derivation verso ICSI5
        T1 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        T2 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        T3 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        T4 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        for bb = 1:length(U)
        for ss = 1:length(sigma_tilde)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            T1(bb,ss,nb,ne) = (1-sigma_tilde(ss))*(U(bb)+1);
            T2(bb,ss,nb,ne) = (U(bb)*sigma_e(bb,ne)+1)*(1-sigma_tilde(ss))*(U(bb)+1) - (U(bb)*sigma_b(bb,nb) + sigma_tilde(ss));
            T3(bb,ss,nb,ne) = (U(bb)*sigma_b(bb,nb)+sigma_tilde(ss))*(U(bb)*sigma_e(bb,ne)+1);
            T4(bb,ss,nb,ne) = -sigma_tilde(ss)*(U(bb)*sigma_e(bb,ne)+1);
        end
        end
        end
        end
        alpha_opt = (sqrt(T1.^2.*T3.^2 + T1.*T2.*T3.*T4 - T1.*T3.*T4.^2) - T1.*T3 )./(T1.*T4);
        alpha_opt(:,1,:,:) = T2(:,1,:,:)./2./T1(:,1,:,:);      % when sigma_tilde = 0, the SR polynom changes, only second order in alpha divided by positive constant.
        
        for bb = 1:length(U)
        for ss = 1:length(sigma_tilde)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            if real(alpha_opt(bb,ss,nb,ne)) > 1
                alpha_opt(bb,ss,nb,ne) = NaN;
            end
        end
        end
        end
        end
        
    case "decod5"
        % Cfr derivation verso ICSI5
        T1 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        T2 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        T3 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        T4 = zeros(length(U),length(sigma_tilde),length(snr_b),length(snr_e));
        for bb = 1:length(U)
        for ss = 1:length(sigma_tilde)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            T1(bb,ss,nb,ne) = (1-sigma_tilde(ss))*(U(bb)+1);
            T2(bb,ss,nb,ne) = (U(bb)*sigma_e(bb,ne)/2+1)*(1-sigma_tilde(ss))*(U(bb)+1) - (U(bb)*sigma_b(bb,nb) + sigma_tilde(ss));
            T3(bb,ss,nb,ne) = (U(bb)*sigma_b(bb,nb)+sigma_tilde(ss))*(U(bb)*sigma_e(bb,ne)/2+1);
            T4(bb,ss,nb,ne) = -sigma_tilde(ss)*(U(bb)*sigma_e(bb,ne)/2+1);
        end
        end
        end
        end
        alpha_opt = (sqrt(T1.^2.*T3.^2 + T1.*T2.*T3.*T4 - T1.*T3.*T4.^2) - T1.*T3 )./(T1.*T4);
        alpha_opt(:,1,:,:) = T2(:,1,:,:)./2./T1(:,1,:,:);      % when sigma_tilde = 0, the SR polynom changes, only second order in alpha divided by positive constant.
        
        for bb = 1:length(U)
        for ss = 1:length(sigma_tilde)
        for nb = 1:length(snr_b)
        for ne = 1:length(snr_e)
            if real(alpha_opt(bb,ss,nb,ne)) > 1
                alpha_opt(bb,ss,nb,ne) = NaN;
            end
        end
        end
        end
        end
end
