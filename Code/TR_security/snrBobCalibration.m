function [snr_b_calibrated , alpha_calibrated] = snrBobCalibration(sr_targetted,U,snr_e,alpha,alpha_step,type)
% Function that returns the SNR that is needed at Bob to reach a targetted
% secrecy rate (SR) as a function of the back of rate (U), the SNR at Eve 
% (snr_e) and alpha. 
% For given targetted SR and a given set of (snr_e,U), we find the value of
% alpha that gives the lowest value of snr_b (i.e. the easiest to obtain in
% reality). We will finally inject this value of alpha in our system since
% it allows to have the lowest SNR at Bob while reaching a given SR for a
% given couple of (U,snr_e)
%
%
% INPUTS:
%   sr_targetted : the targetted secrecy rate
%   U : back of rate
%   snr_e : matrix with the different value of Eve's snr
%   alpha : matrix with the different amounts of energy radiated and
%           dedicated for the data
%   alpha_step : 
%   type : the model investigated
%
% OUTPUT:
%   snr_b : Bob's snr needed to reach sr_targetted for a tuple
%   (U,snr_e,alpha)
%
% By Sidney Golstein
% Last update: 10.07.2020

sigma_e = 1./U./(10.^(snr_e.'./10)); % dim: [length(snr_e), length(U)]
sr_lin = 2^sr_targetted; % linear value of the targetted secrecy rate


snr_b_tmp = zeros(length(sr_targetted), length(snr_e), length(U), length(alpha));
% snr_b = zeros(length(sr_targetted), length(snr_e), length(U));


switch type
    case "model1"

for ss = 1:length(sr_targetted)
    for nn = 1:length(snr_e)
        for uu = 1:length(U)
            for aa = 1:length(alpha)
                snr_b_tmp(ss,nn,uu,aa) = 10*log10( (sr_lin * (U(uu) * sigma_e(nn,uu) + 1 ) - (U(uu) * sigma_e(nn,uu) + 1 - alpha(aa) ) ) / (alpha(aa)*(U(uu)+1)*(U(uu) * sigma_e(nn,uu) + 1 - alpha(aa)) ) ); 
            end
        end
    end
end
    case "model2"
for ss = 1:length(sr_targetted)
    for nn = 1:length(snr_e)
        for uu = 1:length(U)
            for aa = 1:length(alpha)
                snr_b_tmp(ss,nn,uu,aa) = 10*log10( (sr_lin *( U(uu)*(U(uu)+1)* sigma_e(nn,uu)  + alpha(aa)*(U(uu)+1)*(U(uu)+3) + U(uu)*(1-alpha(aa)) ) - (U(uu)*(U(uu)+1)* sigma_e(nn,uu) + (1-alpha(aa))*U(uu) ) ) ... 
                    / (U(uu)*alpha(aa)*(U(uu)+1) * ((U(uu)+1)*sigma_e(nn,uu) + (1-alpha(aa)) ))    );  % Aaaaaaaaaaaaaaaaa
            end 
        end
    end
end        
        
end

[snr_b_calibrated , idx] = min(snr_b_tmp,[],4);
snr_b_calibrated = squeeze(snr_b_calibrated);
alpha_calibrated = squeeze(idx).* alpha_step;

%% Plot section

figure;

sgtitle(['Targetted secrecy capacity : ', num2str(sr_targetted), ' bit(s) per channel use'])

subplot(1,2,1)
surf(snr_e,U,alpha_calibrated.'-1)
%zlim([0 100])
xlabel('Eve SINR (dB)')
ylabel('BOR')
zlabel('$\alpha$ to inject ($\%$)')

subplot(1,2,2)
surf(snr_e,U,snr_b_calibrated.')
xlabel('Eve SINR (dB)')
ylabel('BOR')
zlabel('Bob SINR (dB)')
%title(['Targetted secrecy capacity : ', num2str(sr_targetted), ' bit(s) per channel use'])


%title(['Targetted secrecy capacity : ', num2str(sr_targetted), ' bit(s) per channel use'])


test = 1;


