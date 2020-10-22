% Script that returns the SNR that is needed at Bob, as well as the AN to
% inject in order to reach a targetted secrecy rate (SR) as a function of 
% the back of rate (U), when the SNR of Eve is infinite. 
%
% Therefore, it allows to guarantee a communication SR whatever Eve's
% condition (SNR_E = infinite) under the hypothesis of:
% - SISO communication
% - No spatial nor frequency correlation
% - 3 decoding structures investigated at Eve. 


clear all;
close all;

set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',16)
set(0,'defaultLineLineWidth',2.5)
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'DefaultLineMarkerSize',15);
set(0, 'defaultFigurePosition',  [-1267  44   1256    872])




sr_targetted = linspace(0.01,6,100);        % SR from 0.1 to 5 bit/channel use with 100 points
U = [2 4 8 16 32 64].';
sr_lin = 2.^sr_targetted; % linear value of the targetted secrecy rate

%% Alpha to inject

% Model 1
A1 = sr_lin - 1;
alpha1 = sqrt(A1.^2+A1)-A1;                    % alpha to inject in order to minimize Bob SNR expression when sigma_e = infinite and model = model1
SNRb1 = 10*log10(sqrt(A1.^2+A1)./((U+1).*(-2*A1.^2-2*A1+sqrt(A1.^2+A1).*(2*A1+1))));

% Model 2
B2 = U.*(sr_lin-1);
A2 = sr_lin.*(U+1).*(U+3)-B2;
alpha2 = (-2*B2 + sqrt(4.*A2.*B2+4.*B2.^2))./(2.*A2);                     % alpha to inject in order to minimize Bob SNR expression when sigma_e = infinite and model = model2
SNRb2 = 10*log10((alpha2.*A2 + B2)./((-alpha2.^2+alpha2).*(U+1).*U));

% Model 5
alpha5 = alpha1;                % alpha to inject in order to minimize Bob SNR expression when sigma_e = infinite and model = model5
SNRb5 = SNRb1;



%% Plot section
figure;


subplot(3,2,1)
plot(sr_targetted,1-alpha1.')
title('SDS Decoder')
xlabel('Targetted secrecy rate (bit/channel use)')
ylabel('Required AN energy to inject ($\%$)')
box on; grid on;
legendCell_tmp = cellstr(num2str('Does not depend on the BOR value'));
legendCell =[legendCell_tmp]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
legend(legendCell,'Location','best')

subplot(3,2,2)
plot(sr_targetted,SNRb1.')
title('SDS Decoder')
xlabel('Targetted secrecy rate (bit/channel use)')
ylabel('Required SNR at Bob (dB)')
box on; grid on;
legendCell_tmp = cellstr(num2str(U, 'BOR = %-d'));
legendCell =[legendCell_tmp]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
legend(legendCell,'Location','best')

subplot(1,2,1)
plot(sr_targetted,1-alpha2.')
title('MF Decoder')
xlabel('Targetted secrecy rate (bit/channel use)')
ylabel('Required AN energy to inject ($\%$)')
box on; grid on;
legendCell_tmp = cellstr(num2str(U, 'BOR = %-d'));
legendCell =[legendCell_tmp]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
legend(legendCell,'Location','best')

subplot(1,2,2)
plot(sr_targetted,SNRb2.')
title('MF Decoder')
xlabel('Targetted secrecy rate (bit/channel use)')
ylabel('Required SNR at Bob (dB)')
box on; grid on;
legendCell_tmp = cellstr(num2str(U, 'BOR = %-d'));
legendCell =[legendCell_tmp]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
legend(legendCell,'Location','best')

subplot(3,2,5)
plot(sr_targetted,1-alpha5.')
title('OC Decoder')
xlabel('Targetted secrecy rate (bit/channel use)')
ylabel('Required AN energy to inject ($\%$)');
box on; grid on;
legendCell_tmp = cellstr(num2str('Does not depend on the BOR value'));
legendCell =[legendCell_tmp]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
legend(legendCell,'Location','best')

subplot(3,2,6)
plot(sr_targetted,SNRb5.')
title('OC Decoder')
xlabel('Targetted secrecy rate (bit/channel use)')
ylabel('Required SNR at Bob (dB)')
box on; grid on;
legendCell_tmp = cellstr(num2str(U, 'BOR = %-d'));
legendCell =[legendCell_tmp]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
legend(legendCell,'Location','best')


