function [symb_noisy, noise_power] = TR_security_addNoise(symb, EbN0, M, symb_power)
% ADDNOISE(symb_TX_IQ , EbN0 , M)
%
% Add AWGN noise to a signal for a given Eb/N0 and modulation order.
%
% Arguments:
%     symb: list of the transmitted complex symbols.
%     EbN0: energy per bit to noise power spectral density ratio.
%     M: modulation order.
%
% Return:
%     symb_noise: noisy symbols.
%

k = log2(M);

snr  = EbN0 + 10*log10(k);

noise_PSD = symb_power/10^(snr/10);

noise = sqrt(noise_PSD).*(randn(size(symb,1),size(symb,2))...
        + 1j*randn(size(symb,1),size(symb,2)));

symb_noisy = symb + noise;
noise_power = 1/length(noise)*sum(abs(noise).^2);


end