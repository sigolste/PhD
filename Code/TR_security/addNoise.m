function [noise, noise_energy] = addNoise(symb, snr, symb_energy)
% Add AWGN noise to a signal for a given Eb/N0, a given number of encoded
% bits per symbol and given a reference power
%
% INPUTS:
%     symb: list of the transmitted complex symbols.
%     EbN0: energy per bit to noise power spectral density ratio.
%     snr: Signal to noise ratio
%     symb_energy : reference energy over which the noise power is computed 
%
% OUTPUTS:
%     symb_noise: symbol + noise.
%     noise: noise component of the symbol 
%    

noise_PSD = symb_energy/10^(snr/10); % power spectral density

noise = sqrt(noise_PSD/2).*(randn(size(symb,1),size(symb,2))...
        + 1j*randn(size(symb,1),size(symb,2)));

%noise = noise/sqrt(energy(noise))*sqrt(sigma);
% symb_noisy = symb + noise;
noise_energy = energy(noise);


end