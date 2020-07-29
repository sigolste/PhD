function [sym_RX] = TR_security_OFDM(signal_TX, nb_subcarriers, nb_blocks, length_CP , h )%, EbN0 , M,rrcfilter,sps,span)

%% ------------------------------------------------------------------------
% Perform OFDM for the time reversal system. `
% 
% INPUTS:
%   sym_to_transmit : symbol stream to be transmitted
%   nb_subcarriers  : nb of subcarriers in the channel frequency response
%   nb_blocks       : nb of OFDM blocks to transmit
%   length_CP       : size of cyclic prefix
%   h               : propagation channel impulse response (Alice or Eve)
%   EbN0            : Energy per bit to noise PSD ratio
%   rrcfilter       : Root raised cosine filter
%   sps             : Samples per symbol
%   span            : filter span in symbols
%
% OUTPUT:
%   sym_RX          : Received symbol stream
%
%
% By Pierre Baudoux & Sidney Golstein
% Last Update: March 4th 2019
%%


%%%%%%%%%%%%%%%%%%%%%%%
%%% Inverse Fourier %%%
%%%%%%%%%%%%%%%%%%%%%%%


%signal_before_CP = size(signal_TX)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cyclic Prefix addition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = signal_TX;
prefix = signal_TX(end-length_CP+1:end,:);


signal_TX = [prefix; tmp];
%signal_w_CP = size(signal_TX)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preambule addition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%
%%% Channel & Noise %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Deparallelization to matched filter and do the convolution
% nb_subcar x n_block => fft_size*n_block x 1
signal_TX = reshape(signal_TX, (nb_subcarriers+length_CP)*nb_blocks, 1);


% signal_shaped = upfirdn(signal_TX, rrcfilter, sps);


%%% Channel convolution
% signal_MPC = conv(signal_shaped,h);
signal_MPC = conv(signal_TX,h);
signal_MPC(end-length(h)+2:end,:) = [];


%%% Noise
%signal_MPC_noise = TR_security_addNoise(signal_MPC , EbN0 , M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cyclic Prefix removal %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% RRC filter ; Parallelization to remove CP, do the fft and despread
%signal_MPC_matched = upfirdn(signal_MPC, rrcfilter, 1, sps);
%signal_MPC_matched = signal_MPC_matched( span + 1 : end - span ) ;
%signal_MPC_matched = reshape(signal_MPC_matched, nb_subcarriers+length_CP, nb_blocks);           % OF NOISE: signal_MPC_noise instead of signal_MPC
signal_MPC_reshaped = reshape(signal_MPC, nb_subcarriers+length_CP, nb_blocks);           % OF NOISE: signal_MPC_noise instead of signal_MPC
signal_MPC_reshaped(1:length_CP, :) = [];
%signal_parallelized = size(signal_MPC_noise)

%%%%%%%%%%%%%%%
%%% Fourier %%%
%%%%%%%%%%%%%%%

% column-wise FFT
sym_RX = fft(signal_MPC_reshaped);

