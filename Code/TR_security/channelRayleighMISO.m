function y = channelRayleighMISO(N_TX,len, mu , sigma)
% Returns a rayleigh distributed channel to simualte the multipath effects.
% MISO configuration N_TX x 1 
%
% Inputs:
%   len : size of the channel
%   mu : mean of the distribution
%   sigma : standard deviation of the distribution
%   N_TX : number of TX antennas
%
% Output:
%   y : rayleigh distributed channel 

y = zeros(len*N_TX,len*N_TX);
vec  = 1/sqrt(2)*(sigma*randn(len*N_TX,1) + sigma*1j*randn(len*N_TX,1)) + mu;

y = diag(vec);
% for nn = 1:N_TX
%     y = [y diag(vec(nn*N_TX:,1));
% end
% y = diag(1/sqrt(2)*(sigma*randn(len,1) + sigma*1j*randn(len,1)) + mu);