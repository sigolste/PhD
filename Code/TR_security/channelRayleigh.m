function y = channelRayleigh(len, mu , sigma)
% Returns a rayleigh distributed channel to simualte the multipath effects. 
%
% Inputs:
%   len : size of the channel
%   mu : mean of the distribution
%   sigma : standard deviation of the distribution
%
% Output:
%   y : rayleigh distributed channel 

y = diag(1/sqrt(2)*(sigma*randn(len,1) + sigma*1j*randn(len,1)) + mu);