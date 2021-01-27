function var = modelVariance(alpha,N,U,T,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of variance of Bob/Eve SINR when frequency correlation is 
% introduced.
% This term is needed in order to better approximate the capacity at Bob/
% Eve when frequency correlation among subcarriers is introduced.
% 
% var = E[X^2] - (E[X])^2 where X is Bob/Eve SINR
%
% INPUTS:
%   alpha : percentage of data energy injected in the communication
%   N : Number of symbol sent (N = Q/U)
%   U : Back-off rate
%   T : Choleski decomposition of correlation matrix RHO, s.t. T T^H = RHO
%   type : Scenario that has to be modelized
%
% OUTPUT:
%   var : The modelized variance of Bob/Eve SINR
%
%
%
% Code started : 27.01.2021
% Last update  : 27.01.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preliminary:
% Approximation of E[1+ log(X)] = log(1 + E[X]) - var(X)/(2( 1+E[X] )^2)
% where X = SINR_b ou SINR_e

switch type
    case "bob_noCorrel"
    case "bob_correl"
    case "eve_decod1_noCorrel"
    case "eve_decod2_noCorrel"
    case "eve_decod5_noCorrel"
    case "eve_decod1_Correl"
    case "eve_decod2_Correl"
    case "eve_decod5_Correl"
end