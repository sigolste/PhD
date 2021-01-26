function H8 = modelCorrelH8(U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of ∑_i |Hb,n+iN|^8 when frequency correlation at Bob is 
% introduced.
% This term is needed for the determination of E[|SINR_b|^2]. In
% particular, this term is used for the computation of the variance of Bob
% SINR in order to better approximate the capacity at Bob when frequency
% correlation among subcarriers is introduced.
%
% INPUTS:
%   U : Back-off rate
%
% OUTPUT:
%   H8 : Expected value of ∑_i |Hb,n+iN|^8
%
%
%
% Code started : 21.01.2021
% Last update  : 21.01.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


H8 = 24*U;