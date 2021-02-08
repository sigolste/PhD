function H4 = modelCorrelH4(U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of ∑_i |Hb,n+iN|^4 when frequency correlation at Bob is 
% introduced.
% This term is needed for the determination of E[|SINR_b|] when frequency
% correlation among subcarriers is introduced.
%
% INPUTS:
%   U : Back-off rate
%
% OUTPUT:
%   H8 : Expected value of ∑_i |Hb,n+iN|^4
%
%
%
% Code started : 03.02.2021
% Last update  : 03.02.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


H4 = 2.*U;