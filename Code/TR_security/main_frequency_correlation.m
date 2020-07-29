clear all; close all; clc;

%**************************************************************************
%	Example code to call the function corr_frequency to generates 
%   nb_subca Rayleigh exhibiting frequency correlation assuming an 
%   exponenetial decaying power delay profile with delay spread sigma_tau
%
%   The frequency dependence of the correlation is given by eq. 3.18
%   chapter 3 of P. De Doncker course
%
%   Inputs:
%   nb_subca: number of sub-carrier (or frequency samples within the
%   frequency bandwidth)
%   b_subca: bandwidth of each subcarrier (of frequency separation between 
%   frequency samples)
%   sigma_tau: delay spread
%   nb_realizations: number of samples of each sub-channel 
%
%   Outputs:
%   H: channel matrix whose size is (nb_realizations x nb_subca)
%
%   Code started the 13 September 2019
%   Last update: 13 September 2019
%
%   by Julien Sarrazin
%**************************************************************************

%% Inputs
nb_subca = 256 ;                                % Total number of subcarriers (i.e., total number of RV)
b_subca = .5e6 ;                                % Bandwidth of each subcarrier
sigma_tau = 3e-6 ;                            % Delay spread
nb_realizations = 1e3 ;                         % Number of realizations
delta_f_c = 1 / 2 / pi / sigma_tau ;            % Approximation of coherence bandwidth

%% Random Variable Calculation
H = corr_frequency( nb_subca , b_subca , sigma_tau , nb_realizations ) ;

