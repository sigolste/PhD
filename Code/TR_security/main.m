clear all; close all; clc;

%**************************************************************************
%	Example code to call the function corr_rayleigh to generates correlated 
%   complex normal distributed Random Variables for correlated Rayleigh
%   frequency channels
%
%   The frequency dependence of the correlation is given by eq. 11.16p252
%   in the book: "Antennas and Propagation for Wirless Communication 
%   Systems" by Saunders (Wiley)
%
%   Code started the 29 April 2019
%   Last update: 29 April 2019
%
%   by Julien Sarrazin
%**************************************************************************

%% Inputs
nb_subca = 256 ;                                % Total number of subcarriers (i.e., total number of RV)
bc = 64 ;                                       % Coherence bandwidth (number of subcarriers)
variance = .5 ;                                 % Variance of the sub-channels (here considered constant throught out the total bandwidth)
nb_realizations = 1e6 ;                         % Number of realization of each RV

%% Random Variable Calculation
RV = corr_rayleigh( nb_subca , bc , variance , 1 ).' ;
RV = RV./(abs(mean(RV))) ;
%% Verification
GAMMA_calculated =  cov( RV ) ;


