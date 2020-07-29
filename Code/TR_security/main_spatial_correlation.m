clear all; close all; clc;

%**************************************************************************
%	Example code to call the function corr_spatial to generates 
%   nb_subca independent Rayleigh channel exhibiting spatial correlation 
%   assuming a three-dimensional Rayleigh channel
%
%   The spatial dependence of the correlation is given by eq. 24
%   from Ph. De Doncker, "Spatial Correlation Functions for Fields in Three
%   -Dimensional Rayleigh Channels," PIER, Vol. 40, 55-69, 2003
%
%   Inputs:
%   nb_subca: number of sub-carrier (or frequency samples within the
%   frequency bandwidth)
%   b_subca: bandwidth of each subcarrier (of frequency separation between 
%   frequency samples)
%   r: vector of the considered distances away from Bob
%   fc: carrier frequency
%   nb_realizations: number of samples of each sub-channel 
%
%   Outputs:
%   H: channel matrix whose size is (nb_realizations x nb_subca x nb_distance)
%
%   Code started the 13 September 2019
%   Last update: 13 September 2019
%
%   by Julien Sarrazin
%**************************************************************************

%% Inputs
nb_subca = 256 ;                                % Total number of subcarriers (i.e., total number of RV)
b_subca = .1e6 ;                                 % Bandwidth of each subcarrier
nb_realizations = 1e3 ;                         % Number of realizations
fc = 2e9 ;                                      % Carrier frequency
r = linspace( 0 , 1 , 40 ) ;                    % Eve's distance from Bob . Distance from 0 to 1 with 40 points

%% Random Variable Calculation
H = corr_spatial( nb_subca , b_subca , r , fc , nb_realizations ) ;
