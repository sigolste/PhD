function [H,rho_space,rho_freq] = corr_spatial_frequency(nb_subcar , b_subcar , r , fc , sigma_tau , nb_realizations)
% Generate nb_realitions channels that are spatially correlated between
% them and that have a frequency correlation between their subcarriers


ind_subca = 0 : 1 : nb_subcar - 1 ;                                     % Sub-carrier index
delta_f = ind_subca .* b_subcar ;                                       % Frequency separation between the first subcarrier and the others
rho = sigma_tau ./ ( 1 + 1i * 2 * pi * delta_f * sigma_tau ) ;          % Frequency variation of the correlation (first line of the covaraince matrix)
rho = rho ./ max( abs( rho ) ) ;                                        % Normalization
RHO = toeplitz( rho ) ;                                                 % Covariance matrix
rho_freq = abs(RHO);
T = cholcov( RHO ) ;
    
[H_space ,rho_space] = corr_spatial( nb_subcar , b_subcar , r , fc , 1) ;  % Spatially correlated channels with no frequency correlations. 
H_space = (squeeze(H_space)).'; 
H = H_space*T;                                                              % Apply frequency correlation to the spatially correlated channels
