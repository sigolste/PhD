function [H,rho_space,rho_freq] = corr_spatial_frequency(nb_subcar , b_subcar , r , fc , sigma_tau , nb_realizations)
% Generate nb_realitions channels that are spatially correlated between
% them and that have a frequency correlation between their subcarriers

%% Frequency correlation: 
ind_subca = 0 : 1 : nb_subcar - 1 ;                                         % Sub-carrier index
delta_f = ind_subca .* b_subcar ;                                           % Frequency separation between the first subcarrier and the others
rho_freq = sigma_tau ./ ( 1 + 1i * 2 * pi * delta_f * sigma_tau ) ;         % Frequency variation of the correlation (first line of the covaraince matrix)
rho_freq = rho_freq ./ max( abs( rho_freq ) ) ;                             % Normalization
RHO_freq = toeplitz( rho_freq ) ;                                           % Covariance matrix
rho_freq = abs(RHO_freq);
T_freq = cholcov( RHO_freq ) ;
H_freq = 1/sqrt(2)*( randn( nb_realizations , size( T_freq , 1 ) ) ...
        + 1i * randn( nb_realizations , size( T_freq , 1 ) ) ) * T_freq ;             % nb_subca correlated RV
    
    
c = 3e8 ;                                                               % light velocity   
f = fc + delta_f - delta_f( end ) / 2 ;                                 % Frequency vector
[ F , R ] = ndgrid( f , r ) ;                                           % Grid of frequencies and positions
Rho_spatial = sinc( 2 * F / c .* R ) ;                                  % Spatial correlation for each frequency --> pourquoi pas 2pi??
rho_space = zeros( length( r ) , length( r ) , nb_subcar ) ;           % Initialization of Correlation matrix 

for ii = 1 : nb_subcar
    rho_space( : , : , ii ) = toeplitz( Rho_spatial( ii , : ) ) ;     % Correlation matrix
end
for jj = 1 : nb_subcar
    T_space = cholcov( squeeze( rho_space( : , : , jj) ) ) ;                % Voir ce que ca fait exactement. Crée une loi multinormale pour chaque subcar?
end

H = T_space*H_freq;                                                     

% [H_space ,rho_space] = corr_spatial( nb_subcar , b_subcar , r , fc , 1) ;  % Spatially correlated channels with no frequency correlations. 
% H_space = (squeeze(H_space)).'; 
% H = H_space*T;                                                              % Apply frequency correlation to the spatially correlated channels
