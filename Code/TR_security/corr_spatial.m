function [ H , RHO_SPATIAL ] = corr_spatial( nb_subca , b_subca , r , fc , nb_realizations )

%**************************************************************************
%	Function that generates a complex normal distributed Rayleigh channel
%   with spatial correlation assuming a three-dimensional Rayleigh channel
%
%   The spatial dependence of the correlation is given by eq. 24
%   from Ph. De Doncker, "Spatial Correlation Functions for Fields in Three
%   -Dimensional Rayleigh Channels," PIER, Vol. 40, 55-69, 2003
%
%   Inputs:
%   nb_subca: number of sub-carrier (or frequency samples within the
%   frequency bandwidth)
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
    c = 3e8 ;                                                               % light velocity
    
    ind_subca = 0 : 1 : nb_subca - 1 ;                                      % Sub-carrier index
    delta_f = ind_subca * b_subca ;                                         % Frequency separation between the first subcarrier and the others    
    f = fc + delta_f - delta_f( end ) / 2 ;                                 % Frequency vector
    [ F , R ] = ndgrid( f , r ) ;
    
    Rho_spatial = sinc( 2 * F / c .* R ) ;                                  % Spatial correlation for each frequency --> pourquoi pas 2pi??
    
    RHO_SPATIAL = zeros( length( r ) , length( r ) , nb_subca ) ;           % Initialization of Correlation matrix   
    H = zeros( nb_realizations , nb_subca , length( r ) ) ;                 % Initialization of Channel with uncorrelated-frequency and correlated-spatial variations
    
    for ii = 1 : nb_subca
        RHO_SPATIAL( : , : , ii ) = toeplitz( Rho_spatial( ii , : ) ) ;     % Correlation matrix
    end
    for jj = 1 : nb_subca
        T = cholcov( squeeze( RHO_SPATIAL( : , : , jj) ) ) ;                % Voir ce que ca fait exactement. Crée une loi multinormale pour chaque subcar?
        H( : , jj , : ) = 1/sqrt(2)*( randn( nb_realizations , size( T , 1 ) ) ...
            + 1i * randn( nb_realizations , size( T , 1 ) ) ) * T ;         % Channel with uncorrelated-frequency and correlated-spatial variations
    end
 
    % Plot the spatial variation of the correlation 
    fig = figure ;
    hold on
    for ii = 1 : nb_subca
        plot( r , RHO_SPATIAL( 1 , : , ii ) )
    end
    xlabel('Distance')
    ylabel('spatial correlation $\rho$')
%     addToolbarExplorationButtons(fig)
    
    %% Verification
%     RHO_matrix_calculated = corrcoef( squeeze( H( : , 1 , : ) ) ) ;                               % Should be identical to RHO matrix if enough realizations
%     plot( r , real( RHO_matrix_calculated( 1 , : ) ) )
     
end

