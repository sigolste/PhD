function [ RV ] = corr_rayleigh( nb_subca , bc , variance , nb_realizations )

%**************************************************************************
%	Function that generates correlated complex normal distributed 
%   Random Variables for correlated Rayleigh frequency channels
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

    tau_rms = sqrt( 3 ) / ( 2 * pi * ( bc - 1 ) ) ;                         % RMS delay spread
    delta_f = 0 : 1 : nb_subca - 1 ;                                        % Sub-carrier index
    rho = 1 ./ sqrt( 1 + ( 2 * pi * delta_f * tau_rms ).^2 ) ;              % Frequency variation of the correlation 

    Gamma = variance * rho ;                                                % First line of the covaraince matrix
    GAMMA = toeplitz( Gamma ) ;                                             % Covariance matrix
    T = cholcov( GAMMA ) ;                                                  % Choleski decomposition
    RV = ( randn( nb_realizations , size( T , 1 ) ) ...
        + 1i * randn( nb_realizations , size( T , 1 ) ) ) * T ;             % nb_subca correlated RV
    RV = RV./(abs(mean(RV))) ;
    RV = RV.';

    %% Plot the frequency variation of the correlation 
    figure
    plot( delta_f , rho )
    xlabel('Nb of subcarrier')
    ylabel('correlation $\rho$')

end

