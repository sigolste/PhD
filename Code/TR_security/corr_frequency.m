function [ H , abs_rho, T] = corr_frequency( nb_subca , b_subca , sigma_tau , nb_realizations )

%**************************************************************************
%	Function that generates correlated complex normal distributed 
%   Random Variables Rayleigh channels with frequency correlation fixed by
%   the delay spread of the channel
%
%   The frequency dependence of the correlation is given by eq. 4.17
%   chapter 4 of P. De Doncker course
%
%   Inputs:
%   nb_subca: number of sub-carrier (or frequency samples within the
%   frequency bandwidth)
%   b_subca: bandwidth of each subcarrier (or frequency separation between 
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
%

    ind_subca = 0 : 1 : nb_subca - 1 ;                                      % Sub-carrier index
    delta_f = ind_subca .* b_subca ;                                         % Frequency separation between the first subcarrier and the others
    rho = sigma_tau ./ ( 1 + 1i * 2 * pi * delta_f * sigma_tau ) ;          % Frequency variation of the correlation (first line of the covaraince matrix)
    rho = rho ./ max( abs( rho ) ) ;                                        % Normalization     Revient a avoir au final : rho = 1/(1+2pi j sigma_tau delta_f) car max(abs(rho)) = sigma_tau 
    RHO = toeplitz( rho ) ;                                                 % Covariance matrix
   
    T = cholcov( RHO ) ;                                                    % Choleski decomposition
    T_sqrt = sqrtm(RHO);
    Hw = 1/sqrt(2)*( randn( nb_realizations , size( T , 1 ) ) ...
        + 1i * randn( nb_realizations , size( T , 1 ) ) )  ;             % nb_subca correlated RV
    H = (ctranspose(T)*(Hw).').';
    abs_rho = abs(rho);
    T = ctranspose(T);
    %H = ones(nb_realizations, size(T,1));
    % Plot the frequency variation of the correlation 
%     figure
%     plot( delta_f , abs( rho ) )
%     xlabel('Frequency spacing')
%     ylabel('correlation \rho')
%     
%     %% Verification
%     RHO_matrix_calculated = corrcoef( H ) ;                               % Should be identical to RHO matrix if enough realizations
     
end

