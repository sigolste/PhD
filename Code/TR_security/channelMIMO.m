function [H , H_freq , H_space , H_no_corr] = channelMIMO(fc , c , lambda_c , sigma_tau, Q , N , N_TX , N_RX , coef_space_TX , coef_space_RX , coef_freq)

% MIMO channel generation of dimension [N_RX x N_TX x Q]. 4 possible MIMO
% channels are generated depending on weither or not the spatial and/or
% temporal correlation is taken into account. The spatial correlation is
% implemented thanks to the Kronecker model applied on a 3 dimensional
% Rayleigh Channel. The frequency correlation is fixed by the delay spread
% of the channel. 


b_space_TX = coef_space_TX.*lambda_c;                                     % Distance between subsequent antennas (m)
b_space_RX = coef_space_RX.*lambda_c; 

r_TX = linspace( 0 , (N_TX-1)*b_space_TX , N_TX ) ;                       % Distance between TX antennas in meters
r_RX = linspace( 0 , (N_RX-1)*b_space_RX , N_RX ) ;                       % Distance between RX antennas in meters



dist_TX = r_TX./lambda_c;
dist_RX = r_RX./lambda_c;
        
        
% Frequency parameters                                                  % Delay spread (3us = urban ,  .5us = suburban, .2us = open areas)
delta_f_c = 1 / 2 / pi / sigma_tau ;                                    % Approximation of coherence bandwidth



% Matrices instantiation
H_no_corr = 1/sqrt(2)*( randn( N_RX, N_TX , Q ) ...
    + 1i * randn( N_RX, N_TX, Q ) );                                    % Uncorrelated MIMO channel of size: [N_RX x Q x N_TX] 


H_space = zeros(N_RX,N_TX,Q);                                           % Spatially correlated MIMO channel


H_freq = zeros(N_RX,N_TX,Q);                                            % Frequency correlated MIMO channel



H = zeros(N_RX,N_TX,Q);                                                 % Spatially + frequency correlated MIMO Channel



%delta_f_n_A = coef_freq_A.*delta_f_c;   
delta_f_n = coef_freq.*delta_f_c;   

b_subcar = delta_f_n./N;                                                % Subcarrier bandwidth


ind_subca = 0 : 1 : Q - 1 ;                                             % Subcarrier index
delta_f = ind_subca .* b_subcar ;                                       % Frequency separation between the first subcarrier and the others




%% Spatial correlation

% matrix generation
f_TX = fc + delta_f - delta_f( end ) / 2 ;                               % Frequency vector for AB MISO channel



[ F_TX , R_TX ] = ndgrid( f_TX , r_TX ) ;                                % Grid of frequencies of A MISO channel and positions of A antennas 

[ F_RX , R_RX ] = ndgrid( f_TX , r_RX ) ;                                % Grid of frequencies of B MISO channel and positions of B antennas 

% F_TX = F_RX obviously

% Kronecker model: independant spatial correlation at TX and at RX

Rho_space_TX = sinc( 2 * F_TX / c .* R_TX ) ;                                  % Spatial correlation for each frequency 
Rho_space_RX = sinc( 2 * F_RX / c .* R_RX ) ; 

rho_space_TX = zeros( length( r_TX ) , length( r_TX ) , Q ) ;           % Initialization of Correlation matrix 
rho_space_RX = zeros( length( r_RX ) , length( r_RX ) , Q ) ;           % Initialization of Correlation matrix 

for ii = 1 : Q
    rho_space_TX( : , : , ii ) = toeplitz( Rho_space_TX( ii , : ) ) ;     % Correlation matrix
    rho_space_RX( : , : , ii ) = toeplitz( Rho_space_RX( ii , : ) ) ;     % Correlation matrix
end
for jj = 1 : Q
    T_space_TX = squeeze( rho_space_TX( : , : , jj)  ) ;       % cholcov(...)         % Voir ce que ca fait exactement. Crée une loi multinormale pour chaque subcar?

    T_space_RX = squeeze( rho_space_RX( : , : , jj)  ) ;                              % Voir ce que ca fait exactement. Crée une loi multinormale pour chaque subcar?
end


% Kronecker channel
for ii = 1:Q
    H_space(:,:,ii) = sqrtm(T_space_RX)*H_no_corr(:,:,ii)*sqrtm(T_space_TX);
    %*H_space_B(:,:,ii) = tmp*sqrtm(T_space_A);
end



%% Frequency correlation

rho_freq = sigma_tau ./ ( 1 + 1i * 2 * pi * delta_f * sigma_tau ) ;         % Frequency variation of the correlation (first line of the covaraince matrix)
rho_freq = rho_freq ./ max( abs( rho_freq ) ) ;                             % Normalization
RHO_freq = toeplitz( rho_freq ) ;                                           % Covariance matrix
rho_freq = abs(RHO_freq);
T_freq = cholcov( RHO_freq ) ;



for ii = 1 : N_RX
    for jj = 1:N_TX
        H(ii,jj,:) = squeeze(H_space(ii,jj,:)).'*T_freq;
        H_freq(ii,jj,:) = squeeze(H_no_corr(ii,jj,:)).'*T_freq;
    end
end






