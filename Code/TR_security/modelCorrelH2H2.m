function H2H2 = modelCorrelH2H2(N,U,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of  ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2 when frequency 
% correlation at Bob is introduced                         
% This term is needed for the determination of E[|SINR_b|] when frequency
% correlation among subcarriers is introduced.
%
% INPUTS:
%   N : Number of symbol sent (N = Q/U)
%   U : Back-off rate
%   T : Choleski decomposition of correlation matrix RHO, s.t. T T^H = RHO
%
% OUTPUT:
%   H2H2 : Expected value of 
%   ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2
%                                                                                                                
%
% Code started : 03.02.2021
% Last update  : 03.02.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tmp1 = 0;    % |t_ik|^2 |t_jk|^2 |h_k|^4
tmp2 = 0;    % |t_ik|^2 |t_jl|^2 |h_k|^2 |h_l|^2
tmp3 = 0;         % t_ik t_il^* t_jk t_jl^* |h_k|^2 |h_l|^2

for ii = 0: U-1
    for jj = ii+1 : U-1
        for kk = 1 : 1+ii*N
            tmp1 = tmp1 + 2*abs(T(1+ii*N,kk))^2 * abs(T(1+jj*N,kk))^2;
        end
    end
end
for ii = 0: U-1
    for jj = ii+1 : U-1
        for kk = 1 : 1+ii*N
            for ll = 1 : 1+jj*N
                if ll ~= kk
                    tmp2 = tmp2 + abs(T(1+ii*N,kk))^2*abs(T(1+jj*N,ll))^2;
                end   
            end
        end
    end
end

for ii = 0: U-1
    for jj = ii+1 : U-1
        for kk = 1 : 1+ii*N
            for ll = kk+1 : 1+ii*N
                %[{1+ii*N,kk}  {1+ii*N,ll} , {1+jj*N,kk} {1+jj*N,ll}] % tuples
                tmp3 = tmp3 + 2*real(conj(T(1+ii*N,kk))*T(1+ii*N,ll)* ...
                                     T(1+jj*N,kk) * conj(T(1+jj*N,ll)));
            end
        end
    end
end


H2H2 = 2*(tmp1 + tmp2 + tmp3);       % ∑_i ∑_{j!=i}  |h_i|^2 |h_j|^2