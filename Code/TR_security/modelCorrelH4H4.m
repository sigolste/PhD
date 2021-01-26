function H4H4 = modelCorrelH4H4(N,U,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of ∑_i ∑j!=i |Hb,n+iN|^4 |Hb,n+jN|^4 when frequency 
% correlation at Bob is introduced
% This term is needed for the determination of E[|SINR_b|^2]. In
% particular, this term is used for the computation of the variance of Bob
% SINR in order to better approximate the capacity at Bob when frequency
% correlation among subcarriers is introduced.
%
% INPUTS:
%   N : Number of symbol sent (N = Q/U)
%   U : Back-off rate
%   T : Choleski decomposition of correlation matrix RHO, s.t. T T^H = RHO
%
% OUTPUT:
%   H4H4 : Expected value of ∑_i ∑j!=i |Hb,n+iN|^4 |Hb,n+jN|^4
%
%
%
% Code started : 21.01.2021
% Last update  : 26.01.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ∑_i ∑_j!=i |Hb,n+iN|^4 |Hb,n+jN|^4 = 2 ∑_i ∑_j>i  |Hb,n+iN|^4 |Hb,n+jN|^4


% Derivation of |Hb_n+iN|^6 cf FC 26 recto

% 1 First term of hb_i^4 (with h^4) multiplied by hb_j^4 -> term4
% Derivation cf FC26 verso

tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            tmp1 = tmp1 + 6*abs(T(1+ii*N,kk))^4*abs(T(1+jj*N,kk))^4;
            for ll = 1:1+jj*N
                if ll ~= kk
                    tmp2 = tmp2 + abs(T(1+ii*N,kk))^4*abs(T(1+jj*N,ll))^4;
                    tmp3 = tmp3 + 6*abs(T(1+ii*N,kk))^4*abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^2;
                    for mm = 1:1+jj*N
                        if (ll ~= mm) && (kk ~= mm)
                            tmp4 = tmp4 + abs(T(1+ii*N,kk))^4*abs(T(1+jj*N,ll))^2*abs(T(1+jj*N,mm))^2;
                        end
                    end
                end
            end
        end
    end
end
TERM4_test = 4*(tmp1+tmp2+tmp3+tmp4);

% 2 Second term of hb_i^4 (with h^h h^*) multiplied by hb_j^4 -> term31
% Derivation cf FC26 verso




