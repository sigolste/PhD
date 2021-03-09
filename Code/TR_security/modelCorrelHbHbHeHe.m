function HbHbHeHe = modelCorrelHbHbHeHe(N,U,Tb,Te)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of  ∑_i ∑j!=i Hb,n+iN Hb,n+jN^* He,n+iN^* He,n+jN when 
% frequency  correlation at Bob and Eve is introduced                         
% This term is needed for the determination of E[|SINR_b|] when frequency
% correlation among subcarriers is introduced.
%
% INPUTS:
%   N : Number of symbol sent (N = Q/U)
%   U : Back-off rate
%   Tb : Choleski decomposition of correlation matrix RHO @B: Tb Tb^H = RHO
%   Te : Choleski decomposition of correlation matrix RHO @E: Te Te^H = RHO
%
% OUTPUT:
%   HbHbHeHe : Expected value of: 
%              ∑_i ∑j!=i Hb,n+iN Hb,n+jN^* He,n+iN^* He,n+jN
%                                                                                                                
%
% Code started : 02.03.2021
% Last update  : 02.03.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NB: derivation cfr FC38verso/39recto

tmp1 = 0;
tmp2 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        if ii~= jj
            for kk = 1:1+ii*N
                tmp1 = tmp1 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Te(1+ii*N,kk))*Te(1+jj*N,kk));
                for ll = 1:1+ii*N
                    if ll ~= kk
                        tmp2 = tmp2 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk)))*real(conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    end
                end
               
            end
        end
    end
end
HbHbHeHe = tmp1+tmp2;
