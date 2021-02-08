function H2H2H2H2 = modelCorrelH2H2H2H2(N,U,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of ∑_i ∑j!=i ∑k!=i,j ∑l!=i,j,k |Hb,n+iN|^2 |Hb,n+jN|^2 
%                                             |Hb,n+kN|^2 |Hb,n+lN|^2  
% when frequency correlation at Bob is introduced                         
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
%   H4H2H2 : Expected value of 
% ∑_i ∑j!=i ∑k!=i,j ∑l!=i,j,k |Hb,n+iN|^2 |Hb,n+jN|^2 |Hb,n+kN|^2 
%                                                     |Hb,n+lN|^2                                                             
%
%
%
% Code started : 28.01.2021
% Last update  : 02.02.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Preliminaries:
% A) ∑_i ∑j!=i ∑k!=i,j ∑l!=i,j,k |Hb,n+iN|^2 |Hb,n+jN|^2 |Hb,n+kN|^2 
%                                                        |Hb,n+lN|^2 
%   = 24 ∑_i ∑j>i ∑k>j ∑l>k |Hb,n+iN|^2 |Hb,n+jN|^2 |Hb,n+kN|^2 |Hb,n+lN|^2 
%
% 24 since first term has 4 possibilities, the second has 3, the third 2
% and the last only 1. --> 4!  24 possibilities
%
% B) Derivation of |Hb_n+iN|^2|Hb_n+jN|^2 cf FC 28 recto



% 1 First term of hb_i^2 hb_j^2 (with h^4) multiplied by hb_k^2 hb_l^2  
% -> term4
% Derivation cf FC31 recto
tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = jj+1:U-1
            for ll = kk+1:U-1   %%%%%%%%
                
                for mm = 1:1+ii*N
                    tmp1 = tmp1 + 12*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                     abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,mm))^2;
                    for nn = 1:1+kk*N
                        if nn ~= mm
                            tmp2 = tmp2 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                            abs(T(1+kk*N,nn))^2*abs(T(1+ll*N,nn))^2;
                            tmp3 = tmp3 + 3*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                            abs(T(1+ll*N,mm))^2*abs(T(1+kk*N,nn))^2;
                            for oo = 1:1+ll*N
                                if (oo ~= nn) && (oo ~= mm)
                                    tmp4 = tmp4 + abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                                  abs(T(1+kk*N,nn))^2*abs(T(1+ll*N,oo))^2;
                                end
                            end
                            for oo = nn+1:1+kk*N
                                if oo ~= mm
                                    tmp7 = tmp7 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                                    real(T(1+kk*N,nn)*conj(T(1+ll*N,nn))*conj(T(1+kk*N,oo))*T(1+ll*N,oo));
                                end
                            end
                        end
                    end
                    for nn = 1:1+ll*N
                        if nn ~= mm
                            tmp5 = tmp5 + 3*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                            abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,nn))^2;
                        end
                    end
                    for nn = mm+1:1+kk*N
                        tmp6 = tmp6 + 12*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2*...
                                         real(T(1+kk*N,mm)*conj(T(1+ll*N,mm))*conj(T(1+kk*N,nn))*T(1+ll*N,nn));
                    end
                  
                end             %%%%%%%%
             end
        end
        
    end
end
TERM4_test = 48*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7);

% 2 Second term of hb_i^2 hb_j^2 (with h^2 h^^) multiplied by hb_k^2 hb_l^2  
% -> term22
% Derivation cf FC31 verso
tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
tmp8 = 0;
tmp9 = 0;
tmp10 = 0;
tmp11 = 0;
tmp12 = 0;
tmp13 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = jj+1:U-1
            for ll = kk+1:U-1   %%%%%%%%
                
                for mm = 1:1+ii*N
                    for nn = 1:1+jj*N
                        if nn ~= mm
                            tmp1 = tmp1 + 6*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                            abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,mm))^2;
                            tmp2 = tmp2 + 6*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                            abs(T(1+kk*N,nn))^2*abs(T(1+ll*N,nn))^2;
                            tmp3 = tmp3 + 4*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                            abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,nn))^2;
                           
                            for oo = 1:1+kk*N
                               if (oo ~= nn) && (oo ~= mm)
                                   tmp4 = tmp4 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                   abs(T(1+kk*N,oo))^2*abs(T(1+ll*N,oo))^2;
                                   tmp5 = tmp5 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                   abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,oo))^2;
                                   tmp6 = tmp4;
                                   for pp = 1:1+ll*N
                                       if (pp ~= oo) && (pp ~= nn) && (pp ~= mm)
                                           tmp9 = tmp9 + abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                         abs(T(1+kk*N,oo))^2*abs(T(1+ll*N,pp))^2;
                                       end
                                   end
                                   for pp = oo+1:1+kk*N
                                       tmp12 = tmp12 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                        real(T(1+kk*N,oo)*conj(T(1+ll*N,oo))*conj(T(1+kk*N,pp))*T(1+ll*N,pp));
                                   end
                               end
                            end
                            for oo = 1:1+ll*N
                                if (oo ~= nn) && (oo ~= mm)
                                    tmp7 = tmp7 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                    abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,oo))^2;
                                    tmp8 = tmp8 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                    abs(T(1+kk*N,nn))^2*abs(T(1+ll*N,oo))^2;
                                end
                            end
                        end
                        if nn ~= mm
                            for oo = mm+1:1+kk*N
                                if oo ~= nn
                                    tmp10 = tmp10 + 8*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                real(conj(T(1+kk*N,mm))*T(1+ll*N,mm)*T(1+kk*N,oo)*conj(T(1+ll*N,oo)));
                                end
                            end
                            for oo = nn+1:1+kk*N
                                if oo ~= mm
                                    tmp11 = tmp11 + 8*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                                real(conj(T(1+kk*N,nn))*T(1+ll*N,nn)*T(1+kk*N,oo)*conj(T(1+ll*N,oo)));
                                end
                            end
                        end
                    end
                    for nn = mm+1:1+jj*N
                        tmp13 = tmp13 + 8*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2*...
                                        real(T(1+kk*N,mm)*conj(T(1+ll*N,mm))*conj(T(1+kk*N,nn))*T(1+ll*N,nn));
                    end
              
                end             %%%%%%%%
             end
        end
        
    end
end
TERM22_test = 24*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 + tmp13);


% 3 Third term of hb_i^2 hb_j^2 (with conj(t t^* t t^*) h^2 h^2) multiplied
% by hb_k^2 hb_l^2  -> term22_bis
% Computation cfr FC32 recto (idem as second term with factor 48 in front of
% the formula instead of 24, nn = mm+1:n+iN instead of nn = 1:n+iN and
% replacing |t_{n+i*N,n}|^2*|t_{n+j*N,m}|^2 by 
% Re[t_{n+i*N,m} t_{n+j*N,m}^* t_{n+i*N,n}^* t_{n+j*N,n}]

tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
tmp8 = 0;
tmp9 = 0;
tmp10 = 0;
tmp11 = 0;
tmp12 = 0;
tmp13 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = jj+1:U-1
            for ll = kk+1:U-1   %%%%%%%%
                
                for mm = 1:1+ii*N
                    for nn = mm+1:1+ii*N
                            tmp1 = tmp1 + 6*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                            abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,mm))^2;
                            tmp2 = tmp2 + 6*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                            abs(T(1+kk*N,nn))^2*abs(T(1+ll*N,nn))^2;
                            tmp3 = tmp3 + 4*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                            abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,nn))^2;
                           
                            for oo = 1:1+kk*N
                               if (oo ~= nn) && (oo ~= mm)
                                   tmp4 = tmp4 + 2*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                   abs(T(1+kk*N,oo))^2*abs(T(1+ll*N,oo))^2;
                                   tmp5 = tmp5 + 2*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                   abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,oo))^2;
                                   tmp6 = tmp4;
                                   for pp = 1:1+ll*N
                                       if (pp ~= oo) && (pp ~= nn) && (pp ~= mm)
                                           tmp9 = tmp9 + real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                         abs(T(1+kk*N,oo))^2*abs(T(1+ll*N,pp))^2;
                                       end
                                   end
                                   for pp = oo+1:1+kk*N
                                       tmp12 = tmp12 + 2*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                        real(T(1+kk*N,oo)*conj(T(1+ll*N,oo))*conj(T(1+kk*N,pp))*T(1+ll*N,pp));
                                   end
                               end
                            end
                            for oo = 1:1+ll*N
                                if (oo ~= nn) && (oo ~= mm)
                                    tmp7 = tmp7 + 2*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                    abs(T(1+kk*N,mm))^2*abs(T(1+ll*N,oo))^2;
                                    tmp8 = tmp8 + 2*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                    abs(T(1+kk*N,nn))^2*abs(T(1+ll*N,oo))^2;
                                end
                            end
                            for oo = mm+1:1+kk*N
                                if oo ~= nn
                                    tmp10 = tmp10 + 8*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                real(conj(T(1+kk*N,mm))*T(1+ll*N,mm)*T(1+kk*N,oo)*conj(T(1+ll*N,oo)));
                                end
                            end
                            for oo = nn+1:1+kk*N
                                if oo ~= mm
                                    tmp11 = tmp11 + 8*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                                real(conj(T(1+kk*N,nn))*T(1+ll*N,nn)*T(1+kk*N,oo)*conj(T(1+ll*N,oo)));
                                end
                            end
                    end
                    for nn = mm+1:1+ii*N
                        tmp13 = tmp13 + 8*real(T(1+ii*N,mm)*conj(T(1+jj*N,mm))*conj(T(1+ii*N,nn))*T(1+jj*N,nn))*...
                                        real(T(1+kk*N,mm)*conj(T(1+ll*N,mm))*conj(T(1+kk*N,nn))*T(1+ll*N,nn));
                    end
              
                end             %%%%%%%%
             end
        end
        
    end
end

TERM22_bis_test = 48*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + ...
                      tmp8 + tmp9 + tmp10 + tmp11 + tmp12 + tmp13);


% 4 Fourth term of hb_i^2 hb_j^2 (with h h^* h h^*) multiplied 
% by hb_k^2 hb_l^2  -> term1111 
% Derivation FC32 recto.
tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = jj+1:U-1
            for ll = kk+1:U-1   %%%%%%%%
                for mm = 1:1+ii*N
                    for nn = mm+1:1+ii*N
                        for oo = 1:1+jj*N
                            if (oo ~= nn) && (oo ~= mm)
                                for pp = oo+1:1+jj*N
                                    if (pp ~= nn) && (pp ~= mm)
                                        tmp = tmp + real(conj(T(1+ii*N,mm))*T(1+kk*N,mm)*T(1+ii*N,nn)*conj(T(1+kk*N,nn))* ... 
                                                         conj(T(1+jj*N,oo))*T(1+ll*N,oo)*T(1+jj*N,pp)*conj(T(1+ll*N,pp)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
TERM1111_test = 16*24*tmp;

% 5 Fifth term of hb_i^2 hb_j^2 (with h^2 h^* h) multiplied  
% by hb_k^2 hb_l^2  -> term211
% Derivation FC32 recto&verso
tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = jj+1:U-1
            for ll = kk+1:U-1   %%%%%%%%
                
                for mm = 1:1+ii*N
                    for nn = 1:1+jj*N
                        for oo = nn+1:1+jj*N
                            if (oo ~= mm) && (nn ~= mm)
                                tmp1 = tmp1 + 4*abs(T(1+ii*N,mm))^2*abs(T(1+kk*N,mm))^2* ... 
                                       real(conj(T(1+jj*N,nn)) * T(1+ll*N,nn) * T(1+jj*N,oo) * conj(T(1+ll*N,oo)));
                            end
                        end
                    end
                    for nn = 1:1+kk*N
                        for oo = 1:1+jj*N
                            for pp = oo+1:1+jj*N
                               if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
                                   tmp2 = tmp2 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+kk*N,nn))^2* ... 
                                       real(conj(T(1+jj*N,oo)) * T(1+ll*N,oo) * T(1+jj*N,pp) * conj(T(1+ll*N,pp)));
                               end
                            end
                        end
                    end
                    for nn = 1:1+jj*N
                        for oo = nn+1:1+jj*N
                            if (oo ~= mm) && (nn ~= mm)
                                tmp3 = tmp3 + 4*abs(T(1+ii*N,mm))^2*abs(T(1+ll*N,mm))^2* ... 
                                       real(conj(T(1+jj*N,nn)) * T(1+kk*N,nn) * T(1+jj*N,oo) * conj(T(1+kk*N,oo)));
                            end
                        end
                    end
                    for nn = 1:1+ll*N
                        for oo = 1:1+jj*N
                            for pp = oo+1:1+jj*N
                               if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
                                   tmp4 = tmp4 + 2*abs(T(1+ii*N,mm))^2*abs(T(1+ll*N,nn))^2* ... 
                                       real(conj(T(1+jj*N,oo)) * T(1+kk*N,oo) * T(1+jj*N,pp) * conj(T(1+kk*N,pp)));
                               end
                            end
                        end
                    end
                end
                
            end
        end
    end
end
TERM211_test = 48*(tmp1 + tmp2 + tmp3 + tmp4);

% 6 Sixth term of hb_i^2 hb_j^2 (with h^2 h^* h) multiplied  
% by hb_k^2 hb_l^2  -> term211_bis 

tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = jj+1:U-1
            for ll = kk+1:U-1   %%%%%%%%T(1+ii*N,nn)*conj(T(1+ii*N,oo))
                
                for mm = 1:1+jj*N
                    for nn = 1:1+ii*N
                        for oo = nn+1:1+ii*N
                            if (oo ~= mm) && (nn ~= mm)
                                tmp1 = tmp1 + 4*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,mm))^2* ... 
                                       real(conj(T(1+ii*N,nn)) * T(1+ll*N,nn) * T(1+ii*N,oo) * conj(T(1+ll*N,oo)));
                            end
                        end
                    end
                    for nn = 1:1+kk*N
                        for oo = 1:1+ii*N
                            for pp = oo+1:1+ii*N
                               if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
                                   tmp2 = tmp2 + 2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,nn))^2* ... 
                                       real(conj(T(1+ii*N,oo)) * T(1+ll*N,oo) * T(1+ii*N,pp) * conj(T(1+ll*N,pp)));
                               end
                            end
                        end
                    end
                    for nn = 1:1+ii*N
                        for oo = nn+1:1+ii*N
                            if (oo ~= mm) && (nn ~= mm)
                                tmp3 = tmp3 + 4*abs(T(1+jj*N,mm))^2*abs(T(1+ll*N,mm))^2* ... 
                                       real(conj(T(1+ii*N,nn)) * T(1+kk*N,nn) * T(1+ii*N,oo) * conj(T(1+kk*N,oo)));
                            end
                        end
                    end
                    for nn = 1:1+ll*N
                        for oo = 1:1+ii*N
                            for pp = oo+1:1+ii*N
                               if  (nn ~= mm) && (oo ~= mm) && (oo ~= nn) && (pp ~= mm) && (pp ~= nn)
                                   tmp4 = tmp4 + 2*abs(T(1+jj*N,mm))^2*abs(T(1+ll*N,nn))^2* ... 
                                       real(conj(T(1+ii*N,oo)) * T(1+kk*N,oo) * T(1+ii*N,pp) * conj(T(1+kk*N,pp)));
                               end
                            end
                        end
                    end
                end
                
            end
        end
    end
end
TERM211_bis_test = 48*(tmp1 + tmp2 + tmp3 + tmp4);


% --> ∑_i ∑j!=i ∑k!=i,j ∑l!=i,j,k |Hb,n+iN|^2 |Hb,n+jN|^2 |Hb,n+kN|^2 
%                                                         |Hb,n+lN|^2
% cf FC33 recto
H2H2H2H2 = (TERM4_test + TERM22_test + TERM22_bis_test + TERM1111_test + TERM211_test + TERM211_bis_test );
