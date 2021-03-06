function H4H2H2 = modelCorrelH4H2H2(N,U,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of ∑_i ∑j!=i ∑k!=i,j |Hb,n+iN|^4 |Hb,n+jN|^2 |Hb,n+kN|^2   
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
%   H4H2H2 : Expected value of ∑_i ∑_j!=i ∑_k!=i,j |Hb,n+iN|^4 |Hb,n+jN|^2 
%                                                              |Hb,n+kN|^2 
%
%
%
% Code started : 26.01.2021
% Last update  : 27.01.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Preliminary:
% ∑_i ∑_j!=i ∑_k!=i,j |Hb,n+iN|^4 |Hb,n+jN|^2 |Hb,n+kN|^2 
%           = ∑_i ∑_j!=i ∑_k!=i,j |Hb,n+iN|^2 |Hb,n+jN|^2 |Hb,n+kN|^4 
%           = 2∑_i ∑_j>i ∑_k!=i,j |Hb,n+iN|^2 |Hb,n+jN|^2 |Hb,n+kN|^4 


% Derivation of |Hb_n+iN|^2|Hb_n+jN|^2 cf FC 28 recto

% 1 First term of hb_i^2 hb_j^2 (with h^4) multiplied by hb_k^4  -> term4
% Derivation cf FC28 recto/verso

tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (kk ~= ii) && (kk ~= jj)
                for ll = 1:1+ii*N
                    tmp1 = tmp1 + 6*abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,ll))^2*abs(T(1+kk*N,ll))^4;
                    for mm = 1:1+kk*N
                        if mm ~= ll
                            tmp2 = tmp2 + abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,ll))^2*abs(T(1+kk*N,mm))^4;
                            tmp3 = tmp3 + 6*abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,ll))^2*abs(T(1+kk*N,ll))^2*abs(T(1+kk*N,mm))^2;
                            for nn = 1:1+kk*N
                                if (nn ~= mm) && (nn ~= ll)
                                    tmp4 = tmp4 + abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,ll))^2*abs(T(1+kk*N,mm))^2*abs(T(1+kk*N,nn))^2;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
            
TERM4_test = 8*(tmp1 + tmp2 + tmp3 + tmp4);

% 2 Second term of hb_i^2 hb_j^2 (with t^2 t^2 h^2 h^2) multiplied by 
% hb_k^4  -> term22
% Derivation cf FC29 recto
%
tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (kk ~= ii) && (kk ~= jj)
                
                for ll = 1:1+ii*N
                    for mm = 1:1+jj*N
                        if mm ~= ll
                            tmp1 = tmp1 + 3*abs(T(1+ii*N,ll))^2*abs(T(1+kk*N,ll))^4*abs(T(1+jj*N,mm))^2;
                            tmp2 = tmp2 + 3*abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,mm))^4;
                            tmp3 = tmp3 + 4*abs(T(1+ii*N,ll))^2*abs(T(1+kk*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,mm))^2;
                            for nn = 1:1+kk*N
                                if (nn ~= mm) && (nn ~= ll)
                                    tmp4 = tmp4 + abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,nn))^4;
                                    tmp5 = tmp5 + 4*abs(T(1+ii*N,ll))^2*abs(T(1+kk*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,nn))^2;
                                    tmp6 = tmp6 + 4*abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,mm))^2*abs(T(1+kk*N,nn))^2;
                                    for oo = 1:1+kk*N
                                        if (oo ~= nn) && (oo ~= mm) && (oo ~= ll)
                                            tmp7 = tmp7 + abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+kk*N,nn))^2*abs(T(1+kk*N,oo))^2;
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
end          
TERM22_test = 4*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7);

% 3 Third term of hb_i^2 hb_j^2 (with conj(t t^* t t^*) h^2 h^2) multiplied
% by hb_k^4  -> term22_bis 
% Computation cf FC29 recto (idem as second term with factor 8 in front of
% the formula instead of 4 and mm = ll+1:n+jN instead of mm = 1:n+jN

tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (kk ~= ii) && (kk ~= jj)
                
                for ll = 1:1+ii*N
                    for mm = ll+1:1+jj*N
                        if mm ~= ll
                            tmp1 = tmp1 + 3*real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,ll))^4;
                            tmp2 = tmp2 + 3*real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,mm))^4;
                            tmp3 = tmp3 + 4*real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,ll))^2*abs(T(1+kk*N,mm))^2;
                            for nn = 1:1+kk*N
                                if (nn ~= mm) && (nn ~= ll)
                                    tmp4 = tmp4 + real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,nn))^4;
                                    tmp5 = tmp5 + 4*real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,ll))^2*abs(T(1+kk*N,nn))^2;
                                    tmp6 = tmp6 + 4*real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,mm))^2*abs(T(1+kk*N,nn))^2;
                                    for oo = 1:1+kk*N
                                        if (oo ~= nn) && (oo ~= mm) && (oo ~= ll)
                                            tmp7 = tmp7 + real(conj(T(1+ii*N,ll))*T(1+jj*N,ll)*T(1+ii*N,mm)*conj(T(1+jj*N,mm)))*abs(T(1+kk*N,nn))^2*abs(T(1+kk*N,oo))^2;
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
end          
TERM22_bis_test = 8*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7);



% 4 Fourth term of hb_i^2 hb_j^2 (with h h^* h h^*) multiplied 
% by hb_k^4  -> term1111 
% Derivation FC29 vero.
tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (kk ~= ii) && (kk ~= jj)
                for ll = 1:1+ii*N
                    for mm = ll+1:1+ii*N
                        for nn = 1:1+jj*N
                            for oo = nn+1:1+jj*N
                                if (oo ~= mm) && (oo ~= ll) && (nn ~= mm) && (nn ~= ll)
                                    tmp = tmp + real(conj(T(1+ii*N,ll))*T(1+ii*N,mm)*T(1+jj*N,nn)*conj(T(1+jj*N,oo)) ...
                                         * T(1+kk*N,ll)*conj(T(1+kk*N,mm))*conj(T(1+kk*N,nn))*T(1+jj*N,oo));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

TERM1111_test = 32*tmp;

% 5 Fifth term of hb_i^2 hb_j^2 (with h^2 h^* h) multiplied  
% by hb_k^4  -> term211 
% Derivation FC29 verso / FC30 recto
tmp1 = 0;
tmp2 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (kk ~= ii) && (kk ~= jj)
                for ll = 1:1+ii*N
                    for mm = 1:1+jj*N
                        for nn = mm+1:1+jj*N
                            if (nn ~= ll) && (mm ~= ll)    
                                tmp1 = tmp1 + 2*abs(T(1+ii*N,ll))^2*abs(T(1+kk*N,ll))^2* ...
                                        real(conj(T(1+jj*N,mm))*T(1+kk*N,mm)*T(1+jj*N,nn)*conj(T(1+kk*N,nn)));    
                            end
                        end
                    end
                end
                for ll = 1:1+ii*N
                    for mm = 1:1+kk*N
                        for nn = 1:1+jj*N
                            for oo = nn+1:1+jj*N
                                if (oo ~= mm) && (oo ~= ll) && (nn ~= ll) && (nn ~= ll) && (mm ~= ll)
                                    tmp2 = tmp2 + abs(T(1+ii*N,ll))^2*abs(T(1+kk*N,mm))^2* ...
                                        real(conj(T(1+jj*N,nn))*T(1+kk*N,nn)*T(1+jj*N,oo)*conj(T(1+kk*N,oo)));  
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end      
TERM211_test = 16*(tmp1 + tmp2);


% 6 Sixth term of hb_i^2 hb_j^2 (with h^2 h^* h) multiplied  
% by hb_k^4  -> term211_bis 
tmp1 = 0;
tmp2 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 0:U-1
            if (kk ~= ii) && (kk ~= jj)
                for ll = 1:1+jj*N
                    for mm = 1:1+ii*N
                        for nn = mm+1:1+ii*N
                            if (nn ~= ll) && (mm ~= ll)    
                                tmp1 = tmp1 + 2*abs(T(1+jj*N,ll))^2*abs(T(1+kk*N,ll))^2* ...
                                        real(conj(T(1+ii*N,mm))*T(1+kk*N,mm)*T(1+ii*N,nn)*conj(T(1+kk*N,nn)));    
                            end
                        end
                    end
                end
                for ll = 1:1+jj*N
                    for mm = 1:1+kk*N
                        for nn = 1:1+ii*N
                            for oo = nn+1:1+ii*N
                                if (oo ~= mm) && (oo ~= ll) && (nn ~= ll) && (nn ~= ll) && (mm ~= ll)
                                    tmp2 = tmp2 + abs(T(1+jj*N,ll))^2*abs(T(1+kk*N,mm))^2* ...
                                        real(conj(T(1+ii*N,nn))*T(1+kk*N,nn)*T(1+ii*N,oo)*conj(T(1+kk*N,oo)));  
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end      
TERM211_bis_test = 16*(tmp1 + tmp2);







% --> ∑_i ∑_j!=i ∑_k!=i,j |Hb,n+iN|^4 |Hb,n+jN|^2 |Hb,n+kN|^2 
% cf FC30 recto
H4H2H2 = (TERM4_test + TERM22_test + TERM22_bis_test + TERM1111_test + TERM211_test + TERM211_bis_test );
