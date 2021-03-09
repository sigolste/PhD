function Hb2Hb2He2He2 = modelCorrelHb2Hb2He2He2(N,U,Tb,Te)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of  ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2|He,n+iN|^2 |He,n+jN|^2
% when frequency correlation at Bob and at Eve is introduced                         
% This term is needed for the determination of E[|SINR_e|] when frequency
% correlation among subcarriers is introduced.
%
% INPUTS:
%   N : Number of symbol sent (N = Q/U)
%   U : Back-off rate
%   Tb : Choleski decomposition of correlation matrix RHO: Tb Tb^H = RHO
%
% OUTPUT:
%   H2H2 : Expected value of 
%   ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2
%                                                                                                                
%
% Code started : 03.02.2021
% Last update  : 08.02.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First term of ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2|He,n+iN|^2 |He,n+jN|^2
% Derivatuon cfr FC39verso/40recto
tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
tmp4 = 0;
tmp5 = 0;
tmp6 = 0;
tmp7 = 0;
tmp8 = 0;


 for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj
        for kk = 1:1+ii*N
            tmp1 = tmp1 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,kk))^2;
            for ll = 1:1+ii*N
                if ll ~= kk
                    tmp2 = tmp2 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2;
                    tmp3 = tmp3 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2;
                end
            end
            for ll = 1:1+jj*N
                if ll ~= kk
                    tmp4 = tmp4 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,kk))^2;
                    for mm = 1:1+ii*N
                        if (mm ~= ll) && (mm ~= kk)
                            tmp8 = tmp8 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,mm))^2;
                        end
                    end
                end
            end
            for ll = kk+1:1+ii*N
                tmp5 = tmp5 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,kk))^2;
                tmp6 = tmp6 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,ll))^2;
                for mm = 1:1+ii*N
                    if (mm ~= ll) && (mm ~= kk)
                        tmp7 = tmp7 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,mm))^2;
                    end
                end
            end
        end
        end
    end
end 


H4E_test = tmp1 + tmp2 + tmp3 + tmp4 +  tmp5 + tmp6 + tmp7 + tmp8;





% Secodnd term of ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2|He,n+iN|^2 |He,n+jN|^2
% Derivatuon cfr FC40recto/verso
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
tmp14 = 0;
tmp15 = 0;
tmp16 = 0;
tmp17 = 0;
tmp18 = 0;
tmp19 = 0;
 for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj
            for kk = 1:1+ii*N
                for ll = 1:1+jj*N
                    if ll ~= kk
                        tmp1 = tmp1 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                        tmp2 = tmp2 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                        for mm = 1:1+ii*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp3 = tmp3 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,kk))^2;
                                tmp4 = tmp4 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,ll))^2;
                                for nn = 1:1+jj*N
                                    if (nn ~= mm) && (nn ~= ll) && (nn ~= kk)
                                        tmp5 = tmp5 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,nn))^2;
                                    end
                                end
                            end
                        end
                        for mm = 1:1+jj*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp6 = tmp6 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,mm))^2;
                            end
                        end
                    end
                end
                for ll = 1:1+ii*N
                    if ll ~= kk
                        tmp7 = tmp7 + 2*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                        tmp8 = tmp8 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,kk))^2;
                        for mm = 1:1+jj*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp9 = tmp9 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2;
                                tmp10 = tmp10 + abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2;
                            end
                        end
                    end
                end
                for ll = kk+1:1+ii*N
                    tmp12 = tmp12 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,ll))^2;
                    tmp15 = tmp15 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,kk))^2;
                    for mm = 1:1+jj*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp13 = tmp13 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,kk))^2*abs(Te(1+jj*N,mm))^2;
                                tmp14 = tmp14 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,ll))^2*abs(Te(1+jj*N,mm))^2;
                                for nn = 1:1+ii*N
                                    if (nn ~= mm) && (nn ~= ll) && (nn ~= kk)
                                        tmp19 = tmp19 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,nn))^2*abs(Te(1+jj*N,mm))^2;
                                    end
                                end
                            end
                    end
                    for mm = 1:1+ii*N
                            if (mm ~= ll) && (mm ~= kk)
                                tmp16 = tmp16 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,kk))^2;
                                tmp17 = tmp17 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,ll))^2;
                                for nn = 1:1+jj*N
                                    if (nn ~= mm) && (nn ~= ll) && (nn ~= kk)
                                        tmp18 = tmp18 + 2*real(Tb(1+ii*N,kk)*conj(Tb(1+ii*N,ll))*conj(Tb(1+jj*N,kk))*Tb(1+jj*N,ll))*abs(Te(1+ii*N,mm))^2*abs(Te(1+jj*N,nn))^2;
                                    end
                                end
                            end
                    end
                end
            end
        end
    end
 end
H22E_test = tmp1 + tmp2 + tmp3 + tmp4 +  tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 +  tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;


% Third term of ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^2|He,n+iN|^2 |He,n+jN|^2
% Derivation cfr FC40verso/41recto

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
tmp14 = 0;
 for ii = 0:U-1
    for jj = 0:U-1
        if ii ~= jj 
            for kk = 1:1+ii*N
                for ll = kk+1:1+ii*N
                    tmp1 = tmp1 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp2 = tmp2 + 4*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp3 = tmp3 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp4 = tmp4 + 2*abs(Tb(1+ii*N,ll))^2*abs(Tb(1+jj*N,kk))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    tmp5 = tmp5 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll))*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                    for mm = 1:1+ii*N
                        if (mm ~= kk) && (mm ~= ll)
                            tmp6 = tmp6 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,mm)*conj(Te(1+jj*N,mm))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                            for nn = mm+1:1+ii*N
                                if (nn ~= kk) && (nn ~= ll)
                                    tmp7 = tmp7 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll))*real(Te(1+ii*N,mm)*conj(Te(1+jj*N,mm))*conj(Te(1+ii*N,nn))*Te(1+jj*N,nn));
                                end
                            end
                        end
                    end
                    for mm = ll+1:1+ii*N
                        if (mm ~= kk) 
                            tmp8 = tmp8 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,ll))*Tb(1+jj*N,ll))*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                            tmp9 = tmp9 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,mm))*Tb(1+jj*N,mm))*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,ll))*Te(1+jj*N,ll));
                        end
                    end
                end
                for ll = 1:1+ii*N
                    if (ll ~= kk)
                        for mm = ll+1:1+ii*N
                            if mm ~= kk
                                tmp10 = tmp10 + 4*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,kk))^2*real(Te(1+ii*N,ll)*conj(Te(1+jj*N,ll))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                                tmp11 = tmp11 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,ll)*conj(Te(1+jj*N,ll))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                                tmp12 = tmp12 + 4*real(Tb(1+ii*N,kk)*conj(Tb(1+jj*N,kk))*conj(Tb(1+ii*N,mm))*Tb(1+jj*N,mm))*real(Te(1+ii*N,ll)*conj(Te(1+jj*N,ll))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                            end
                        end
                    end
                end
                for ll = 1:1+jj*N
                    if (ll ~= kk)
                        for mm = ll+1:1+ii*N
                            if mm ~= kk
                                tmp13 = tmp13 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,kk)*conj(Te(1+jj*N,kk))*conj(Te(1+ii*N,mm))*Te(1+jj*N,mm));
                            end
                        end
                        for mm = 1:1+ii*N
                            if (mm ~= kk) && (mm ~= ll)
                                for nn = mm+1 : 1+ii*N
                                    if (nn ~= kk) && (nn ~= ll)
                                        tmp14 = tmp14 + 2*abs(Tb(1+ii*N,kk))^2*abs(Tb(1+jj*N,ll))^2*real(Te(1+ii*N,mm)*conj(Te(1+jj*N,mm))*conj(Te(1+ii*N,nn))*Te(1+jj*N,nn));
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

H22E_bis_test = tmp1 + tmp2 + tmp3 + tmp4 +  tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp10 + tmp11 + tmp12 +  tmp13 + tmp14;


Hb2Hb2He2He2 = H4E_test + H22E_test + H22E_bis_test;



