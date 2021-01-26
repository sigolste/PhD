function H2H6 = modelCorrelH2H6(N,U,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modelization of ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^6 when frequency 
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
%   H2H6 : Expected value of ∑_i ∑j!=i |Hb,n+iN|^2 |Hb,n+jN|^6
%
%
%
% Code started : 21.01.2021
% Last update  : 21.01.2021
%
% © Sidney Golstein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ∑_i ∑_j!=i |Hb,n+iN|^6 |Hb,n+jN|^2 = 
% ∑_i ∑_j>i  |Hb,n+iN|^6 |Hb,n+jN|^2 + ∑_i ∑_j>i |Hb,n+iN|^2 |Hb,n+jN|^6
%
% These two terms are computed separately in Sec 1. and 2. respectively


%% 1. ∑_i ∑_j>i |Hb,n+iN|^6 |Hb,n+jN|^2 

% Derivation of |Hb_n+iN|^6 cf CF 23 recto

% 1.1 First term of hb_i^6 (with h^6) multiplied by hb_j^2 -> term62
% Derivation cf CF23 verso

tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            tmp = tmp + 24*abs(T(1+ii*N,kk))^6*abs(T(1+jj*N,kk))^2;
            for ll = 1:1+jj*N
                if ll ~= kk
                    tmp = tmp + 6*abs(T(1+ii*N,kk))^6*abs(T(1+jj*N,ll))^2;
                end
            end
        end
    end
end
TERM62_test = tmp; % OK idem que TERM62


% Second term of hb_i^6 (with h^2 * h^4) multiply by hb_j^2 -> term242
% Derivation cf CF23 verso

% TERM242 = mean(term242);
tmp1 = 0;
tmp2 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            for ll = 1:1+ii*N
                if ll ~= kk
                    tmp1 = tmp1 + 4*abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^4*abs(T(1+jj*N,kk))^2 + 6*abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^4*abs(T(1+jj*N,ll))^2; %4
                    for mm = 1:1+jj*N
                        if (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
                            tmp2 = tmp2 + 2*abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^4*abs(T(1+jj*N,mm))^2;
                        end
                    end
                end
            end
        end
    end
     
end
tmp = tmp1 + tmp2;
TERM242_test = 9*tmp; %-> ok idem que TERM242

% Third term of hb_i^6 (with h^2 * h^2 * h^2) multiplied by hb_j^2 -> term2222
% Derivation cf CF24 recto

% TERM2222 = mean(term2222);
tmp1 = 0;
tmp2 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            for ll = 1:1+ii*N
                for mm = 1:1+ii*N
                    if (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
                        tmp1 = tmp1 + 2*(abs(T(1+ii*N,kk))^2*abs(T(1+jj*N,kk))^2*abs(T(1+ii*N,ll))^2*abs(T(1+ii*N,mm))^2 ...
                            + abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,ll))^2*abs(T(1+ii*N,mm))^2 ...
                            + abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,mm))^2);
                    end
                    for nn = 1:1+jj*N
                        if (nn ~= kk) && (nn ~= ll) &&  (nn ~= mm) && (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
                            tmp2 = tmp2 + abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^2*abs(T(1+ii*N,mm))^2*abs(T(1+jj*N,nn))^2;
                        end
                    end
                end
            end
        end
    end
end
tmp = tmp1 + tmp2;
TERM2222_test = 6*tmp;  %-> ok idem que TERM2222



% Fourth term of hb_i^6 (with h^4 * h * h^*) multiplied by hb_j^2 -> term41112
% Derivation cf CF24 recto

%TERM4112 = mean(term4112);
tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            for ll = 1:1+ii*N
                for mm = ll+1:1+ii*N
                    if (kk ~= ll) && (kk ~= mm)
                        tmp = tmp + abs(T(1+ii*N,kk))^4*real(T(1+ii*N,ll)*conj(T(1+ii*N,mm)) ...
                                                            *conj(T(1+jj*N,ll))*T(1+jj*N,mm));
                    end
                end
            end
        end
    end
end
TERM4112_test = 36*tmp;  %-> ok idem que TERM4112

% Fifth and seventh term of hb_i^6 with hb_j^2 will be zero

% Sixth term of hb_i^6 (with h^2 h^2 h h^*) multiplied by hb_j^2 -> term22112
% Derivation cf CF24 verso
% TERM22112 = mean(term22112);

tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1           % j > i
        for kk = 1:1+ii*N
            for ll =1:1+ii*N
                for mm = 1:1+ii*N
                    for nn = mm+1:1+ii*N
                        if (mm ~= ll) && (mm ~= kk) && (kk ~= ll) && (nn ~= ll) && (nn ~= kk) 
                            tmp = tmp + abs(T(1+ii*N,kk))^2*abs(T(1+ii*N,ll))^2 ...
                                *real(T(1+ii*N,mm)*conj(T(1+ii*N,nn)) ...
                                *conj(T(1+jj*N,mm))*T(1+jj*N,nn));
                        end
                    end 
                end 
            end
        end
    end
end
TERM22112_test = 72*tmp;


% --> ∑_i ∑_j>i |hb_{n+iN}|^6 |hb_{n+jN}]|^2 
% FINALH6H2 = mean(finalH6H2);
FINALH6H2_TEST = TERM62_test + TERM242_test + TERM2222_test + TERM4112_test + TERM22112_test;



%% 2. ∑_i ∑_j>i |Hb,n+iN|^2 |Hb,n+jN|^6 

% % First term of hb_i^2 multiplied by  hb_j^6 (with h^6) -> term26
% Derivation cf CF24 verso
% TERM26 = mean(term26); 

tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            tmp = tmp + 24*abs(T(1+jj*N,kk))^6*abs(T(1+ii*N,kk))^2;
        end
        for kk = 1:1+jj*N
            for ll = 1:1+ii*N
                if ll ~= kk
                    tmp = tmp + 6*abs(T(1+jj*N,kk))^6*abs(T(1+ii*N,ll))^2;
                end
            end
        end
    end
end
TERM26_test = tmp; % OK idem que TERM26


% % Second term of hb_i^2 multiplied by  hb_j^6 (with h^2 * h^4) -> term224
% Derivation cf CF25 recto
% TERM224 = mean(term224);
tmp1 = 0;
tmp2 = 0;
tmp3 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            for ll = 1:1+jj*N
                if ll ~= kk
                    tmp1 = tmp1 + 4*abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^4*abs(T(1+ii*N,kk))^2 ; %4
                end
            end
        end
        for kk = 1:1+jj*N
            for ll = 1:1+ii*N
                if ll ~= kk
                    tmp2 = tmp2 + 6*abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^4*abs(T(1+ii*N,ll))^2; %4
                end
            end
        end
        for kk = 1:1+jj*N
            for ll = 1:1+jj*N
                for mm = 1:1+ii*N
                    if (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
                        tmp3 = tmp3 + 2*abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^4*abs(T(1+ii*N,mm))^2;
                    end
                end
            end
        end
    end
     
end
tmp = tmp1 + tmp2 + tmp3;
TERM224_test = 9*tmp; %-> ok idem que TERM224


% % Third term of hb_i^2 multiplied by  hb_j^6 (with h^2 * h^2 * h^2) -> term2222_bis
% Derivation cf CF25 verso
% TERM2222_BIS = mean(term2222_bis);
tmp1 = 0;
% tmp2 = 0;
% tmp3 = 0;
tmp4 = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+ii*N
            for ll = 1:1+jj*N
                for mm = 1:1+jj*N
                    if (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
                        tmp1 = tmp1 + 2*abs(T(1+jj*N,kk))^2*abs(T(1+ii*N,kk))^2*abs(T(1+jj*N,ll))^2*abs(T(1+jj*N,mm))^2 ;
                    end
                end
            end
        end
%         for kk = 1:1+jj*N
%             for ll = 1:1+ii*N
%                 for mm = 1:1+jj*N
%                     if (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
%                         tmp2 = tmp2 + 2*abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^2*abs(T(1+ii*N,ll))^2*abs(T(1+jj*N,mm))^2 ;
%                     end
%                 end
%             end
%         end
%         for kk = 1:1+jj*N
%             for ll = 1:1+jj*N
%                 for mm = 1:1+ii*N
%                     if (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
%                         tmp3 = tmp3 + 2*abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+ii*N,mm))^2 ;
%                     end
%                 end
%             end
%         end
        for kk = 1:1+jj*N
            for ll = 1:1+jj*N
                for mm = 1:1+jj*N
                    for nn = 1:1+ii*N
                        if (nn ~= kk) && (nn ~= ll) &&  (nn ~= mm) && (mm ~= kk) && (mm ~= ll) &&  (kk ~= ll)
                            tmp4 = tmp4 + abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^2*abs(T(1+jj*N,mm))^2*abs(T(1+ii*N,nn))^2 ;
                        end
                    end
                end
            end
        end
    end
end
tmp = 3*tmp1 + tmp4; % + tmp2 + tmp3 + tmp4;
TERM2222_BIS_test = 6*tmp;  %-> ok idem que TERM2222_bis


% Fourth term of hb_i^2 multiplied by hb_j^6  (with h^4 * h * h^*)  -> term2411
% Derivation cf CF25 verso
% TERM2411 = mean(term2411);
tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1
        for kk = 1:1+jj*N
            for ll = 1:1+ii*N
                for mm = ll+1:1+ii*N
                    if (kk ~= ll) && (kk ~= mm)
                        tmp = tmp + abs(T(1+jj*N,kk))^4*real(T(1+jj*N,ll)*conj(T(1+jj*N,mm)) ...
                                                            *conj(T(1+ii*N,ll))*T(1+ii*N,mm));
                    end
                end
            end
        end
    end
end
TERM2411_test = 36*tmp;  %-> ok idem que TERM4112


% Fifth and seventh term of hb_i^2 with hb_j^6 will be zero

% Sixth term of hb_i^2 multiplied by hb_j^6 (with h^2 h^2 h h^*) -> term22211
% Derivation cf CF25 verso
% TERM22211 = mean(term22211);
tmp = 0;
for ii = 0:U-1
    for jj = ii+1:U-1           % j > i
        for kk = 1:1+jj*N
            for ll =1:1+jj*N
                for mm = 1:1+ii*N
                    for nn = mm+1:1+ii*N
                        if (mm ~= ll) && (mm ~= kk) && (kk ~= ll) && (nn ~= ll) && (nn ~= kk) 
                            tmp = tmp + abs(T(1+jj*N,kk))^2*abs(T(1+jj*N,ll))^2 ...
                                *real(T(1+jj*N,mm)*conj(T(1+jj*N,nn)) ...
                                *conj(T(1+ii*N,mm))*T(1+ii*N,nn));
                        end
                    end 
                end 
            end
        end
    end
end
TERM22211_test = 72*tmp;


% --> ∑_i ∑_j>i |hb_{n+iN}|^6 |hb_{n+jN}]|^2 
% FINALH2H6 = mean(finalH2H6);
FINALH2H6_TEST = TERM26_test + TERM224_test + TERM2222_BIS_test + TERM2411_test + TERM22211_test;


% --> ∑_i ∑_j!=i |hb_{n+iN}|^6 |hb_{n+jN}]|^2 
% H2H6 = mean(h2h6);
H2H6 = FINALH6H2_TEST + FINALH2H6_TEST;
