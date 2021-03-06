function sinr = sinrModelingFrequencyCorrelation(alpha,U,N,Tb,Te,snr_b,snr_e,type)

sigma_b = 1./U/10^(snr_b/10);    % expected noise energy @Bob
sigma_e = 1./U/10^(snr_e/10);    % expected noise energy @Eve
% sigma_an = 1./U;

alpha = alpha.';
switch type 
    case "bob_correl"
        H2H2 = modelCorrelH2H2(N,U,Tb);
        H4 = modelCorrelH4(U);
        e_sym = alpha./(U^2)*(H4 + H2H2);
        e_noise = sigma_b;
        sinr = e_sym./e_noise;
        
    case "eve_decod1_correl"
        sinr = alpha./(U.*sigma_e + (1-alpha));
        
    case "eve_decod1_correlE_correlB"   
        %Correlation at Bob and at Eve, SDS decoder : cfr FC39 recto
        HbHbHeHe = modelCorrelHbHbHeHe(N,U,Tb,Te);
        e_sym = alpha./(U^2)*(U + HbHbHeHe);
        e_noise = sigma_e;
        e_an = (1-alpha)./U;
        sinr = e_sym./(e_noise + e_an);
        
    case "eve_decod2_correlE_correlB"
        %Correlation at Bob and at Eve, MF decoder : cfr FC39 verso - FC
        Hb2Hb2He2He2 = modelCorrelHb2Hb2He2He2(N,U,Tb,Te); % Derivation FC40/FC41
        e_noise = sigma_e;
        e_sym = alpha./(U^2)*(4*U + Hb2Hb2He2He2);
         % TERM AN:
        tmp1 = 0;    % |t_ik|^2 |t_jk|^2 |h_k|^4
        tmp2 = 0;    % |t_ik|^2 |t_jl|^2 |h_k|^2 |h_l|^2
        tmp3 = 0;         % t_ik t_il^* t_jk t_jl^* |h_k|^2 |h_l|^2

        for ii = 0: U-1
            for jj = ii+1 : U-1
                for kk = 1 : 1+ii*N
                     %[{1+ii*N,kk}  {1+jj*N,kk}] % tuples
                    tmp1 = tmp1 + 2*abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,kk))^2;
                end
            end
        end
        for ii = 0: U-2
            for jj = ii+1 : U-1
                for kk = 1 : 1+ii*N
                    for ll = 1 : 1+jj*N
                        if ll ~= kk
                            %[{1+ii*N,kk}  {1+jj*N,ll}] % tuples
                            tmp2 = tmp2 + abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,ll))^2   ;
                        end   
                    end
                end
            end
        end

        for ii = 0: U-2
            for jj = ii+1 : U-1
                for kk = 1 : 1+ii*N
                    for ll = kk+1 : 1+ii*N
                        %[{1+ii*N,kk}  {1+ii*N,ll} , {1+jj*N,kk} {1+jj*N,ll}] % tuples
                        tmp3 = tmp3 + 2*real(conj(Tb(1+ii*N,kk)) * Tb(1+ii*N,ll) * Tb(1+jj*N,kk) * conj(Tb(1+jj*N,ll)))   ;
                    end
                end
            end
        end


        sum_double_square_an = 2*(tmp1 + tmp2 + tmp3);                  % ∑_i ∑_{j!=i}  |h_i|^2 |h_j|^2
        sum_fourth_an = 2*U;
        
        v1_2 = 1./U;
        v1_4 = 2./U.^2 ./ (1./(U.^2).*(sum_fourth_an + sum_double_square_an));
        
        e_an = (1-alpha).*U./(U-1).*(v1_2 - v1_4 );%1./(U+1);%*
        sinr = e_an;%e_sym./(e_noise+e_an);
        
        
    case "eve_decod2_correl"
        % Consideration: normalisation of decod 2 by sqrt((U+1)/(U+3)). 
        % Else, multiply all by (U+3)/(U+1). Also, we consider that energy 
        % of an signal is 1/U
        
        H2H2 = modelCorrelH2H2(N,U,Tb);          % ∑_i ∑_{j!=i}  |h_i|^2 |h_j|^2
        H4 = modelCorrelH4(U);                  % ∑_i |h_i|^4 
        e_sym = alpha./U^2.*(H4+H2H2);
        e_noise = sigma_e;
        
        % TERM AN:
        tmp1 = 0;    % |t_ik|^2 |t_jk|^2 |h_k|^4
        tmp2 = 0;    % |t_ik|^2 |t_jl|^2 |h_k|^2 |h_l|^2
        tmp3 = 0;         % t_ik t_il^* t_jk t_jl^* |h_k|^2 |h_l|^2

        for ii = 0: U-1
            for jj = ii+1 : U-1
                for kk = 1 : 1+ii*N
                     %[{1+ii*N,kk}  {1+jj*N,kk}] % tuples
                    tmp1 = tmp1 + 2*abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,kk))^2;
                end
            end
        end
        for ii = 0: U-2
            for jj = ii+1 : U-1
                for kk = 1 : 1+ii*N
                    for ll = 1 : 1+jj*N
                        if ll ~= kk
                            %[{1+ii*N,kk}  {1+jj*N,ll}] % tuples
                            tmp2 = tmp2 + abs(Tb(1+ii*N,kk))^2 * abs(Tb(1+jj*N,ll))^2   ;
                        end   
                    end
                end
            end
        end

        for ii = 0: U-2
            for jj = ii+1 : U-1
                for kk = 1 : 1+ii*N
                    for ll = kk+1 : 1+ii*N
                        %[{1+ii*N,kk}  {1+ii*N,ll} , {1+jj*N,kk} {1+jj*N,ll}] % tuples
                        tmp3 = tmp3 + 2*real(conj(Tb(1+ii*N,kk)) * Tb(1+ii*N,ll) * Tb(1+jj*N,kk) * conj(Tb(1+jj*N,ll)))   ;
                    end
                end
            end
        end


        sum_double_square_an = 2*(tmp1 + tmp2 + tmp3);                  % ∑_i ∑_{j!=i}  |h_i|^2 |h_j|^2
        sum_fourth_an = 2*U;
        
        v1_2 = 1./U;
        v1_4 = 2./U.^2 ./ (1./(U.^2).*(sum_fourth_an + sum_double_square_an));
        
         e_an = (1-alpha).*U./(U-1).*(v1_2 - v1_4 );%1./(U+1);%*
        
        sinr = e_sym./(e_noise + e_an);
        
     case "eve_decod5"
        e_sym_decod5_model_e    = 2*alpha./U;
        e_noise_decod5_model_e  =sigma_e;
        e_an_decod5_model_e     = 2./U.*(1-alpha);
        
        sinr = e_sym_decod5_model_e./(e_noise_decod5_model_e + e_an_decod5_model_e);
end