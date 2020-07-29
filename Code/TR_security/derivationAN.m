function [term_1 , term_19] = derivationAN(Hb,He , Q , despread_matrix, e)
        A = despread_matrix*diag(Hb);
        B = ctranspose(A)*A;
        [eigVectors eigValues] = eigs(B);
        [U,S,V] = svd(A);
        sigma = S(1:N,1:N); %
        sigmaH = ctranspose(sigma); %
        UH = ctranspose(U); %
        VH = ctranspose(V);
        V1H = VH(1:N,:); %
        V1 = ctranspose(V1H); %
        V2H = VH(N+1:end,:);
        V2 = ctranspose(V2H);
        
        
        W_tilde = 1/sqrt(2)*(randn(Q-N,1) + 1j*randn(Q-N,1));
        W_tildeH = ctranspose(W_tilde); %
        W = V2*W_tilde;
        
        %%%%% V2 = null(A) ; % CA PERMET DE PAS FAIRE LA SVD
        sigmacarre = sigma^2; %
        V1carre = V1*V1H; %
        V2carre = V2*V2H; %
        
%         V1 = eigVectors(:,Q-N+1:end);
% %         V1 = [conj(A(1,1))/conj(A(1,3)) 0 ; 0 conj(A(2,2))/conj(A(2,4)); 1 0 ; 0 1]; 
% 
%         
%         V2 = eigVectors(:,1:Q-N);
        for ii = 1:Q
            v1ih = V1(ii,:);
            v1i = ctranspose(v1ih);
            v2ih = V2(ii,:);
            v2i = ctranspose(v2ih);
            tmp1(ii) = v1ih*v1i;
            tmp2(ii) = v2ih*v2i;
            tmp3(ii) = tmp1(ii)^2;
            tmp4(ii) = tmp2(ii)^2;
       end
%        term_19 = tmp3;
        term_1 = 1/N*trace(U*sigma*V1H*abs(He).^2*V2*V2H*abs(He).^2*V1*sigmaH*UH); % OK bien egal à l energie AN apres decod 2

         for ii = 1:Q
            term_19(ii) = nonzeros(V1(ii,:));
         end
       
         
        % BOR = 2
        
        if Q/N == 2
            % Components taken by column
            %%%% V1
            COMP1 = conj(A(1,1))/sqrt(abs(A(1,1)).^2+abs(A(1,3)).^2);
            COMP2 = conj(A(1,3))/sqrt(abs(A(1,1))^2 + abs(A(1,3))^2);

            COMP3 = conj(A(2,2))/sqrt(abs(A(2,2))^2 + abs(A(2,4))^2);
            COMP4 = conj(A(2,4))/sqrt(abs(A(2,2))^2 + abs(A(2,4))^2);

            %%%% V2
            COMP5 = -A(2,4)/(sqrt(abs(A(2,4)/A(2,2))^2+1)*A(2,2));
            COMP6 = 1/(sqrt(abs(A(2,4)/A(2,2))^2+1));

            COMP7 = -A(1,3)/(sqrt(abs(A(1,3)/A(1,1))^2+1)*A(1,1));
            COMP8 = 1/(sqrt(abs(A(1,3)/A(1,1))^2+1));
            
        end
%       W = W*sqrt(e/energy(W));    
        
%         
        % term 1
% 
%         eq = zeros(Q,1);
%         term_2 = zeros(N,N);
%         
%         % term 2
%         for ii = 1:Q
%             eq(ii) = 1;
%             eqT = eq.';
%             tmp = U*sigma*V1H*eq*eqT*V2*V2H*eq*eqT*V1*sigmaH*UH;
%             term_2 = term_2+tmp;
%             eq = zeros(Q,1);
%         end
%         term_2 = 1/N*trace(term_2);
%         
%         
%         % term 3
%         term_3 = zeros(N,N);
%         ei = zeros(Q,1);
%         ej = zeros(Q,1);
%         for ii = 1:Q
%             ei(ii) = 1;
%             eiT = ei.';
%             for jj = 1:Q
%                 ej(jj) = 1;
%                 ejT = ej.';
%                 tmp = abs(He(ii,ii)).^2*abs(He(jj,jj)).^2*(U*sigma*V1H*ei*eiT*V2*V2H*ej*ejT*V1*sigmaH*UH);
%                 term_3 = term_3+tmp;
%                 ej = zeros(Q,1);
%             end
%             ei = zeros(Q,1);
%         end
%         term_3 = 1/N*trace(term_3);
%         
%         
%         % term 4
%         term_4 = zeros(N,N);
%         ei = zeros(Q,1);
%         ej = zeros(Q,1);
%         for ii = 1:Q
%             ei(ii) = 1;
%             eiT = ei.';
%             tmp1 = abs(He(ii,ii)).^4*(U*sigma*V1H*ei*eiT*V2*V2H*ei*eiT*V1*sigmaH*UH);
%             term_4 = term_4+tmp1;
%             for jj = 1:Q
%                 if jj ~= ii
%                     ej(jj) = 1;
%                     ejT = ej.';
%                     tmp2 = abs(He(ii,ii)).^2*abs(He(jj,jj)).^2*(U*sigma*V1H*ei*eiT*V2*V2H*ej*ejT*V1*sigmaH*UH);
%                     term_4 = term_4+tmp2;
%                     ej = zeros(Q,1);
%                 end
%             end
%             ei = zeros(Q,1);
%         end
%         term_4 = 1/N*trace(term_4);
%         
%         
%         
%         
%         term_5 = zeros(N,N);
%         ei = zeros(Q,1);
%         ej = zeros(Q,1);
%         for ii = 1:Q
%             ei(ii) = 1;
%             eiT = ei.';
%             tmp1 = 2*(U*sigma*V1H*ei*eiT*V2*V2H*ei*eiT*V1*sigmaH*UH);
%             term_5 = term_5+tmp1;
%             for jj = 1:Q
%                 if jj ~= ii
%                     ej(jj) = 1;
%                     ejT = ej.';
%                     tmp2 = (U*sigma*V1H*ei*eiT*V2*V2H*ej*ejT*V1*sigmaH*UH);
%                     term_5 = term_5+tmp2;
%                     ej = zeros(Q,1);
%                 end
%             end
%             ei = zeros(Q,1);
%         end
%         term_5 = 1/N*trace(term_5);
%         
%         
%         term_6 = zeros(N,N);
%         tmp2 = zeros(Q,Q);
%         ei = zeros(Q,1);
%         for ii = 1:Q
%             ei(ii) = 1;
%             eiT = ei.';
%             tmp = ei*eiT*V2*V2H*ei*eiT;
%             tmp2 = tmp2+tmp;
%             ei = zeros(Q,1);
%         end
%         term_6 = 1/N*trace(sigmacarre*V1H*tmp2*V1);
%         
%         
%         term_7 = 0;
%         ei = zeros(Q,1);
%         for ii = 1:Q
%             ei(ii) = 1;
%             eiT = ei.';
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             tmp = trace(sigmacarre*v1i*eiT*V2*V2H*ei*v1ih);
%             term_7 = term_7+tmp;
%             ei = zeros(Q,1);
%         end
%         term_7 = 1/N*term_7;
%         
%         
%         
%         term_8 = 0;
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp = trace(sigmacarre*v1i*v2ih*v2i*v1ih);
%             term_8 = term_8+tmp;
%         end
%         term_8 = 1/N*term_8;
%         
%         
%         
%         term_9 = 0;
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp = v2ih*v2i*v1ih*sigmacarre*v1i;
%             term_9 = term_9+tmp;
%         end
%         term_9 = 1/N*term_9;
%         
%         term_10 =0;
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp = v2ih*v2i*v1ih*v1i*mean(diag(sigmacarre));
%             term_10 = term_10+tmp;
%         end
%         term_10 = 1/N*term_10;




%         term_11 = 0;
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp1 = v2ih*v2i;
%             tmp2 = v1ih*v1i;
%             tmp3 = mean(diag(sigmacarre));
%             term_11 = term_11+tmp1*tmp2*tmp3;
%         end
%         term_11 = 1/N*term_11;
        
        
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             tmp(ii) = v1ih*v1i;
%         end
%         term_12 = mean(tmp);
%         
%         
%         for ii = 1:Q
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp(ii) = v2ih*v2i;
%         end
%         term_13 = mean(tmp);
        
        
%         term_14 = 0;
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp1 = v1ih*v1i;
%             tmp2 = v2ih*v2i;
%             term_14 = term_14+tmp1*tmp2;
%         end
%         term_14 = 1/N*term_14;
        
        
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp1(ii) = v1ih*v1i;
%             tmp2(ii) = v2ih*v2i;
%         end
%         term_15 = 1/N*Q*mean(tmp1.*tmp2);
%         
%         figure;
%         plot(tmp1,'r'); hold on;
%         plot(tmp2,'b'); hold on;
%         plot(tmp1+tmp2,'g') ; hold on;
%         plot(tmp1.*tmp2,'m');
% %         
%         
%         
%         for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp1(ii) = v1ih*v1i;
%             tmp2(ii) = v2ih*v2i;
%         end
%         term_16 = 1/N*Q*mean(tmp1)*mean(tmp2);
        
        
        
%         for ii = 1:Q
%             vih = V(ii,:);
%             vi = ctranspose(vih);
%             tmp = vih*vi;
% 
%         end
%         term_17 = 1/N*Q*mean(tmp);
        


%        
%        for ii = 1:Q
%             v1ih = V1(ii,:);
%             v1i = ctranspose(v1ih);
%             v2ih = V2(ii,:);
%             v2i = ctranspose(v2ih);
%             tmp1(ii) = v1ih*v1i;
%             tmp3(ii) = tmp1(ii)^2;
%        end
%        term_18 = 1/N*Q*mean(tmp1-tmp3);
        
      
%        term_19 = mean(abs(nonzeros(V1)).^4); %norm(v1ih)^4;
         
         
%          for ii = 1:Q
%             term_19(ii) = nonzeros(V1(ii,:));
%          end

       
%        figure;
%        plot(tmp1,'r','Marker','o'); hold on;
%        plot(tmp2,'b','Marker','o'); hold on;
%        plot(tmp1+tmp2,'g','Marker','o') ; hold on;
%        plot(tmp1.*tmp2,'m','Marker','o');
%        legend('v1^2','v2^2','v1^2 + v2^2', 'v1^2 * v2^2')
        

%         figure;
%         plot(tmp3,'b','Marker','o'); hold on;
%         plot(tmp4,'r','Marker','o'); hold on;
%         plot(tmp3+tmp4,'g','Marker','o')
%         legend('v1^4','v2^4','v1^4+v2^4')
%         xlabel('subcarrier')

%         test = V1H*V1;
%         test2 = test^2;
%         term_19 =   V1;
%         term_19 = W;
       
% 
% 
%        
%        for ii = 1:N
%            term_20(:,ii) = nonzeros(V1(:,ii));
% 
%        end
%        
       
              
%       term_19 =  mean(tmp3);%nonzeros(V1);
%        term_20 = mean(tmp3);%1/N*Q*(mean(tmp1)-mean(tmp3));
%         term_19 = 1/N*Q*(mean(tmp1)-mean(tmp3));
%        
%        figure;
%        plot(tmp1,'r','Marker','o'); hold on;
%        plot(tmp2,'b','Marker','o'); hold on;
%        plot(tmp1+tmp2,'g','Marker','o') ; hold on;
%        plot(tmp1.*tmp2,'m','Marker','o');
%        legend('v1^2','v2^2','v1^2 + v2^2', 'v1^2 * v2^2')
%        xlabel('subcarrier')