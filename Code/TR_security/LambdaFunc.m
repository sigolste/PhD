function G_M_q = LambdaFunc(M, q, D)

% D: # expansion terms

temp = zeros(1,D-q+1);
d = 1;
for l = q:D-1
    if l==0 && q==0
        LahNumber_l_q = 1;
    elseif l>0 && q==0
        LahNumber_l_q = 0;
    elseif l>0 && q==1
        LahNumber_l_q = factorial(l);
    else
        LahNumber_l_q = nchoosek(l-1,q-1)*factorial(l)/factorial(q);
    end
    temp1   = (-1).^q.*sqrt(pi).*gamma(2*M).*gamma(1/2+l-M).*LahNumber_l_q;
    temp2   = 2.^(M-q).*gamma(1/2-M).*gamma(1/2+l+M).*factorial(l);
    temp(d) = temp1./temp2;
    d = d+1;
end

G_M_q = sum(temp);




% D: # expansion terms
% 
% temp = zeros(1,D-q+1);
% d = 1;
% for l = q:D-1
%     if l==0 && q==0
%         LahNumber_l_q = 1;
%     elseif l>0 && q==0
%         LahNumber_l_q = 0;
%     elseif l>0 && q==1
%         LahNumber_l_q = factorial(l);
%     else
%         LahNumber_l_q = nchoosek(l-1,q-1)*factorial(l)/factorial(q);
%     end
%     temp1   = (-1)^q*sqrt(pi)*gamma(2*M)*gamma(1/2+l-M)*LahNumber_l_q;
%     temp2   = 2^(M-q)*gamma(1/2-M)*gamma(1/2+l+M)*factorial(l);
%     temp(d) = temp1/temp2;
%     d = d+1;
% end
% 
% G_M_q = sum(temp);