function W = TR_security_generate_AN_waterfilling(H , nb_subcar , BOR, despread_matrix, alpha)
% Function that generates the AN vector (W) in order to implement security
% for the TR communication scheme. 
%
% W has no effect at the legitimate receiver position, but corrupt the
% correct detection at any eavesdropper position.
% The channel between TX and the legitimate RX is H. It is of dimensions 
% [nb_subcar,1]. We need to keep the despreading process into account such
% that, after despreading: S^H.*H.*W = 0
% 
%
% INPUTS
%   H: legeitimate receiver's channel --> [H] = [nb_subcar,1]
%   nb_subcar: number of OFDM subcarriers
%   BOR: Back Of Rate used for the communication
%   despread_matrix: Despreading matrix --> [S] = [nb_subcar/BOR,nb_subcar]
%   alpha
%
% OUTPUT:
%   W: The artificial noise vector
%
% By Sidney Golstein
% Last Update:  5 March 2019



        % Resolve the underdetermined system --> leads to peaks of power
        % To keep into account the BOR effect, as well as the despreading, 
        % the following steps are implemented:
        %   1. Create W, random variable with ,nb_subcar/BOR less elements
        %   2. Create G = H.*W, with [G] = [nb_subcar,1] 
        %   3. Split G in G1 and G2, with G1 determined and G2 to be found.
        %   4. Split S in S1 and S2, both determined
        %   5. Find the last elements of W thanks to S1 S2 G1 G2.
        %
        
        alpha_tilde = sqrt(ones(nb_subcar,1) - alpha); 
        
        alpha_tilde_1 = alpha_tilde(1:nb_subcar-nb_subcar/BOR,1);
        alpha_tilde_2 = alpha_tilde(nb_subcar-nb_subcar/BOR+1:end,1);
        
        % Split the despread matrix
        S1 = despread_matrix(:,1:nb_subcar-nb_subcar/BOR);
        S2 = despread_matrix(:,nb_subcar-nb_subcar/BOR+1:end);


        W = randn(nb_subcar - nb_subcar/BOR ,1) + 1i*randn(nb_subcar - nb_subcar/BOR ,1);

        H_tmp = H(1:nb_subcar - nb_subcar/BOR,1);

        G1 = H_tmp.*alpha_tilde_1.*W;
        G2 = H(nb_subcar - nb_subcar/BOR +1 : end ,1).*alpha_tilde_2; 

        tmp = S1*G1;


        % Step 3
        to_add = -tmp./(S2*G2);
        W = [W;to_add];

        % Step 4
        W = reshape(W,[],1);

        energy_AN = 1/nb_subcar*sum(abs(W).^2,1);

        W = W/sqrt(energy_AN*BOR);
        
  
        
end




  






