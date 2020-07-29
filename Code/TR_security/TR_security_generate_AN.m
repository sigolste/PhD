function W = TR_security_generate_AN(H , nb_subcar , BOR, despread_matrix, type)
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
%
% OUTPUT:
%   W: The artificial noise vector
%
% By Sidney Golstein
% Last Update:  5 March 2019

N = nb_subcar/BOR;
switch type
    case "underdetermined"
        % Resolve the underdetermined system --> leads to peaks of power
        % To keep into account the BOR effect, as well as the despreading, 
        % the following steps are implemented:
        %   1. Create W, random variable with ,nb_subcar/BOR less elements
        %   2. Create G = H.*W, with [G] = [nb_subcar,1] 
        %   3. Split G in G1 and G2, with G1 determined and G2 to be found.
        %   4. Split S in S1 and S2, both determined
        %   5. Find the last elements of W thanks to S1 S2 G1 G2.
        %
        
        % Split the despread matrix
        S1 = despread_matrix(:,1:nb_subcar-nb_subcar/BOR);
        S2 = despread_matrix(:,nb_subcar-nb_subcar/BOR+1:end);


        W = randn(nb_subcar - nb_subcar/BOR ,1) + 1i*randn(nb_subcar - nb_subcar/BOR ,1);

        H_tmp = H(1:nb_subcar - nb_subcar/BOR,1);

        G1 = H_tmp.*W;
        G2 = H(nb_subcar - nb_subcar/BOR +1 : end ,1); 

        tmp = S1*G1;


        % Step 3
        to_add = -tmp./(S2*G2);
        W = [W;to_add];

        % Step 4
        W = reshape(W,[],1);

        energy_AN = 1/nb_subcar*sum(abs(W).^2,1);

%         W = W/sqrt(energy_AN*BOR);
        
    case "determined"
        % Resolve determined system to find W. W is random since the system
        % is made determined in a random fashion. If nb_subcar = Q,
        % nb_subcar/BOR = N ; the system is made determined by adding Q-N
        % rows in the matrix S^H H , and Q-N elements in the null vector
        % (B). The new matrix S^H H is denoted H_bis
        % Therefore, the rectangular system H_bis W = B can be solve via: 
        % W = H_bis \ B
        
        H = despread_matrix*diag(H);                                            % nb_subcar/BOR x nb_subcar
        tmp = sum(sum(H));
        mean_energy = 1/numel(H)*abs(tmp)^2;
        H_add = mean_energy.*( randn(nb_subcar-nb_subcar/BOR , nb_subcar) ...   
                + 1j*randn(nb_subcar-nb_subcar/BOR , nb_subcar) );              % nb_subcar-nb_subcar/BOR x nb_subcar
        H_bis = [H ; H_add];                                                    % nb_subcar x nb_subcar
        B = zeros(nb_subcar/BOR,1);
        B_add = mean_energy.*( randn(nb_subcar-nb_subcar/BOR , 1) ...
                + 1j*randn(nb_subcar-nb_subcar/BOR , 1));
        B = [B ; B_add];
        W = H_bis\B;       
        
    case "svd"
        A = despread_matrix*diag(H);
        [U,S,V] = svd(A);
        VH = ctranspose(V);
        %V1H = VH(1:N,:);
        V2H = VH(N+1:end,:);
        V2 = ctranspose(V2H);
        W_tilde = 1/sqrt(2)*(randn(nb_subcar-N,1) + 1j*randn(nb_subcar-N,1));
        W = V2*W_tilde;
end
end




  






