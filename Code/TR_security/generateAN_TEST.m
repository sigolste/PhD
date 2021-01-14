function [W,V1,S] = generateAN_TEST(Hb , Q , U, despread_matrix, e , type)
% Function that generates the AN vector (W) in order to implement security
% for the TR communication scheme. 
%
% W has no effect at the legitimate receiver position, but corrupt the
% correct detection at any eavesdropper position.
% The channel between TX and the legitimate RX is H. It is of dimensions 
% [Q,1]. We need to keep the despreading process into account such
% that, after despreading: S^H.*H.*W = 0
% 
%
% INPUTS
%   H: legeitimate receiver's channel --> [H] = [Q,1]
%   Q: number of OFDM subcarriers
%   U: Back Of Rate used for the communication
%   despread_matrix: Despreading matrix --> [S] = [Q/U,Q]
%   e : energy of the artifical noise (AN) vector
%
% OUTPUT:
%   W: (AN) vector
%
% By Sidney Golstein
% Last Update:  30 April 2020
Hb = diag(Hb);
N = Q/U;
switch type
    case "underdetermined"
        % Resolve the underdetermined system --> leads to peaks of power
        % To keep into account the U effect, as well as the despreading, 
        % the following steps are implemented:
        %   1. Create W, random variable with ,Q/U less elements
        %   2. Create G = H.*W, with [G] = [Q,1] 
        %   3. Split G in G1 and G2, with G1 determined and G2 to be found.
        %   4. Split S in S1 and S2, both determined
        %   5. Find the last elements of W thanks to S1 S2 G1 G2.
        %
        
        % Split the despread matrix
        S1 = despread_matrix(:,1:Q-Q/U);
        S2 = despread_matrix(:,Q-Q/U+1:end);


        W = randn(Q - Q/U ,1) + 1i*randn(Q - Q/U ,1);

        H_tmp = H(1:Q - Q/U,1);

        G1 = H_tmp.*W;
        G2 = H(Q - Q/U +1 : end ,1); 

        tmp = S1*G1;


        % Step 3
        to_add = -tmp./(S2*G2);
        W = [W;to_add];

        % Step 4
        W = reshape(W,[],1);

        W = W*sqrt(e/energy(W)); 

%         W = W/sqrt(energy_AN*U);
        
    case "determined"
        % Resolve determined system to find W. W is random since the system
        % is made determined in a random fashion. If Q = Q,
        % Q/U = N ; the system is made determined by adding Q-N
        % rows in the matrix S^H H , and Q-N elements in the null vector
        % (B). The new matrix S^H H is denoted H_bis
        % Therefore, the rectangular system H_bis W = B can be solve via: 
        % W = H_bis \ B
        
        H = despread_matrix*diag(H);                                            % Q/U x Q
        tmp = sum(sum(H));
        mean_energy = 1/numel(H)*abs(tmp)^2;
        H_add = mean_energy.*( randn(Q-Q/U , Q) ...   
                + 1j*randn(Q-Q/U , Q) );              % Q-Q/U x Q
        H_bis = [H ; H_add];                                                    % Q x Q
        B = zeros(Q/U,1);
        B_add = mean_energy.*( randn(Q-Q/U , 1) ...
                + 1j*randn(Q-Q/U , 1));
        B = [B ; B_add];
        W = H_bis\B;       
        W = W*sqrt(e/energy(W));
        
    case "svd"
        A = despread_matrix*diag(Hb);
        [U,S,V] = svd(A);
        S = S(1:N,1:N);         % Non zeroes singular values
%         [b cc d] = svd(A); 
%         ccc = cc(1:N,1:N);
%         cccc = energy(ccc);
        V2 = null(A);
        W_tilde = 1/sqrt(2)*(randn(Q-N,1) + 1j*randn(Q-N,1));
        W = V2*W_tilde;
        W = W*sqrt(e/energy(W));
        V1 = V(:,1:N);
        %%%%% V2 = null(A) ; % CA PERMET DE PAS FAIRE LA SVD

end
end




  






