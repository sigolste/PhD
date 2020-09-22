function [W] = generateANMISO(N_TX,Hb , Q , U, despread_matrix, e , type)
% Function that generates the AN vector (W) in order to implement security
% for the TR MISO communication scheme
%
% W has no effect at the legitimate receiver position, but corrupt the
% correct detection at any eavesdropper position.
% The channel between TX and the legitimate RX is H. It's dimension is
% [Q,Q,N_TX]. We need to keep the despreading process into account such
% that, after despreading: S^H.*H.*W = 0
% 
%
% INPUTS
%   N_RX: number of TX antennas 
%   H: legitimate receiver's channel --> [H] = [Q,Q,N_TX]
%   Q: number of OFDM subcarriers
%   U: Back Of Rate used for the communication
%   despread_matrix: Despreading matrix --> [S] = [Q/U,Q]
%   e : energy of the artifical noise (AN) vector
%
% OUTPUT:
%   W: (AN) vector --> [W] = [Q,Q,N_TX]
%
% By Sidney Golstein
% Last Update:  19 September 2020
% Hb = diag(Hb);
N = Q/U;
switch type
    case "svd"
        A = despread_matrix*Hb;
        V2 = null(A);
        W_tilde = 1/sqrt(2)*(randn(Q*N_TX-N,1) + 1j*randn(Q*N_TX-N,1));
        W = V2*W_tilde;
        W = W*sqrt(e/energy(W));
%         W = zeros(Q*N_TX,1);
% %         for nn = 1:N_TX
%         A = despread_matrix*squeeze(Hb(:,:,nn));
%         V2 = null(A);
%         W_tilde = 1/sqrt(2)*(randn(Q-N,1) + 1j*randn(Q-N,1));
%         tmp = V2*W_tilde;
% %         size(tmp)
% %         e
% %         energy(tmp)
%         tmp = tmp*sqrt(e/energy(tmp));
%         W(:,:,nn) = diag(tmp);
        %%%%% V2 = null(A) ; % CA PERMET DE PAS FAIRE LA SVD
        end

end





  






