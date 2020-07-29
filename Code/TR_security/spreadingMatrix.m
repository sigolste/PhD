function [matrix_spread,matrix_despread] = spreadingMatrix(Q,N,U)

matrix_spread   = zeros( Q , N );
code_spread     = reshape( 2.*randi([0 1], 1, Q).' - 1 , [] , U )/sqrt(U);   % Diagonal elements of the spreading matrix. BPSK spreading code

for i = 1:size(code_spread,2)
    matrix_spread(i*N-N+1:i*end/U,:) = diag(code_spread(:,i));
end
matrix_despread = ctranspose(matrix_spread); 
