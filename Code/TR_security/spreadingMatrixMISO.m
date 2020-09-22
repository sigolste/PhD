function [matrix_spread,matrix_despread] = spreadingMatrixMISO(Q,Na,N,U)

len = Q*Na;
matrix_spread   = []; %zeros( Q , N );
code_spread     = reshape( 2.*randi([0 1], 1, len).' - 1 , [] , 1 )/sqrt(U);   % Diagonal elements of the spreading matrix. BPSK spreading code

for i = 1:U*Na
    matrix_spread = [matrix_spread ; diag(code_spread(N*(i-1)+1:N*(i-1)+N))]; %upsample(code_spread(:,i),U,i-1);%  i*N-N+1:i*end/U,:) = diag(code_spread(:,i));
end

matrix_despread = ctranspose(matrix_spread); 