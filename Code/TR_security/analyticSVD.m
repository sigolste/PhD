% [U S V] = svd(A). V can be found vie the eigenvalue decomposition of A^*A
% In particular, the eigenvectors of A^*A are the columns of V (up to a
% complex constant).
%
% A is a NxQ matrix

syms z1 z2 z3 z4

A_analytic = [z1 0 z2 0; 0 z3 0 z4];

B_analytic = ctranspose(A_analytic)*A_analytic;

[eigVectors,eigValues] = eig(B_analytic);