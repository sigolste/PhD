function e = energy(x)
% Function that returns the mean energy of x. If x is a vector, it returns
% the mean energy per element. If x is a diagonal matrix, it returns the
% mean energy per element in the diagonal

if size(x,2) == size(x,1)   %
    e = mean(abs(diag(x).^2));
else
    e = mean(abs(x).^2);
end
    