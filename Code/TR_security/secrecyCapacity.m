function sr = secrecyCapacity(x,y)
% Reurns the secrecy rate value given two SINRs x and y
sr = log2(1+x)-log2(1+y);