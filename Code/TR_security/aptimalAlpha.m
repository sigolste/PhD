function alpha_opt = optimalAlpha(T1,T2,T3,T4)
% Compute the maximum of SR if the SR has the form:
% SR = (-alpha^2 T1 + alpha T2 + T3)/(alpha T4 + T3)
% This is the case for decod1 and decod2. Only T1 T2 T3 T4 change

alpha_opt = (sqrt(T1^2*T3^2+T1*T2*T3*T4-T1*T3*T4^2)-T1*T3)/(T1*T4);