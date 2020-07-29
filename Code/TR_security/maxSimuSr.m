function [maximum , alpha] = maxSimuSr(sr,alpha_step)

[maximum,idx] = max(sr);
alpha = idx*alpha_step/100;