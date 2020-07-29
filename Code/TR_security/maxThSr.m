function maxSR_th = maxThSr(U , sr , alpha_opt , alpha_step )
% Return the theoretical max of SR for different BOR values
nb_step = 100/alpha_step + 1;
id = round(alpha_opt.*nb_step);

for bb = 1:length(U)
    maxSR_th(bb) = sr(id(bb),bb);
end