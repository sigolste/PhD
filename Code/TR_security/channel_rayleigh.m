function y = channel_rayleigh(nb_subcar)
y = diag(1/sqrt(2)*(randn(nb_subcar,1) + 1j*randn(nb_subcar,1)));