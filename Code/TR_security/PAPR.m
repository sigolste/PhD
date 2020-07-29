function res = PAPR(signal)

signal_peak = max(signal.*conj(signal));
signal_rms = mean(signal.*conj(signal));
res = 10*log10(signal_peak/signal_rms);