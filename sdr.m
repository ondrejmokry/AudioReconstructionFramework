function  [sdr_val] = sdr (u, v)
% function that computes SDR, where u is considered the original
% signal and v is the processed signal

sdr_val = 20*log10(norm(u)/norm(u-v));

end