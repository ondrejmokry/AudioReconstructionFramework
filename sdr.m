function  [sdr_val] = sdr (u, v)
% SDR computes SDR, where u is considered the original
% signal and v is the processed signal.
%
% Date: 29/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

sdr_val = 20*log10(norm(u)/norm(u-v));

end