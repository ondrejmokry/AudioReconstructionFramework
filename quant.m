function [xq, delta] = quant(x, w)
% QUANT perfoms signal quantization of the input signal. Number of the
% quantization levels is computed according to the wordlength.
%
% Input arguments
%       x       vector of input signal
%       w       wordlength

% quantization step
delta = 2^(-w+1); 

% mid-riser quantization
xq = sign(sign(x)+eps) .* delta .* (floor(abs(x)/delta)+1/2);

% fix the extreme cases for 1 or -1 back to the correct quantization level
xq(xq > 1) = xq(xq > 1) - delta;
xq(xq < -1) = xq(xq < -1) + delta;

end