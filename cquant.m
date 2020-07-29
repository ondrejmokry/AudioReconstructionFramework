function [xq, delta] = cquant(x, w)
% CQUANT perfoms signal quantization of the input signal in the complex
% setting. Number of the quantization levels is computed according to the
% wordlength w.
%
% Input arguments
%       x       vector of input signal
%       w       wordlength
%
% Date: 16/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% ensure the vector to be quantized attains values from -1 to 1
maxval = max(abs([real(x); imag(x)]));
x      = x/maxval;

% quantize the vector
[ xq, delta ] = quant([real(x); imag(x)], w);
xq            = xq(1:length(x)) + 1i*xq(length(x)+1:end);

% return values corresponding to the original range
xq    = maxval * xq;
delta = maxval * delta;

end