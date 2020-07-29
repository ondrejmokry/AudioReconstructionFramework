function p = cproj(x, xq, delta)
% CPROJ perfoms projection of vector x onto the set of feasible
% solutions for the dequantization problem in the complex setting.
%
% Input arguments
%       x       (complex) vector of input signal
%       xq      quantized signal
%       delta   quantization step
%
% Date: 29/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

p = proj([real(x); imag(x)], [real(xq); imag(xq)], delta);
p = p(1:length(x)) + 1i*p(length(x)+1:end);
    
end