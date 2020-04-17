function proj = proj(x, xq, delta)
% PROJ perfoms projection of vector x onto the set of feasible
% solutions for the dequantization problem.
%
% Input arguments
%       x       vector of input signal
%       xq      quantized signal
%       delta   quantization step

overstep_above = (x - xq) > delta/2;
overstep_below = (xq - x) > delta/2;

proj = x;

proj(overstep_above) = xq(overstep_above) + delta/2 - eps;
proj(overstep_below) = xq(overstep_below) - delta/2 + eps;

end