function p = cproj(x, xq, delta)

    p = proj([real(x); imag(x)], [real(xq); imag(xq)], delta);
    p = p(1:length(x)) + 1i*p(length(x)+1:end);
    
end