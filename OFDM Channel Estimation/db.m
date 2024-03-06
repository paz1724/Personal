function x = db(x,a)
    x = 10 * log10(abs(x));
    
    if nargin == 2
        minVal = min(x);
        if minVal < 0 
            x = x - minVal;
        end
    end
end