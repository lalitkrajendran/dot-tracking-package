function y = nans(size1, size2, size3)
% function that returns a nan array of the given dimensions
    if nargin == 1
        size2 = size1;
        size3 = 1;
    elseif nargin == 2
        size3 = 1;
    end
    
    y = ones(size1, size2, size3) * NaN;    
end