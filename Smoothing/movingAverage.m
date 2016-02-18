function y = movingAverage(x,w)
% y = movingAverage(x,w)
%
% Implements a simple moving average smoothing algorithm
%
% INPUTS:
%   x = vector of uniformly spaced data
%   w = cut-off frequency.   0 < w < 1    (1 == half sample rate)
%
% OUTPUTS:
%   y = smoothed data
%

n = length(x);
k = round(1/(2*w));  %Half-width of the smoothing window

if n <= 2*k
   error('Filter frequency is too low! Need: length(x) > 2*round(1/w)+1');
end

y = zeros(size(x));
for i=1:n
    idxLow = i-k; if idxLow < 1, idxLow = 1; end
    idxUpp = i+k; if idxUpp > n, idxUpp = n; end
    idx = idxLow:idxUpp;
    y(i) = mean(x(idx));
end

end