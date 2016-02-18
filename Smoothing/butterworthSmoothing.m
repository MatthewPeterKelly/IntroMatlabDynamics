function y = butterworthSmoothing(x,w)
% y = butterworthSmoothing(x,w)
%
% Implements a simple second-order butterworth smoothing filter. The
% boundary conditions are implemented by using a first-order filter, and
% then by using the raw data value at the very edge.
%
% INPUTS:
%   x = vector of uniformly spaced data
%   w = cut-off frequency.   0 < w < 1    (1 == half sample rate)
%
% OUTPUTS:
%   y = smoothed data
%

% Forward Pass: (output lags data)
yForward = forwardFilterPass(x,w);

% Backward Pass:  (output leads data)
yBackward = flip(forwardFilterPass(flip(x),w));

% Combine:   (output matches phase with data);
y = 0.5*(yForward + yBackward);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = forwardFilterPass(x,w)

n = length(x);
if n < 3
    error('Insufficient Data! Need length(x) >= 3')
end

%%%% Compute the butterworth coefficients:

% Useful constants:
q = sqrt(2);
c = tan(0.5*pi*(1-w));

% Second-Order
b0 = 1/(1 + q*c + c*c);
b1 = 2*b0;
b2 = b0;

a0 = 1;
a1 = -2*(c*c-1.0)*b0;
a2 = (1-q*c+c*c)*b0;

% First-Order:
c0 = 1/(1+c);
c1 = c0;

d0 = 1.0;
d1 = -(1-2*c0);

%%%% Run filter:
y = zeros(size(x));
y(1) = x(1);
y(2) = (x(2)*c0 + x(1)*c1 - y(1)*d1)/d0;
for i=3:n
    x0 = x(i);
    x1 = x(i-1);
    x2 = x(i-2);
    y1 = y(i-1);
    y2 = y(i-2);
    y0 = (x0*b0 + x1*b1 + x2*b2 - y1*a1 - y2*a2)/a0;
    y(i) = y0;
end

end

