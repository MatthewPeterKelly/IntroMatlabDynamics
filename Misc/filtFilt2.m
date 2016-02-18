function y = filtFilt2(x,wn)
% y = filtFilt2(x,wn)
%
% Implements a simple second-order butterworth smoothing filter.
%
% INPUTS:
%   x = vector of uniformly spaced data
%   wn = cut-off frequency.   0 < wn < 1    (1 == half sample rate)
%
% OUTPUTS:
%   y = smoothed data
%

error('This code has a bug. Still in development.');

%%%% Compute the butterworth coefficients:
q = sqrt(2);
c = tan(0.5*pi*(1-wn));

b0 = 1/(1 + q*c + c*c);
b1 = 2*b0;
b2 = b0;

a0 = 1;
a1 = -2*(c*c-1.0)*b0;
a2 = (1-q*c+c*c)*b0;

%%%% Set-Up for filtering
n = length(x);
y = zeros(size(x));

%%%% Forward Pass:
y(1:2) = x(1:2);
for i=3:n
    x0 = x(i);
    x1 = x(i-1);
    x2 = x(i-2);
    y1 = y(i-1);
    y2 = y(i-2);
    y0 = (x0*b0 + x1*b1 + x2*b2 - y1*a1 - y2*a2)/a0;
    y(i) = y0;
end

%%%% Backward Pass:
z = flip(y);
for i=3:n
    y0 = y(i);
    y1 = y(i-1);
    y2 = y(i-2);
    z1 = z(i-1);
    z2 = z(i-2);
    z0 = (y0*b0 + y1*b1 + y2*b2 - z1*a1 - z2*a2)/a0;
    z(i) = z0;
end
y = flip(z);

end