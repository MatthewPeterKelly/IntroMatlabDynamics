% TEST_filtFilt2.m
%
% This script tests the simple filtering algorithm:

% Generate a random set of data to be filtered
t = linspace(0,10,500);
x = sin(t) + 0.1*randn(size(t));

% Filter data:
wn = 0.1;
y = filtFilt2(x,wn);

% Plot:
figure(1); clf; hold on;
plot(t,x,'.');
plot(t,y);


