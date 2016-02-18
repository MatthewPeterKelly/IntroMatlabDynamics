% TEST_smoothingMethods.m
%
% This script demonstrates how two simple smoothing algorithms work
%
% Butterworth Filter:   (second-order)
%   - The butterworth filter is a standard digital filter used in signal
%   processing. It is recursive, meaning that the smoothed value at any
%   point is a function of all previous points. 
%   - A butterworth filter typically has a time lag, but it can be
%   converted into a smoothing algorith by simply running the filter
%   forward over the data, and then backward over the data, and then taking
%   the average of the two methods.
%   - In matlab, you can run a butterworth filter by using the commands:
%       >> help butter
%       >> help filtfilt
%
% Moving Average:  
%   - The moving average is a simple smoothing filter, where you just take
%   the average value of the data in some moving window.
%
% Smoothing parameter Wn:
%   - Both filters have a single parameter: the cut-off frequency (wn).
%   This parameter adjusts how much smoothing should take place. wn -> 0
%   corresponds to infinite smoothing: the average value of the data set.
%   In the limit as wn -> 1, there is no smoothing. The filtered data is
%   the same as the raw data.
%   - wn is normalized by the sample rate, as shown below. For example,
%   suppose that your data is taken at 1000 Hz, and you would like to
%   remove all frequency content above 50 Hz. Then you would calculate:
%   wn = 2*(50 Hz cut-off)*(0.001 sec sample period) = 0.1
%
% Boundaries:
%   Notice that both smoothing algorithms perform poorly near the
%   boundaries of the data set. This is because both algorithms are only
%   fully defined for the interior points. You must make some sort of
%   approximation at the boundaries. There are many fancy things to do
%   here, but I've implemented two relatively simple schemes. In the
%   butterworth filter, I assume that the filtered data matches the raw
%   data exactly at the boundary points, and then use a first-order filter
%   to compute the points at i=2 and i = (n-1). Then the remaining points
%   us a second order filter. In the moving average, I simply use fewer
%   points as the filter approaches the edges of the data.
%


clc; clear;

% Generate a random set of data to be smoothed
t = linspace(0,15,1000)';
xTruth = exp(-0.1*t).*sin(t);
xNoise = xTruth + 0.1*randn(size(t));

% Set the smoothing paramter:
wn = 0.02;   %  = 2*(cutoff frequency)*(sample period); 

% Smooth Data using each of the two methods:
yButter = butterworthSmoothing(xNoise,wn);   % Butterworth Smoothing
yMovAvg = movingAverage(xNoise,wn);  % Moving average

% Plot:
figure(1); clf; hold on;
plot(t,xTruth,'k-');
plot(t,xNoise,'b.');
plot(t,yButter,'r--');
plot(t,yMovAvg,'m:');
legend('truth','data','butterworth','moving average')
xlabel('time')
ylabel('value')
title('Comparison of Smoothing Filters')

