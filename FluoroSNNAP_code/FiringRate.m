function [w,t] = FiringRate(x,fps,T)
% Given a set of spikes (frame # when spike occured), and the max number of frames,
% compute the inhomogeneous firing rate - 1/ISI with zero-order hold

% Convert frames to seconds

s = x/fps;
mx = length(x);
r = zeros(1,T);
t = 0:1/fps:T/fps-1/fps;

for i=2:mx
    r(x(i-1):x(i)) = 1/(s(i)-s(i-1));
end

% Gaussian filter to smooth
time = -4:.01:4;
f = exp(-(time.^2)/(2*1^2));
f = f/sum(sum(f));
w = conv(r,f,'same');
w = w*max(r)/max(w);