function [rise,fall] = TransientKinetics(x,fps)
% Given a calcium transient, return rise and fall time
% First need to convert raw trace to a deltaF/F0

x0 = x(1);
dx_x0 = x;

%% Define time and upscale
t = 0:1/fps:length(x)/fps-1/fps;
tfine = linspace(0,t(end),10*length(t));
xfine = interp1(t,dx_x0,tfine);


%%Rise and fall time
% The rise time t1/2: the time between the onset of the spike initiation 
% and the half-peak response. The decay time t1/2: the time of half decay 
% of a single exponential fit of the recovery from the peak response to the baseline.
[xmax,idx_max] = max(xfine);

x_rise_half = xmax/2;
[~,idx_half] = min(abs(xfine(1:idx_max)-x_rise_half));

if(~isempty(idx_half))
    rise = tfine(idx_max-idx_half+1);
else
    rise = NaN;
end

% compute fall time - half-life of exponential decay fit
tfall = tfine(idx_max:end)-tfine(idx_max);
xfall = xfine(idx_max:end);
[fit,gof] = exponentialFit(tfall,xfall);
if(gof.rsquare>.9 && fit.b>.1)
    fall = log(2)/fit.b;
else
    xfall = xfall-min(xfall);
    x_fall_half = xfall(1)/2;
    [~,idx_fall] = min(abs(xfall-x_fall_half));
    if(~isempty(idx_fall))
        fall = tfall(idx_fall);
    else
        fall = NaN;
    end
end