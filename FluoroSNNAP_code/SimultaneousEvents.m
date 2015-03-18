function c=SimultaneousEvents(x,y,lag)
% Given two spike trains x and y, the number of times that an event appears
% in x shortly after it appears in y within x specified time lag
c = 0;

mx = length(x); % Number of spikes in x
my = length(y);  % Number of spikes in y

% For x given spike at t=t0 in y, figure out if x spike occured in x within
% t=t0 to t=t0+lag. If the spike occured at t=t0 in both x and y, c=c+.5,
% else c= c+1 
for i=1:my
    for j=1:mx
        if(x(j)-y(i) <=lag &  x(j)-y(i) >0)
            c = c+1;
        elseif(x(j)==y(i))
            c = c+.5;
        end
    end
end