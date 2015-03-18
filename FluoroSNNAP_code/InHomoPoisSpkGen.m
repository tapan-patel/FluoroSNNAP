function SpkTime = InHomoPoisSpkGen(r)


% Takes in a firing rate r as a function of time t, and generate Inhomogeneous
% Poisson spike train. 

% Set deltaT
if(max(r)*.1>0.05)
    deltaT = 0.01;
else
    deltaT = 0.1;
end

SpkTime = find(r*deltaT>rand(1,length(r)));
