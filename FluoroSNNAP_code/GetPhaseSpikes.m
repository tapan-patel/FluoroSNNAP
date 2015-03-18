function phase = GetPhaseSpikes(spikes,tfinal)

spikes = sort(spikes); % Make sure spike times are in ascending order
phase = 2*pi*rand(1,tfinal);
k =0; % Current spike
numspikes = length(spikes);
try
    for t=1:tfinal
        if(sum(t==spikes))
            k = k+1;
        end
        
        switch k
            case 0
                phase(t) = 2*pi*(t-spikes(1))/(spikes(1));
            case numspikes
                %phase(t) = 0;
            otherwise
                phase(t) = 2*pi*(t-spikes(k))/(spikes(k+1)-spikes(k)) + 2*pi*k;
        end
    end
end
% phase = mod(phase,2*pi);