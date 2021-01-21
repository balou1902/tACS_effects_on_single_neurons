function [PLV, times] = getPLV(spikes,signal_tacs,start,ending,dt,indx)

% start or ending = in samples
% indx = to see if we want the PLV before tacs (1), during tacs (2), or
% after tacs (3)

switch indx
    case 1 % pre
        z = hilbert(signal_tacs(1:start,1));
        index_true= round(spikes/dt);
    case 2 % on
        z = hilbert(signal_tacs(start:ending,1));
        index_true= round(spikes/dt - start );
    case 3 % post
        z = hilbert(signal_tacs(ending:end,1));
        index_true= round(spikes/dt - ending );
end

phase_tacs = ((angle(z)));

tmp = phase_tacs(index_true); % radian
PLV = abs(mean(exp(1i*tmp)));


times = tmp;    


end
