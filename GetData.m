function [results] = GetData(v,t,stim_start,stim_end,tstop,tacs_peak)


% v = matrix_somaV(idx,:);
% tacs = matrix_tacs(idx,:);
[~,index] = findpeaks(v,'MinPeakHeight',0);

spikes = t(index); % find all spikes
spikes_on = spikes(spikes>stim_start & spikes<stim_end); % spikes while tacs on
spikes_off = setdiff(spikes,spikes_on); % spikes while tacs off

% [~,index] = findpeaks(tacs); % in samples
% tacs_peak = t(index); % in ms
% tacs_peak(end) = []; % not a real peak - we delete it
tacs_T = mean(diff(tacs_peak)); % time period (ms)
tacs_T_sample = mean(diff(index(1:end-1))); % time period (samples)

% Calculate vector strength, Rayleigh statistics and spike rate while tACS on
x = sum(cos( 2*pi*(spikes_on - tacs_peak(1)) / tacs_T));
y = sum(sin( 2*pi*(spikes_on - tacs_peak(1)) / tacs_T));

VS_on = 1/length(spikes_on) * sqrt(x.^2 + y.^2);
Rayleigh_on = length(spikes_on) * VS_on.^2;
spike_rate_on = length(spikes_on)/(stim_end - stim_start)*1000; % Spike rate on computation (spikes per second)

clear x y
% Calculate vector strength, Rayleigh statistics and spike rate while tACS off
x = sum(cos( 2*pi*(spikes_off - tacs_peak(1)) / tacs_T));
y = sum(sin( 2*pi*(spikes_off - tacs_peak(1)) / tacs_T));

VS_off = 1/length(spikes_off) * sqrt(x.^2 + y.^2);
Rayleigh_off = length(spikes_off) * VS_off.^2;
spike_rate_off = length(spikes_off)/(tstop - (stim_end - stim_start))*1000; % Spike rate off computation (spikes per second)

% spikes_on_samples = spikes_on/dt;
spikes_off_before = spikes_off(spikes_off<stim_start);
spikes_off_after = spikes_off(spikes_off>stim_end);

% save data
results.spikes = spikes';
results.spikes_on = spikes_on';
results.spikes_off = spikes_off';
results.x = x;
results.y = y;
results.VS_on = VS_on;
results.Rayleigh_on = Rayleigh_on;
results.spike_rate_on = spike_rate_on;
results.VS_off = VS_off;
results.Rayleigh_off = Rayleigh_off;
results.spike_rate_off = spike_rate_off;
results.spikes_off_before = spikes_off_before';
results.spikes_off_after = spikes_off_after';

end

