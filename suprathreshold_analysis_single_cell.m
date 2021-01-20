clear all
close all
clc

% Description:
%   Code to analyze the somatic membrane potential of a single cell when
%   stimulated by tACS
% code created in June 2020


% Correspondance: htran@umn.edu

%% -------------------------------------------------------------------------
clear all
close all
clc

% We create a vector of different tACS amplitudes
% Realistic amplitudes in humans: 0.1 - 3 mV/mm (Can be higher for
% non-human primates)

amplitude = [0:0.1:2 2.2:0.2:5 6.25:1.25:10 12.5:2.5:25];
Namp = length(amplitude);

% we load the different data 
for i=1:Namp
    cd(['tacs_' num2str(i, '%.1d')])
    somaV(i,:)= load('somaV.txt');    
    if i==Namp
        t=load('t.txt');
        dt = t(2) - t(1);
        tacs = load('TACS.txt');
        tacs = tacs./max(tacs); % we normalize
    end   
    cd('..\')
end
disp('All data loaded.')

% ---- tACS parameters
DEL = 120000; % ms
DUR = 2*DEL; % ms 
tstop = 4*DEL; % ms
tacs_start = DEL; % ms
tacs_end =  DEL+DUR; % ms
tacs_start_sample = DEL/dt; % samples
tacs_end_sample = (DEL + DUR)/dt; % samples
freq = 10;

% Creation of a sham tACS
sham_tacs =  sin(2*pi*freq*t/1000);
[~, idx] = findpeaks(tacs, 'MinPeakHeight', 0);

ShamPeaks = t(idx); % tacs peaks in ms 
tacs_peak = ShamPeaks(ShamPeaks>tacs_start & ShamPeaks<tacs_end);

%% 

% Computation of number of spikes
for i=1:Namp   
    clear idx
    [~, idx] = findpeaks(somaV(i,:), 'MinPeakHeight', 0);
    tmp = t(idx); % tacs peaks in seconds
    spikes{i,1} = tmp(tmp<tacs_start ); % PRE
    spikes{i,2} = tmp(tmp>tacs_start & tmp<tacs_end); % TACS ON    
    spikes{i,3} = tmp(tmp>tacs_end ); % POST
    
end

% Computation of the firing rate
for i=1:Namp
    spikes_rate(i,1) = length(spikes{i,2})/DUR*1000; % TACS ON
    spikes_rate(i,2) = (length(spikes{i,1})+length(spikes{i,1}))/DUR*1000; % TACS OFF (pre + post)
end

% we plot the firing rate with the amplitude
figure('color','w')
plot(amplitude,spikes_rate(:,1),'o--')
hold on
plot(amplitude,spikes_rate(:,2),'o--')
grid minor
legend('tacs on','tacs off', 'location', 'northwest')
xticks(amplitude)
xlabel('EF amplitude (V/m)')
title('Firing rate')
ylabel('spikes per second')
xlim([0 5])

% Computation of the phase-locking value (PLV) (
for i=1:Namp
     PLV(i,1) = getPLV(spikes{i,1},sham_tacs,tacs_start_sample,tacs_end_sample,dt,1);
     [PLV(i,2) times{i,1}] = getPLV(spikes{i,2},sham_tacs,tacs_start_sample,tacs_end_sample,dt,2);
     PLV(i,3) = getPLV(spikes{i,3},sham_tacs,tacs_start_sample,tacs_end_sample,dt,3);
end


%% Computation of the ISI - InterSpike Interval

x=[];
g=[];
clear isi

for i=1:36
    isi{i,1} = diff(spikes{i,2});
    x = [ x;  isi{i,1} ];
    g = [g ; amplitude(i)*ones(length(isi{i,1}), 1)];
end

figure('color','w')
boxplot(x, g,'notch','on')
grid minor
ylabel(' ms')
xlabel(' EF amplitude (V/m)')
title('ISI - tacs on')
xlabel([0 5])


% we compare the ISI for one amplitude
% we choose which amplitude to study
idx_amp = 5;

spk_off = [diff(spikes{idx_amp,1}) ; diff(spikes{idx_amp,3})];
spk_on = diff(spikes{idx_amp,2});

x = [spk_off; spk_on];
g = [zeros(length(spk_off), 1); ones(length(spk_on), 1)];

figure('color','w')
boxplot(x,g,'notch','on','labels',{'tacs OFF', 'tacs ON'});
grid minor
ylabel(' time (ms)')
xlabel(' EF amplitude (V/m)')
title(['Amplitude of ' num2str(amplitude(idx_amp), '%.1f') ' V/m'])


%% Additional statistic computation: rayleigh and Vector Strength

for i=1:Namp
    [results{i,1}] = GetData(somaV(i,:),t,tacs_start,tacs_end,tstop,tacs_peak);
    Rayleigh(i,1) = results{i,1}.Rayleigh_on;
    VS(i,1) = results{i,1}.VS_on;  
end

figure('color','w')
subplot(2,1,1)
plot(amplitude,VS,'o--')
grid minor
title('Vector Srength')
subplot(2,1,2)
plot(amplitude,Rayleigh,'o--')
grid minor
xlabel(' EF amplitude (V/m)')
title('Rayleigh statistics')



%% Computation of polar histograms for a specific amplitude

% Requirement: Circular Statistics Toolbox 
% link to download the toolbox: 
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

% all the 5 cells
figure('color','w')
for i=1:5
   subplot(1,5,i)
   polarhistogram(mat_times{1,i}{idx_amp},36)
   title(['Cell #' num2str(i) ' - PLV: ' num2str(round(PLV(idx_amp,i),2)) ])
end
suptitle(['Phase histogram - amplitude: ' num2str(amplitude(idx_amp)) 'V/m'])

% one specific cell
idx_cell = 5
figure('color','w')
polarhistogram(mat_times{1,idx_cell}{idx_amp},36)
title(['Cell #' num2str(idx_cell) ' - PLV: ' num2str(round(PLV(idx_amp,i),2))...
       ' - amplitude:'  num2str(amplitude(idx_amp)) 'V/m' ])

% to get the mean angle
alpha = mat_times{1,idx_cell}{idx_amp,1};
alpha_bar = circ_mean(alpha) % mean alpha
R = circ_r(alpha) % mean length
S = circ_var(alpha) % the spread in a data set.

% ================ testing for circular uniformity 

% 1) Rayleigh test
% H0: the population is distributed uniformly around the circle
% H1: the population is not distributed uniformly around the circle

p = circ_rtest(alpha)

% small p indicates a significant departure from uniformity and indicates
% to reject the null hypothesis

% 2) Omnibus test
alpha_deg = circ_rad2ang(alpha);
p = circ_otest(alpha)

% 3) Rao test
p = circ_raotest(alpha)
%For the examples, we find that, at the 0:05 significance level, the null hypothesis cannot be
%rejected for either sample (P > 0.05).
