clear all
close all
clc

% description:
%   Code to analyze the somatic membrane potential of all the cells belonging to a same cell type when
%   stimulated by tACS
% code created in July 2020


% correspondance: htran@umn.edu

%% -------------------------------------------------------------------------

amplitude = [0:0.1:2 2.2:0.2:5 6.25:1.25:10 12.5:2.5:25];
Namp = length(amplitude);

% we have 5 different cell types
L1GNC = [1:5];
L23PC = [6:10];
L4LBC = [11:15];
L5PC = [16:20];
L6PC = [21:25];

% we load all the data from the different folders
% change here the name of cell type and the matlab path
tic
for i=1:5
       cd(['D:\results_cell' num2str(L6PC(i), '%.1d') ])
    disp('==================')
    disp(['Cell loaded: ' num2str(i)])        
    for j=1:Namp
        j
        cd(['tacs_' num2str(j, '%.1d')])
        if j==1
            t=load('t.txt');
            dt = t(2) - t(1);
            tacs = load('TACS.txt');
            tacs = tacs./max(tacs); % we normalize
        end        
        tmp_soma(j,:)= load('somaV.txt');
        cd('..\')
    end    
    somaV{1,i} = tmp_soma;    
    cd('..\')
    clear tmp_soma    
end
toc

% ---- tACS parameters
DEL = 120000; % ms
DUR = 2*DEL; % ms 
tstop = 4*DEL; % ms

tacs_start = DEL; % ms
tacs_end =  DEL+DUR; % ms
tacs_start_sample = DEL/dt; % samples
tacs_end_sample = (DEL + DUR)/dt; % samples
freq = 10;
T = 1/freq; %second 

% creation of a sham tACS
sham_tacs =  sin(2*pi*freq*t/1000);
[~, idx] = findpeaks(tacs, 'MinPeakHeight', 0);
ShamPeaks = t(idx); % tacs peaks in ms 

tacs_peak = ShamPeaks(ShamPeaks>tacs_start & ShamPeaks<tacs_end);

%% 

% computation of number of spikes and the firing rate
for j=1:5
    tmpV = somaV{1,j};
    
    % number of spikes
    for i=1:Namp
        clear idx
        [~, idx] = findpeaks(tmpV(i,:), 'MinPeakHeight', 0);
        tmp = t(idx); % tacs peaks in seconds        
        spikes{i,1} = tmp(tmp<tacs_start ); % PRE
        spikes{i,2} = tmp(tmp>tacs_start & tmp<tacs_end); % TACS ON
        spikes{i,3} = tmp(tmp>tacs_end ); % POST
    end
    
    % firing rate
    for i=1:Namp
        spikes_rate(i,1) = length(spikes{i,2})/DUR*1000; % TACS ON
        spikes_rate(i,2) = (length(spikes{i,1})+length(spikes{i,3}))/DUR*1000; % TACS OFF (pre + post)
    end   
    SR{1,j} = spikes_rate;
    SPK{1,j} = spikes;
    clear spikes spikes_rate
    
end

% we concatenate the data for a better visualization
FR=[];
for i=1:5
    FR = [FR  SR{1,i}(:,1)];
end

figure('color','w')
plot(amplitude(1:46),FR' - FR(1,:)'*ones(1,46),'o--')
plot(amplitude(1:Namp),FR','o--')
xticks(amplitude)
xlabel('EF amplitude (V/m)')
title('Normalized firing rate')
ylabel('\Delta')
legend('cell #1', 'cell #2','cell #3','cell #4','cell #5','location','southeast')
grid minor

% linear regression on the firing rate
for icell=1:5
slope_FR(icell,:) = polyfit(amplitude',FR(1:46,icell),1);
end

% linear regression on the PLV
for icell=1:5
slope_PLV(icell,:) = polyfit(amplitude',PLV(1:46,icell),1);
end

figure('color','w')
scatter(slope_FR(:,1), slope_PLV(:,1),'filled')
grid minor
xlabel('Slope of firing rate')
ylabel('Slope of PLV')


%% we compute the ISI during the tACS and when tACS is off

for icell = 1:5
    tmp = SPK{1,icell};
    ISI{1,icell} = cellfun(@diff,tmp,'UniformOutput',false);
end

% here, you can change which cell and which amplitude you want to
% investigate and plot

icell = 5
iamp = 21
edges = [0:5:350]; % change here for a better visualization

figure('color','w')
hold on
histogram(ISI{1,icell}{iamp,1},edges)%,'Normalization','probability')
histogram(ISI{1,icell}{iamp,2},edges)%,'Normalization','probability')
histogram(ISI{1,icell}{iamp,3},edges)%,'Normalization','probability')
grid minor
% legend('pre', 'on', 'post')
xlabel('time (ms)')
ylabel('Occurences')
title(['Cell #' num2str(icell) ' - ' num2str(amplitude(iamp), '%.1f') ' V/m'])

clear mat_tmp
for icell=1:5
    for iamp=1:44
        mat_tmp{iamp,1} = [ ISI{1,icell}{iamp,1} ; ISI{1,icell}{iamp,3} ]; % OFF
        mat_tmp{iamp,2} = ISI{1,icell}{iamp,2}; % ON
    end
    isi_vs{1,icell} = mat_tmp;
    
end

%% we plot all the cell for each amplitude

clear V viON viOFF
edges = [0:5:350];
iamp =1;

for icell=1:5
    clear von voff
    figure('color','w','units','normalized','outerposition',[0 0.5 0.45 0.5])
    hold on
    % tACS OFF
    hoff = histogram(isi_vs{1,icell}{iamp,1},edges,'Normalization','probability')
    voff = hoff.Values;
    % tACS ON
    hon = histogram(isi_vs{1,icell}{iamp,2},edges,'Normalization','probability')
    BinEdges = hon.BinEdges;
    von = hon.Values;
    grid minor
    legend('OFF', 'ON')
    xlabel('ms')
    ylabel('Occurences')
    
    viON(icell,:) = von;
    viOFF(icell,:) = voff;
    V(icell,:) = von-voff;   
end

%% 1/ we plot all the amplitudes for each cell - WE PLOT THE ISI DIFFERENCE 

edges = [0:4:350]; % ms
id_amp = [1 6 11 16 21];
clear V 

% change here the color you want (comment/uncomment)
color = 'b'; % L1 GNC
% color =[0.3010, 0.7450, 0.9330];
% color =[0.4660, 0.6740, 0.1880];
% color =[0.4940, 0.1840, 0.556]; % L5PC
% color =[0.8500, 0.3250, 0.0980];


for icell = 1:5
    for iamp=1:length(id_amp)
        clear von voff tmp_on tmp_off
        % OFF
        hoff = histogram(isi_vs{1,icell}{id_amp(iamp),1},edges,'Normalization','probability')
        voff = hoff.Values;
        % ON
        hon = histogram(isi_vs{1,icell}{id_amp(iamp),2},edges,'Normalization','probability')
        BinEdges = hon.BinEdges;
        von = hon.Values;
        
        V(iamp,:) = von - voff;
    end
    close
    h=figure('color','w','units','normalized','outerposition',[0 0 1 1])
    for iamp=1:5
        subplot(5,1,iamp)
        bar(BinEdges(2:end),V(iamp,:),'FaceColor',color)
        grid minor
        
        if iamp == 1
            y = ylabel('0 mV/mm')
        elseif iamp ==2
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        elseif iamp ==3
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        elseif iamp ==4
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        elseif iamp ==5
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        end
        
        set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        set(get(gca,'YLabel'),'Rotation',0)
        set(gca,'xticklabel',{[]})
        set(gca,'linewidth',1,'fontweight','bold','fontsize',15);
        ylim([-0.02 0.015])
        
        if iamp==5
            xticklabels({'0','50','100','150','200','250','300','350'})
            xticks([0 50 100 150 200 250 300 350])
            xlabel('ms')
        end
        
    end
    suptitle(['L1 NGC - cell #' num2str(icell)])
    print(h,sprintf('L1_NGC_cell_%d.png',icell),'-dpng','-r500');
    
end


% we plot the ISI
% change here the path
cd('D:\ISI_allCells\Per cell\ISI_full_distribution')
id_amp = [1 6 11 16 21];
clear V 

for icell = 1:5
    
    for iamp=1:length(id_amp)
        clear von voff tmp_on tmp_off
        hold on
        % OFF
        histogram(isi_vs{1,icell}{id_amp(iamp),1},edges,'Normalization','probability')
        % ON
        histogram(isi_vs{1,icell}{id_amp(iamp),2},edges,'Normalization','probability')
    end
    close
    
    
    h=figure('color','w','units','normalized','outerposition',[0 0 1 1])
    for iamp=1:5
        subplot(5,1,iamp)
        hold on
        if iamp ==1
            histogram(isi_vs{1,icell}{id_amp(iamp),1},edges,'Normalization','probability','FaceColor','k')
        end
        histogram(isi_vs{1,icell}{id_amp(iamp),2},edges,'Normalization','probability','FaceColor',color)
        if iamp ==1
            legend('tACS off','tACS on')
        end
        grid minor
        
        if iamp == 1
            y = ylabel('0 mV/mm')
            
        elseif iamp ==2
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        elseif iamp ==3
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        elseif iamp ==4
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        elseif iamp ==5
            y = ylabel([num2str(amplitude(id_amp(iamp)), '%.1f') ' mV/mm'])
        end
        
        set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        set(get(gca,'YLabel'),'Rotation',0)
        set(gca,'xticklabel',{[]})
        set(gca,'linewidth',1,'fontweight','bold','fontsize',15);
        %     ylim([-0.02 0.015])
        
        if iamp==5
            xticklabels({'0','50','100','150','200','250','300','350'})
            xticks([0 50 100 150 200 250 300 350])
            xlabel('ms')
        end
        
    end
    
    suptitle([' L1 NGC - cell #' num2str(icell)])
    
    print(h,sprintf('L1_NGC_cell_%d.png',icell),'-dpng','-r500');
    
end


% all cell when tacs = 0 mV/mm
figure('color','w','units','normalized','outerposition',[0 0 1 1])   

for icell=1:5   
    hold on
    histogram(isi_vs{1,icell}{1,2},edges,'Normalization','probability')    
end
grid minor
y=ylabel('0 mV/mm')
xticklabels({'0','50','100','150','200','250','300','350'})
xticks([0 50 100 150 200 250 300 350])
xlabel('ms')
legend('cell #1','cell #2','cell #3','cell #4','cell #5')

set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'YLabel'),'Rotation',0)
set(gca,'linewidth',1,'fontweight','bold','fontsize',15);
suptitle('Normalized ISI distributions for L1 NGC')

%print('ISI_L1NGC.png','-dpng','-r500');  

%%
% iamp = 44
% icell = 5
% edges = [0:5:350];

figure('color','w','units','normalized','outerposition',[0 0 1 1])
for icell=1:5
subplot(5,1,icell)
    hold on
histogram(isi_vs{1,icell}{iamp,1},edges,'Normalization','probability') % OFF
histogram(isi_vs{1,icell}{iamp,2},edges,'Normalization','probability') % ON
% grid minor
legend('OFF','ON')
% bar([12:12:348],V(icell,:))

ylim([-0.015 0.015])
set(gca,'linewidth',3,'fontweight','bold','fontsize',15)

end
suptitle([' Amplitude of - ' num2str(amplitude(iamp), '%.1f') ' V/m - L1 NGC'])



%% case 2 = we compared the same amplitude for all the cells PRE / ON / OFF

figure('color','w')
subplot(1,3,1)
for icell=1:5
    hold on
    histogram(ISI{1,icell}{iamp,1},edges,'Normalization','probability')
end
grid minor
ylabel('Number of spikes')
title('PRE')
subplot(1,3,2)
for icell=1:5
    hold on
    histogram(ISI{1,icell}{iamp,2},edges,'Normalization','probability')
end
grid minor
xlabel('time (ms)')
title('ON')
subplot(1,3,3)
for icell=1:5
    hold on
    histogram(ISI{1,icell}{iamp,3},edges,'Normalization','probability')
end
grid minor
title('POST')
legend('cell #1', 'cell #2','cell #3', 'cell #4','cell #5')
suptitle([' Amplitude of - ' num2str(amplitude(iamp), '%.1f') ' V/m'])

%% we concatenate ISI off and on
clear tmpOFF tmpON countsON countsOFF delta ED

for icell=1:5    
    for iamp=1:44
        mat_tmp{iamp,2} = ISI{1,icell}{iamp,2}; % ON
        mat_tmp{iamp,1} = [ ISI{1,icell}{iamp,1} ; ISI{1,icell}{iamp,3} ]; % OFF
    end    
    vec{1,icell} = mat_tmp;   
    for iamp=1:44
        [tmpOFF(iamp,:), ED(iamp,:)] = histcounts(vec{1,icell}{iamp,1},edges); % OFF
        [tmpON(iamp,:), ~] = histcounts(vec{1,icell}{iamp,2},edges); % ON
    end    
    countsOFF{1,icell} = tmpOFF;
    countsON{1,icell} = tmpON;
    clear mat_tmp
end

% ED(2:end,:) = [];

for icell=1:5
    delta{1,icell} = countsON{1,icell} - countsOFF{1,icell};
end

% delta is the substraction (countsON - countsOFF) which computes the
% diffences in the number of spikes in each bin for each cell for each
% amplitude.

%% here, we plot the different between the ISI when tACS is on and off

% the idea is that tACS has a high amplitude, neurons are likely to going
% to fire an action potential at peaks of tACS

clear data
icell = 1;
iamp = 33;

figure('color','w')
for iamp = 1
    x = delta{1,icell}(iamp,:);
    numpoints = length(x);
    xValues = [1:numpoints];
    
    % Find where data is positive or negative.
    positiveIndexes = x >= 0;
    negativeIndexes = x <  0;
   
    % Plot column plot (bar chart).
    bar(xValues(positiveIndexes), x(positiveIndexes),'g','BarWidth', 1)
    hold on
    bar(xValues(negativeIndexes), x(negativeIndexes),'r', 'BarWidth', 1)
    pause(0.2)
    title(['Cell #' num2str(icell) ' - ' num2str(amplitude(iamp), '%.1f') ' V/m'])
%     clf('reset')
end


%% Number of spikes

mat_spk = zeros(Namp,5);

for j=1:5
    for i=1:Namp
        tmp = SPK{1,j}(i,2);
        mat_spk(i,j) = length(tmp{1,1});    
    end 
end

% to get the firing rate evolution
mat_spk =  mat_spk - mat_spk(1,:);


figure('color','w')
plot(amplitude,mat_spk','o--')
grid minor
xlabel('EF amplitude (V/m)')
% ylabel('Number of spikes')
ylabel('Spike number rate')
xlim([0 5])

figure('color','w')
bar(amplitude,mat_spk)
grid minor
suptitle('Number of spikes')
ylabel('Number of spikes')
xlabel('Amplitudes (V/m)')
xlim([-0.05 5])

%%

tic

for j=1:5
    j
    spikes = SPK{1,j};
    
    % spikes_rates
    for i=1:Namp
        [PLV(i,1)  pre{i,1}]= getPLV(spikes{i,1},sham_tacs,tacs_start_sample,tacs_end_sample,dt,1);
        [PLV(i,2) times{i,1}] = getPLV(spikes{i,2},sham_tacs,tacs_start_sample,tacs_end_sample,dt,2);
        [PLV(i,3)  post{i,1}] = getPLV(spikes{i,3},sham_tacs,tacs_start_sample,tacs_end_sample,dt,3);
    end
    
    mat_plv{1,j} = PLV;
    mat_times{1,j} = times;
    mat_timesPre{1,j} = pre;
    mat_timesPost{1,j} = post;
    
    
%     clear PLV times
end
toc

%%

PLV=[];
for i=1:5
    PLV = [PLV  mat_plv{1,i}(:,2)];
end


figure('color','w')
plot(amplitude(1:44),PLV','o--')
xticks(amplitude)
grid minor
xlabel('EF amplitude (V/m)')
title('PLV')
xlim([0 5])
% ylim([0 0.1])
legend('cell #1', 'cell #2','cell #3', 'cell #4','cell #5')
% legend('L23PC - cell 1', 'L23PC - cell 2','L23PC - cell 3','L23PC - cell 4','L23PC - cell 5','location','southeast')
%  legend('L6PC - cell 1', 'L6PC - cell 2','L6PC - cell 3','L6PC - cell 4','L6PC - cell 5','location','southeast')
% legend('L5PC - cell 1', 'L5PC - cell 2','L5PC - cell 3','L5PC - cell 4','L5PC - cell 5','location','southeast')

subplot(1,2,2)
plot(amplitude(1:33),PLV' - PLV(1,:)'*ones(1,33),'o--')
xticks(amplitude)
grid minor
xlabel('EF amplitude (V/m)')
title('Normalized PLV')
xlim([0 5])
ylim([-1 0.9])




%%
clc

tmp_pre =[]
tmp_on=[]
tmp_post=[];

iamp = 21

for icell=1:5
   tmp_pre = [tmp_pre mat_plv{1,icell}(iamp,1)];
   tmp_on = [tmp_on mat_plv{1,icell}(iamp,2)]
   tmp_post = [tmp_post mat_plv{1,icell}(iamp,3)];
end

figure('color','w')
hold on
scatter(0.25*ones(size(tmp_pre)), tmp_pre, '+')
scatter(0.75*ones(size(tmp_on)), tmp_on, '+')
scatter(1.25*ones(size(tmp_post)),tmp_post, '+')
grid minor
% ylim([0 0.08])
xlim([0 1.5])
ylabel('PLV ')
xticks([0.25 0.75 1.25])
xticklabels({'pre', 'on','post'})
set(gca,'fontweight','bold','Fontsize',14)
% title(['PLV for amplitude of ', num2str(amplitude(iamp)) 'mV'])

















%% Different types of polar histograms

idx_amp =31;

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

   
   
%% Polar histogram: we run statistical test to see if the distribution is uniform or not

idx_cell=5

% 1) Rayleigh test
% H0: the population is distributed uniformly around the circle
% H1: the population is not distributed uniformly around the circle
clc

% small p indicates a significant departure from uniformity and indicates
% to reject the null hypothesis

alpha = mat_times{1,idx_cell}{1};
[p1 z1] = circ_rtest(alpha);

alpha = mat_times{1,idx_cell}{6};
[p2 z2]= circ_rtest(alpha);

alpha = mat_times{1,idx_cell}{11};
[p3 z3] = circ_rtest(alpha);

alpha = mat_times{1,idx_cell}{16};
[p4 z4] = circ_rtest(alpha);

alpha = mat_times{1,idx_cell}{21};
[p5 z5]= circ_rtest(alpha);

alpha = mat_times{1,idx_cell}{26};
[p6 z6]= circ_rtest(alpha);

alpha = mat_times{1,idx_cell}{31};
[p7 z7] = circ_rtest(alpha);

% we record all the p-value and z-value in vectors
pvalue(1,:) = [p1 p2 p3 p4 p5 p6]
pvalue(2,:) = [z1 z2 z3 z4 z5 z6]



%% we plot the polar hisotgrams

% change here (comment / uncomment) the color you want
% c = [0.3010, 0.7450, 0.9330]%L23
% c = [0.4660, 0.6740, 0.1880] %L4
% c = [0.4940, 0.1840, 0.556] %L5
% c = [0.8500, 0.3250, 0.0980]% L6
c='b';

% As a title for each plot, the p-value is displayed.

figure('color','w')

subplot(1,7,1)
polarhistogram(mat_times{1,idx_cell}{1},36,'FaceColor',c);
title(num2str(p1))

subplot(1,7,2)
polarhistogram(mat_times{1,idx_cell}{6},36,'FaceColor',c);
title(num2str(p2))

subplot(1,7,3)
polarhistogram(mat_times{1,idx_cell}{11},36,'FaceColor',c);
title(num2str(p3))

subplot(1,7,4)
polarhistogram(mat_times{1,idx_cell}{16},36,'FaceColor',c);
title(num2str(p4))

subplot(1,7,5)
polarhistogram(mat_times{1,idx_cell}{21},36,'FaceColor',c);
title(num2str(p5))

subplot(1,7,6)
polarhistogram(mat_times{1,idx_cell}{26},36,'FaceColor',c);
title(num2str(p6))

subplot(1,7,7)
polarhistogram(mat_times{1,idx_cell}{31},36,'FaceColor',c);
title(num2str(p7))

% print('phaseHist_L6','-depsc')



%%  Mean angle and mean length for polar histogram

% change here the index of the cell 
idx_cell = 1

% to get the mean angle
alpha = mat_times{1,idx_cell}{idx_amp,1};
alpha_bar = circ_mean(alpha) % mean alpha
R = circ_r(alpha) % mean length
S = circ_var(alpha) % the spread in a data set.

% tmp_amp = [1 6 11 16 21]
for icell = 1:5
    for i=1:44
        alpha = mat_times{1,icell}{i,1};
               
        alpha_bar(i,icell) =  circ_mean(alpha);
        R(i,icell) = circ_r(alpha);
        S(i,icell) = circ_var(alpha);
        clear alpha        
    end
end




%% 
% figure('color','w')
% plot(t/1000,somaV{1,5}(1,:), 'linewidth',1.7)
% hold on
% % plot(t/1000,tacs-80, 'linewidth',2)
% grid minor
% xlim([110 114])
% box off
% axis off
% print('baseline','-depsc')
% 
% figure('color','w')
% % figure('color','w')
% plot(t/1000,somaV{1,5}(11,:), 'linewidth',1.7)
% hold on
% % plot(t/1000,tacs*40-140,'linewidth',4)
% grid minor
% xlim([133 137])
% box off
% axis off
% % print('strongTACS','-depsc')
% 
% figure('color','w')
% plot(t/1000,somaV{1,5}(6,:), 'linewidth',1.7)
% hold on
% % plot(t/1000,tacs-80, 'linewidth',2)
% grid minor
% xlim([400 404])
% box off
% axis off
% print('weakTACS','-depsc')