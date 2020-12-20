clear all
close all
clc

% Description:
%   The code is created to investigate the subthreshold neuron activity for
%   different tACS amplitudes. 
% Code created on 7/30/2020

% Correspondence: htran@umn.edu

% % -------------------------------------------------------------------------

clear all
close all
clc

% change here you path of the working directory
cd('D:\subthreshold_investigation')
amplitude = [0:0.1:2 2.2:0.2:5 6.25:1.25:10 12.5:2.5:25];
Namp = length(amplitude);
Ncell = 25;

% we define the 5 different neurons types
L1GNC = [1:5];
L23PC = [6:10];
L4LBC = [11:15];
L5PC = [16:20];
L6PC = [21:25];

% we load the data
tic
for icell = 1:Ncell
    % change here the path of the working directory
    cd(['D:\subthreshold_investigation\results_cell' num2str(icell, '%.1d')])
    disp('==================')
    disp(['Cell loaded: ' num2str(icell)])
    for j=1:Namp
        cd(['tacs_' num2str(j, '%.1d')])       
        if j==1
            t=load('t.txt');
            dt = t(2) - t(1);
            tacs = load('TACS.txt');
            tacs = tacs./max(tacs); % we normalize the tacs
        end       
        tmp_soma(j,:)= load('somaV.txt');
        cd('..\')
    end
    somaV{1,icell} = tmp_soma;
    cd('..\')
    clear tmp_soma
end
toc
%%

clear val tmpBase valuesPeaks baseline

% we computed the baseline and the max for each cell for each amplitude
val=[];
tic
for icell =1:25
    tmp = somaV{1,icell};      
    for i=1:length(amplitude)    
       tmpBase(i,1) = mean(tmp(i,1900/dt:2000/dt))    ;
       [valtmp,~] = findpeaks(tmp(i,1400/dt:1500/dt), 'MinPeakHeight',tmpBase(i,1));
       val(i,1)=mean(valtmp);
    end      
    baseline(:,icell) = tmpBase; % mV
    valuesPeaks(:,icell) = val; % mV
end
toc 

clear tmp val tmpBase valtmp

valuesPeaks(1,:) = baseline(1,:);

%% We compute the difference in the membrane voltage from the baseline

% delta is a matrix (N_amplitudes X N_cells)
delta= valuesPeaks - baseline;
delta(1,6:10) = G(2,6:10) - 0.012; % G is a matrix to get rid of NaN
delta(1,11:15) = G(2,11:15) - 0.005;

% conversion mV in V
delta = delta*0.001;

figure('color','w')
for i=6:15
    hold on
plot(amplitude(1:26), delta(1:26,i),'o--')
end
grid minor
xticks(amplitude)
xlabel('EF amplitude (V/m)')
ylabel('Polarization (V)')

%% Linear regression on the polarization matrix

for icell=1:Ncell
    [P(icell,:), S(icell,1)] = polyfit(amplitude(1:36)',delta(1:36,icell),1);
end
% P(:,1) correspond to the slope (polarization length in meters)

%% we plot the polarization length for each cell type
for icell=1:Ncell
    y = delta(1:36,icell);
    R2(icell,1) = 1 - (S(icell,1).normr/norm(y - mean(y)))^2;
end

% we compute the coefficient of determination
R2 = reshape(R2,[5,5]);
PP = reshape(P(:,1),[5,5]); % slope in meters
PP = PP*1000; % we convert meter into millimeters

figure('color','w')
boxplot(PP,'Labels',{'L1GNC','L23PC','L4LBC','L5PC','L6PC'})
grid minor
ylabel('Polarization Length \lambda (mm)')

%% we group cell according to their type
%  Reference: [radman2009role]

interneuron = [PP(:,1) ; PP(:,3)];
L23pyr = PP(:,2);
L56pyr = [PP(:,4) ; PP(:,5)];

mat = [interneuron ; L23pyr ; L56pyr];
grp =[zeros(1,10), ones(1,5), 2*ones(1,10)];

figure('color','w')
boxplot(mat, grp,'PlotStyle','compact','LabelOrientation','horizontal','Labels',{'Interneurons','L23 PC','L56 PC'})
ylabel('Polarization Length \lambda (mm)')
title('Cortical cell type polarization sensitivty')

%% We plot the somatic membrane polarization for a cell for each type

clc
close all
figure('color','w','Position', [10 10 700 700])
hold on
%L1 NGC
scatter(amplitude(1:2:35), delta(1:2:35,2)*1000,75,'o','MarkerEdgeColor','b','linewidth',4)
%L2/3 PC
scatter(amplitude(1:2:35), (delta(1:2:35,8)*1000+0.01) - 0.03598,75 ,'o','MarkerEdgeColor',[0.3010, 0.7450, 0.9330],'linewidth',4)
%L4 PC
scatter(amplitude(1:2:35), (delta(1:2:35,14)*1000)- 0.007373,75,'o','MarkerEdgeColor',[0.4660, 0.6740, 0.1880],'linewidth',4)
%L5 PC
scatter(amplitude(1:2:35), delta(1:2:35,16)*1000,75,'o','MarkerEdgeColor',[0.4940, 0.1840, 0.556],'linewidth',4)
%L6 PC
scatter(amplitude(1:2:35), delta(1:2:35,21)*1000,75,'o','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'linewidth',4)
% xticks(amplitude(1:2:35))
xlim([0 3])
set(gca,'linewidth',5,'fontweight','bold','fontsize',20);
yticks([0.1 0.2 0.3 0.4])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

%% we plot the polarization length for all the cells

size =4;
figure('color','w','Position', [10 10 700 700])
hold on
%L1 NGC
scatter(1.2,PP(1,1),300,'b','+','linewidth',size)
scatter(0.9,PP(2,1),300,'b','o','linewidth',size)
scatter(1,PP(3,1),300,'b','x','linewidth',size)
scatter(1.1,PP(4,1),300,'b','s','linewidth',size)
scatter(0.8,PP(5,1),300,'b','d','linewidth',size)
%L2/3 PC
scatter(2.1,PP(1,2),300,[0.3010, 0.7450, 0.9330],'+','linewidth',size)
scatter(2.2,PP(2,2),300,[0.3010, 0.7450, 0.9330],'o','linewidth',size)
scatter(2,PP(3,2),300,[0.3010, 0.7450, 0.9330],'x','linewidth',size)
scatter(1.9,PP(4,2),300,[0.3010, 0.7450, 0.9330],'s','linewidth',size)
scatter(1.8,PP(5,2),300,[0.3010, 0.7450, 0.9330],'d','linewidth',size)
%L4 LBC
scatter(3.1,PP(1,3),300,[0.4660, 0.6740, 0.1880],'+','linewidth',size)
scatter(3.2,PP(2,3),300,[0.4660, 0.6740, 0.1880],'o','linewidth',size)
scatter(3,PP(3,3),300,[0.4660, 0.6740, 0.1880],'x','linewidth',size)
scatter(2.9,PP(4,3),300,[0.4660, 0.6740, 0.1880],'s','linewidth',size)
scatter(2.8,PP(5,3),300,[0.4660, 0.6740, 0.1880],'d','linewidth',size)
%L5 PC
scatter(4.1,PP(1,4),300,[0.4940, 0.1840, 0.556],'+','linewidth',size)
scatter(4.2,PP(2,4),300,[0.4940, 0.1840, 0.5560],'o','linewidth',size)
scatter(3.9,PP(3,4),300,[0.4940, 0.1840, 0.5560],'x','linewidth',size)
scatter(3.8,PP(4,4),300,[0.4940, 0.1840, 0.5560],'s','linewidth',size)
scatter(4,PP(5,4),300,[0.4940, 0.1840, 0.5560],'d','linewidth',size)
%L6 PC
scatter(5,PP(1,5),300,[0.8300, 0.3250, 0.0980],'+','linewidth',size)
scatter(4.9,PP(2,5),300,[0.8300, 0.3250, 0.0980],'o','linewidth',size)
scatter(4.8,PP(3,5),300,[0.8300, 0.3250, 0.0980],'x','linewidth',size)
scatter(5.2,PP(4,5),300,[0.8300, 0.3250, 0.0980],'s','linewidth',size)
scatter(5.1,PP(5,5),300,[0.8300, 0.3250, 0.0980],'d','linewidth',size)

xticks(1:5)
xticklabels({'L1 NGC','L2/3 PC','L4 LBC','L5 PC','L6 PC'})
xlim([0 6])
y=ylabel('Polarization length (mm)')
set(gca,'linewidth',5,'fontweight','bold','fontsize',20);


%% Statistical tests  / t-test

%  interneurons VS L23
[h,p,ci,stats] = ttest2([PP(:,1) ; PP(:,3)], PP(:,2));

% interneurons VS L5
[h,p,ci,stats] = ttest2([PP(:,1) ; PP(:,3)], PP(:,4));

% interneurons VS L6
[h,p,ci,stats] = ttest2([PP(:,1) ; PP(:,3)], PP(:,5));
