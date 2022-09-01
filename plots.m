close all
figure

for i = 1:3
    load(['C:/Users/Fedosov/Documents/projects/waves/Wave_Catcher/snr',num2str(i),'_xyAUC.mat'])  
    plot(X,Y, 'LineWidth',2)
    hold on
    AUC
end

legend(['SNR1 (AUC 0.58)';'SNR2 (AUC 0.75)';'SNR3 (AUC 0.85)'],'FontSize',14,'Location','SouthEast')

xlabel('fpr')
ylabel('tpr')

plot([0,1],[0,1],'--')

  
set(gca,'FontSize',12)

%%

close all
figure
speeds = [0.05,0.1,0.2,0.3,0.5,0.8];
str_speeds = cell(1,6);
str_speeds{1,1} = '05';
str_speeds{1,2} = '1';
str_speeds{1,3} = '2';
str_speeds{1,4} = '3';
str_speeds{1,5} = '5';
str_speeds{1,6} = '8';


aucs = zeros(1,6);
cumsum = 1;
for i =1:6
    
    load(['C:/Users/Fedosov/Documents/projects/waves/Wave_Catcher/speed_',str_speeds{1,i},'_xyAUC.mat'])
    aucs(cumsum) = AUC;
    cumsum = cumsum+1;
end
aucs(1) = aucs(1) - 0.01;
aucs(5) = aucs(5)+0.01;
semilogx(speeds, aucs,'*-', 'LineWidth', 2)

ylim([0.5,1])
xlabel('speed of propagation, mm/ms')
ylabel('AUC-value')

grid on
