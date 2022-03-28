
positionlist = [0.035, 0.3, 0.28, 0.60;... %75R-->P, outcome o
    0.035, 0.05, 0.28, 0.2; ...%75R-->P, outcome x
    0.365, 0.3, 0.28, 0.60;...
    0.365, 0.05, 0.28, 0.2;...
    0.695, 0.05, 0.28, 0.85];
triallist= {'trial_80Rwo_ind','trial_80Pno_ind','trial_20R_ind'};
titlelist = {'75R', '75P', '0', '75R ¡æ 75P','75p ¡æ 75R','0'};
mkdir([figure_folder, '\C_raw']);
oo=64;
C_event = reshape(C_raw(oo,set_eventframe(:)),[],12*fr);
coloraxis = [prctile(C_event(:),0.1), prctile(C_event(:),99.95)];
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 40 15]);
colormap(othercolor(27))
subplot(1,2,1)
a=[];
for ii = 1:3;
    eval(['trialind=',triallist{ii},';',]);
    eval(['lickall=',licklist{ii},';',]);
    trial_first = sum(trialind(30:sum(rev_b)));
    hold on
    a =[a;C_event(find(trialind(1:sum(rev_b)),10),:)];
    
end
imagesc(a);caxis(coloraxis);
line([0.5*fr 0.5*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %0.5 baseline
line([2*fr 2*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %2 cue
line([3.5*fr 3.5*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %3.5 delay1
line([4*fr 4*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %4 motor
line([5.5*fr 5.5*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %5.5 delay2
line([8*fr 8*fr], [0 30],'color','w', 'linestyle', '-.','linewidth',1.5); %8 motor off
line([1 12*fr], [trial_first trial_first],'color','w', 'linewidth',2);
xlim([1 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})
ylim([1 sum(trialind)]); yticks(sum(trialind));

subplot(1,2,2)
a=[];
for ii = 1:3;
    eval(['trialind=',triallist{ii},';',]);
    eval(['lickall=',licklist{ii},';',]);
    hold on
    a =[a;C_event(find(trialind(sum(rev_b)):end,10),:)];
    
end
imagesc(a);caxis(coloraxis);
line([0.5*fr 0.5*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %0.5 baseline
line([2*fr 2*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %2 cue
line([3.5*fr 3.5*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %3.5 delay1
line([4*fr 4*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %4 motor
line([5.5*fr 5.5*fr], [0 30],'color','w', 'linestyle', ':','linewidth',1.5); %5.5 delay2
line([8*fr 8*fr], [0 30],'color','w', 'linestyle', '-.','linewidth',1.5); %8 motor off
line([1 12*fr], [trial_first trial_first],'color','w', 'linewidth',2);
xlim([1 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})
ylim([1 sum(trialind)]); yticks(sum(trialind));

saveas(f1,[figure_folder,'\C_raw\','PSTH1_',num2str(oo),'_craw_ex.tif']);
close all
