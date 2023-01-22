%% setting
behavior_mat = 'dHP08_rwprobrev02_20210203_155560.mat';%need correction
cell_mat = '04-Jan_21_18_03.mat';
fr=30; % frame rate


%% initiation 1
event2mat_ca2('event_signal.csv') %make event mat file
cellsort_ca %for make final cell group for Craw
data_set_rec(behavior_mat)
plot_reversal_rec(behavior_mat)
%% initiation 2
figure_folder = [cd,'\figure']; %save figure folder
load('event.mat') % for event file
cmap1 = [30 45 232; 120 177 255; 232 126 58; 255 30 70]/255; %blue skyblue pink red
%cell cut with event frequency
S = full(S);
rec_totaltime = size(S,2) ./fr;
Sfrequency = sum(S>0,2)./rec_totaltime;
neuron_over15= Sfrequency>0.015;
neuron_over25= Sfrequency>0.025;
neuron_over50 = Sfrequency>0.05;
C_raw0 = C_raw; %save neuron
S_raw0 = S;

C_raw=C_raw0(neuron_over15,:);
S_raw=S_raw0(neuron_over15,:);
clear neuron
save('initiation.mat')
%% start
% load('initiation.mat')
day = datetime('today');
d = datestr(day,'yymmdd');
load(behavior_mat) %load behavior file
frame_time = sync(sync(:,2)==1,1);
%% align with event1 - calculate automatic remove framed
% preprocess(:,1)=all frame (:,2)=cut by automatic (:,3) = gpio
preprocess_cut = readtable('preprocess.csv','ReadRowNames',false);
preprocess_cut = table2array(preprocess_cut);

for ii= 1:size(preprocess_cut,1);
    preprocess_cut(ii,3) = sum(sync(:,2)==1);
end

cut_1st = preprocess_cut(:,1)-preprocess_cut(:,2); %cut automatic revome
cut_end = preprocess_cut(:,1)-preprocess_cut(:,3); %cut dropped sync

for ii = 1:size(preprocess_cut,1);
    if ii ==1; cut = 1:cut_1st(ii);
    else cut = [cut, sum(preprocess_cut(1:ii-1,3))+1:sum(preprocess_cut(1:ii-1,3))+cut_1st(ii)];
    end
end


frame_time = sync(sync(:,2)==1,1);
rec_start_time = frame_time(1,1);
frame_time(cut) = [];

% Cut C_raw & S for skipped gpio
C_raw = C_raw0; S_raw = S_raw0;
for ii = 1:size(preprocess_cut,1);
    if ii ==1; cut = preprocess_cut(1,2)-(1:cut_end(ii))+1;
    else cut = [cut, sum(preprocess_cut(1:ii,2)):sum(preprocess_cut(1:ii,2))-cut_end(ii)+1];
    end
end
C_raw(:,cut) =[];S_raw(:,cut)=[];
C_raw00 = C_raw;
C_raw=C_raw(neuron_over15,:);
S_raw=S_raw0(neuron_over15,:);

nframe = size(C_raw,2);

%% align with event2 - sync & gpio & trig
% Trig
trig_on = trig(find(diff(trig(:,2))==1)+1,1);
trig_off = trig(find(diff(trig(:,2))==-1)+1,1);
if isempty(trig_on); trig_on = trig(1,1); trig_off = trig(2,1); end; % for all rec.

% GPIO recorded all session
GPIO1_in = GPIO1(:,2)>40000;
diff_GPIO1 = diff(GPIO1_in);

% Gpio on & off time
event_on_all = GPIO1(find(diff_GPIO1==1)+1,1);
event_off_all = GPIO1(find(diff_GPIO1==-1)+1,1);

% remove arduino off noise
if length(event_on_all)>nTrial
    delete_last=abs(event_off_all-event_on_all-0.5)>0.2;
    event_on_all(delete_last)=[];
    event_off_all(delete_last)=[];
end

rec_region = [1, nTrial];

%find GPIO1 input frame
GPIO1onframe = NaN(length(event_on_all),1); GPIO1offframe = NaN(length(event_off_all),1);
for jj = 1:length(event_on_all)
    Lind = find(frame_time<event_on_all(jj,1));
    Hind = find(event_on_all(jj,1)<=frame_time);
    GPIO1onframe(jj,1) = intersect(Lind,Hind-1);
end
for kk = 1:length(event_off_all);
    Lind = find(frame_time<event_off_all(kk,1));
    Hind = find(event_off_all(kk,1)<=frame_time);
    GPIO1offframe(kk,1) = intersect(Lind,Hind-1);
end


%% align with event1 - behavior & sync for example
trialtime_beh = diff(stateTime(:,1))./1000000; %total time for each trial
rec_start_trial = []; rec_end_trial = [];
trialtime_GPIO1set = diff(event_on_all);
for ll = 1:length(trig_on) %find recording start trial
    trialtime_GPIO1 = trialtime_GPIO1set(sum(rec_region(1:ll)):sum(rec_region(1:ll+1))-2);
    mm_tmp = 10;
    for mm = 1:nTrial-mm_tmp
        timebeh = trialtime_beh(mm:mm+mm_tmp-1,1);
        timediff1 = timebeh-trialtime_GPIO1(1:mm_tmp,1);
        timediff2 = timebeh-trialtime_GPIO1(end-mm_tmp+1:end,1);
        if abs(timediff1)<0.1
            rec_start_trial= [rec_start_trial;mm];
        end
        if abs(timediff2)<0.1
            rec_end_trial = [rec_end_trial; mm+mm_tmp];
        end
    end
end

% find anova diff
% set_rec_tmp=[];
% for jj = 1:200;
%     for ii = 1:200;
%         [p,tbl,stats] = anovan(lickNum(ii:end-jj+1), odorCue(ii:end-jj+1),'display','off');
%         a =multcompare(stats,'display','off');
%         if (sum(outcomeContingency(:,1))==3 & a(:,6)<0.05) 
%             set_rec_tmp = [set_rec_tmp;ii,jj];
%         elseif (sum(outcomeContingency(:,1))==4);
%             ai = find(trial_kind==1);
%             set_rec_tmp = [1 1];
%         end
%     end
% end
% [~,i_m] = min(sum(set_rec_tmp,2));
% set_rectrial = [set_rec_tmp(i_m,1):length(lickNum)+1-set_rec_tmp(i_m,2)];

set_rectrial = [1:nTrial];

trial_gpio = 1:length(GPIO1onframe);

[X,Y] = meshgrid(0:12*fr-1,GPIO1onframe); % 0~12
set_eventframe = X+Y; %meshgrid for event frame
used_cue = unique(odorCue);
reward_cue = used_cue(outcomeContingency(2,used_cue+1)==1);
% r25=1 75rp=2 75pr=3 ze=4 for 4cue
% r80=1 r50=2 r20=3 p80=4 for 3cue
trial_kind= lickEachCue.cue_kind;
outcome_rec = waterReward(set_rectrial);
odorCue_rec = odorCue(set_rectrial);
if length(rev_index)>4;
    rev_b= [ones(1,rev_index(1,5)),zeros(1,nTrial-rev_index(1,5))];
else; rev_b = ones(1,nTrial);end
rev_b = rev_b(set_rectrial);

%% Draw all craw
imagesc(C_raw)
set(gca,'CLim',[-5 15])

line([GPIO1onframe(odorCue_rec==0) GPIO1onframe(odorCue_rec==0)],[1 size(C_raw,1)],'color','b','linewidth',2,'linestyle',':');
line([GPIO1onframe(odorCue_rec==1) GPIO1onframe(odorCue_rec==1)],[1 size(C_raw,1)],'color','r','linewidth',2,'linestyle',':');
%line([GPIO1onframe(odorCue_rec==2) GPIO1onframe(odorCue_rec==2)],[1 size(C_raw,1)],'color','c','linewidth',2,'linestyle',':');
line([GPIO1onframe(odorCue_rec==3) GPIO1onframe(odorCue_rec==3)],[1 size(C_raw,1)],'color','k','linewidth',2,'linestyle',':');

%% First lick
if length(rev_index)>4; %rev on
    lick_trial = [timeline{1,1};timeline{1,2}];
else lick_trial = timeline{1,1}; end;
%firstlick = [1;find(diff(lick_trial(:,4))~=0)+1];
% r80=1 r50=2 r20=3 p80=4 for 3cue
firstlick_set = lick_trial;
trial_80R_ind = (find(trial_kind==1)-1)==odorCue;
trial_80Rwo_ind = trial_80R_ind & waterReward;
trial_80Rwx_ind = trial_80R_ind & ~waterReward;
trial_50R_ind = (find(trial_kind==2)-1)==odorCue;
trial_50Rwo_ind = trial_50R_ind & waterReward;
trial_50Rwx_ind = trial_50R_ind & ~waterReward;
trial_20R_ind = (find(trial_kind==3)-1)==odorCue; %% also 0 cue
trial_20Rwo_ind = trial_20R_ind & waterReward;
trial_20Rwx_ind = trial_20R_ind & ~waterReward;
trial_80P_ind = (find(trial_kind==4)-1)==odorCue;
trial_80Pno_ind = trial_80P_ind & waterReward;
trial_80Pnx_ind = trial_80P_ind & ~waterReward;




firstlick_set(:,1) = firstlick_set(:,1)*fr/1000000;

firstlick_hab = firstlick_set(ismember(firstlick_set(:,2),find(trial_80R_ind)),[1,2,4]);
firstlick_80ro = firstlick_set(ismember(firstlick_set(:,2),find(trial_80Rwo_ind)),[1,2,4]);
firstlick_80rx = firstlick_set(ismember(firstlick_set(:,2),find(trial_80Rwx_ind)),[1,2,4]);
firstlick_80po = firstlick_set(ismember(firstlick_set(:,2),find(trial_80Pno_ind)),[1,2,4]);
firstlick_80px = firstlick_set(ismember(firstlick_set(:,2),find(trial_80Pnx_ind)),[1,2,4]);
firstlick_20ro = firstlick_set(ismember(firstlick_set(:,2),find(trial_20Rwo_ind)),[1,2,4]);
firstlick_20rx = firstlick_set(ismember(firstlick_set(:,2),find(trial_20Rwx_ind)),[1,2,4]);
firstlick_50ro = firstlick_set(ismember(firstlick_set(:,2),find(trial_50Rwo_ind)),[1,2,4]);
firstlick_50rx = firstlick_set(ismember(firstlick_set(:,2),find(trial_50Rwx_ind)),[1,2,4]);

licksetlist = {'firstlick_hab','firstlick_80ro','firstlick_80rx','firstlick_80po','firstlick_80px','firstlick_50ro','firstlick_50rx','firstlick_20ro','firstlick_20rx'};
for ii = 1:9
    eval(['seta=',licksetlist{ii},';',]);
    seta = seta(ismember(seta(:,2),set_rectrial),:);
    
    lickarrange = NaN(size(seta,1),1);
    diffa = diff(seta);
    for jj = 1:size(seta,1);
        if jj==1;
            lickarrange(jj)=1;
            tmp=1;
        elseif diffa(jj-1,2)~=0
            tmp = tmp+1;
            lickarrange(jj) = tmp;
        else
            lickarrange(jj) = tmp;
        end
    end
    seta(:,4) = lickarrange;
    eval([licksetlist{ii},'=seta',';',]);
    
    
end

save('behavior_align.mat')

%% Draw figure - psth by C - habituation
if stateTime(:,4)==0;
    linelist = diff(stateTime(1,[1 2 3 5]), 1,2)./1000000;
else; linelist = diff(stateTime(1,:),1,2)./1000000;
end
linelist = cumsum(linelist);
mkdir([figure_folder, '\C_raw']);

for oo = 1:size(C_raw,1)
    C_event = reshape(C_raw(oo,set_eventframe(:)),[],12*fr);
    coloraxis = [prctile(C_event(:),0.05), prctile(C_event(:),99.98)];
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 15]);
    colormap(othercolor(28))
    hold on
    imagesc(C_event); caxis(coloraxis);
    scatter(firstlick_hab(:,1),firstlick_hab(:,4),2.5,'k','filled','linewidth',1.5)
    for lines = 1:length(linelist-1);
        line([linelist(lines)*fr linelist(lines)*fr], [0 sum(trial_80R_ind)],'color','w', 'linestyle', ':','linewidth',1.8); %0.5 baseline
    end
    line([linelist(end)*fr linelist(end)*fr], [0 sum(trial_80R_ind)],'color','c', 'linestyle', ':','linewidth',1.8); %0.5 baseline
    line([(2.5+linelist(end))*fr (2.5+linelist(end))*fr], [0 sum(trial_80R_ind)],'color','w', 'linestyle', '-.','linewidth',2); %8 motor off
    if length(linelist)==3;
        xlim([1 9*fr]); xticks([(linelist+[diff(linelist),2]./2)*fr]); xticklabels({'Cue', 'Delay','Outcome'})
    else
        xlim([1 12*fr]); xticks([(linelist+[diff(linelist),2]./2)*fr]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})
        line([3.5*fr 3.5*fr], [0 sum(trial_80R_ind)],'color','w', 'linestyle', ':','linewidth',1.8)
    end
    ylim([1 sum(trial_80R_ind)]); yticks(sum(trial_80R_ind))
    %title('75R --> 75P')
    title('Habituation')
    saveas(f1,[figure_folder,'\C_raw\','Hab',num2str(oo),'_craw.tif']);
    close all
end
% mean
mkdir([figure_folder, '\C_raw_mean']);
for oo = 1:size(C_raw,1)
    C_event = reshape(C_raw(oo,set_eventframe(:)),[],12*fr);
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
    hold on
    Ctemp = mean(C_event(:,1:0.5*fr),2);
    CM = mean(Ctemp); CSTD=std(Ctemp);
    plot(mean(movmean((C_event(trial_80Rwo_ind,:)-CM)./CSTD,[4,4],2)),'linewidth',1.2,'color','b')
    plot(mean(movmean((C_event(trial_80Rwx_ind,:)-CM)./CSTD,[4,4],2)),':','linewidth',1.2,'color', 'b')
    lim = axis;
    for lines = 1:length(linelist-1);
        line([linelist(lines)*fr linelist(lines)*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.8); %0.5 baseline
    end
    line([(2.5+linelist(end))*fr (2.5+linelist(end))*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.','linewidth',2); %8 motor off
    if length(linelist)==3;
        xlim([1 9*fr]); xticks([(linelist+[diff(linelist),2]./2)*fr]); xticklabels({'Cue', 'Delay','Outcome'})
    else
        xlim([1 12*fr]); xticks([(linelist+[diff(linelist),2]./2)*fr]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})
        line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.8)
    end
    hold off
    saveas(f1,[figure_folder, '\C_raw_mean\','PSTHmean_',num2str(oo),'_craw.tif']);
    close all
end

%% Draw figure - psth by C - Rwpn,Rwprob

if stateTime(:,4)==0;
    linelist = diff(stateTime(1,[1 2 3 5]), 1,2)./1000000;
else; linelist = diff(stateTime(1,:),1,2)./1000000;
end
linelist = cumsum(linelist);
positionlist_rp = [0.035, 0.3, 0.28, 0.60;... %75R-->P, outcome o
    0.035, 0.05, 0.28, 0.2; ...%75R-->P, outcome x
    0.365, 0.3, 0.28, 0.60;...
    0.365, 0.05, 0.28, 0.2;...
    0.695, 0.05, 0.28, 0.85];
triallist_rp= {'trial_80Rwo_ind','trial_80Rwx_ind','trial_80Pno_ind','trial_80Pnx_ind','trial_20R_ind'};
licklist_rp = {'firstlick_80ro','firstlick_80rx','firstlick_80po','firstlick_80px','firstlick_20rx'};
titlelist_rp = {'75R', '75P', '0', '75R ¡æ 75P','75p ¡æ 75R','0'};
positionlist_prob = [0.035, 0.3, 0.28, 0.60;... %75R-->P, outcome o
    0.035, 0.05, 0.28, 0.2; ...%75R-->P, outcome x
    0.365, 0.7, 0.28, 0.20;...
    0.365, 0.05, 0.28, 0.6;...
    0.695, 0.05, 0.28, 0.85];
positionlist_prob_rev = [0.035, 0.5, 0.28, 0.4;... %75R-->P, outcome o
    0.035, 0.05, 0.28, 0.4; ...%75R-->P, outcome x
    0.365, 0.5, 0.28, 0.4;... 
    0.365, 0.05, 0.28, 0.4;...
    0.695, 0.05, 0.28, 0.85];
triallist_prob= {'trial_80Rwo_ind','trial_80Rwx_ind','trial_50Rwo_ind','trial_50Rwx_ind','trial_20R_ind'};
licklist_prob = {'firstlick_80ro','firstlick_80rx','firstlick_50ro','firstlick_50rx','firstlick_20rx'};
titlelist_prob = {'75R', '25R', '0', '75R ¡æ 25R','25R ¡æ 75R','0'};
if ismember(trial_kind,[0 1 3 4]);
    positionlist = positionlist_rp; triallist = triallist_rp; licklist = licklist_rp; titlelist = titlelist_rp;
elseif ismember(trial_kind, [0 1 2 3]);
    positionlist = positionlist_prob;triallist = triallist_prob; licklist = licklist_prob; titlelist = titlelist_prob;
    if sum(rev_b)<420; positionlist = positionlist_prob_rev; end
end
mkdir([figure_folder, '\C_raw']);

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 40 15]);
for oo = 1:size(C_raw,1)
    C_event = reshape(C_raw(oo,set_eventframe(:)),[],12*fr);
    coloraxis = [prctile(C_event(:),0.1), prctile(C_event(:),99.95)];
    
    colormap(othercolor(27))
    for ii = 1:5;
        subplot('Position', positionlist(ii,:)) %x y w h
        eval(['trialind=',triallist{ii},';',]);
        eval(['lickall=',licklist{ii},';',]);
        trial_first = sum(trialind(1:sum(rev_b)));
        hold on
        imagesc(C_event(trialind,:)); caxis(coloraxis);
        scatter(lickall(:,1),lickall(:,4),5,'k','filled','linewidth',1.5)
        line([0.5*fr 0.5*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %0.5 baseline
        line([2*fr 2*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %2 cue
        line([3.5*fr 3.5*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %3.5 delay1
        line([4*fr 4*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %4 motor
        line([5.5*fr 5.5*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %5.5 delay2
        line([8*fr 8*fr], [0 sum(trialind)],'color','w', 'linestyle', '-.','linewidth',1.5); %8 motor off
        line([1 12*fr], [trial_first trial_first],'color','w', 'linewidth',2);
        xlim([1 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})
        ylim([1 sum(trialind)]); yticks(sum(trialind));
        if rem(ii,2)==1;
            title(titlelist{fix(ii/2)+1+3*(sum(rev_b)<nTrial)});
        end
        
    end
    saveas(f1,[figure_folder,'\C_raw\','PSTH1_',num2str(oo),'_craw_ex.tif']);
    clf
end

%% Draw figure - psth by S - Rwpn
if stateTime(:,4)==0;
    linelist = diff(stateTime(1,[1 2 3 5]), 1,2)./1000000;
else; linelist = diff(stateTime(1,:),1,2)./1000000;
end
linelist = cumsum(linelist);
positionlist = [0.035, 0.3, 0.28, 0.60;... %75R-->P, outcome o
    0.035, 0.05, 0.28, 0.2; ...%75R-->P, outcome x
    0.365, 0.3, 0.28, 0.60;...
    0.365, 0.05, 0.28, 0.2;...
    0.695, 0.05, 0.28, 0.85];
triallist= {'trial_80Rwo_ind','trial_80Rwx_ind','trial_80Pno_ind','trial_80Pnx_ind','trial_20R_ind'};
licklist = {'firstlick_80ro','firstlick_80rx','firstlick_80po','firstlick_80px','firstlick_20rx'};
titlelist = {'75R', '75P', '0', '75R ¡æ 75P','75p ¡æ 75R','0'};
mkdir([figure_folder, '\S_raw']);
for oo = 1:size(S_raw,1)
    S_event = reshape(S_raw(oo,set_eventframe(:)),[],12*fr);
    coloraxis = [prctile(S_event(:),0), prctile(S_event(:),100)];
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 40 15]);
    
    colormap(hot)
    for ii = 1:5;
        subplot('Position', positionlist(ii,:)) %x y w h
        eval(['trialind=',triallist{ii},';',]);
        eval(['lickall=',licklist{ii},';',]);
        trial_first = sum(trialind(1:sum(rev_b)));
        hold on
        imagesc(S_event(trialind,:)); caxis(coloraxis);  
        scatter(lickall(:,1),lickall(:,4),1.5,'w','filled','linewidth',1.5)
        line([0.5*fr 0.5*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %0.5 baseline
        line([2*fr 2*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %2 cue
        line([3.5*fr 3.5*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %3.5 delay1
        line([4*fr 4*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %4 motor
        line([5.5*fr 5.5*fr], [0 sum(trialind)],'color','w', 'linestyle', ':','linewidth',1.5); %5.5 delay2
        line([8*fr 8*fr], [0 sum(trialind)],'color','w', 'linestyle', '-.','linewidth',1.5); %8 motor off
        line([1 12*fr], [trial_first trial_first],'color','w', 'linewidth',2);
        xlim([1 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})
        ylim([1 sum(trialind)]); yticks(sum(trialind));
        if rem(ii,2)==1;
            title(titlelist{fix(ii/2)+1+3*(sum(rev_b)<nTrial)});
        end
    end
    saveas(f1,[figure_folder,'\S_raw\','PSTH1_',num2str(oo),'.tif']);
    close all
end

%% z-score cell activity
% clear
% load('behavior_align.mat','C_raw','set_eventframe','fr',...
%     'waterReward','trial_80R_ind','trial_50R_ind','trial_20R_ind','trial_80P_ind',...
%     'trial_80Rwo_ind','trial_80Rwx_ind','trial_50Rwo_ind','trial_50Rwx_ind','trial_20Rwo_ind','trial_20Rwx_ind','trial_80Pno_ind','trial_80Pnx_ind')

C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_b = mean(C_event(:,0*fr+1:0.5*fr,:),2); C_b = permute(C_b,[3 1 2]);
% C_z_base = (C_raw-mean(C_b,2))./std(C_b,0,2);
% C_z_base = C_z_base (:,set_eventframe(:));
% C_z_base = reshape(C_z_base', [],12*fr,size(C_z_base,1));
C_z_all = (C_raw-mean(C_raw(:,set_eventframe(1,1):set_eventframe(end,end)),2))./std(C_raw(:,set_eventframe(1,1):set_eventframe(end,end)),0,2);
C_z_all = C_z_all(:,set_eventframe(:));
C_z_all = reshape(C_z_all', [],12*fr,size(C_z_all,1));
c_tmp = C_z_all;
save('behavior_align.mat','C_z_all','-append')

C_z_cue=cat(3,permute(mean(c_tmp(trial_80R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_50R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_20R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80P_ind,0*fr+1:5.5*fr,:),1),[3 2 1]));
C_z_out=cat(3,permute(mean(c_tmp(trial_80Rwo_ind,5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80Rwx_ind,5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_50Rwo_ind,5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_50Rwx_ind,5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_20R_ind ,5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80Pno_ind,5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80Pnx_ind,5*fr+1:11*fr,:),1),[3 2 1]));

st=strsplit(cd,'\');
eval(['C_z_cue_',st{4},'=C_z_cue; C_z_out_',st{4},'=C_z_out;'])
eval(['save C_z_',st{4},'.mat C_z_cue_',st{4},' C_z_out_',st{4}])

zzset = [mean(C_z_all(waterReward==1,:,:),1);mean(C_z_all(waterReward==0,:,:),1);...
    mean(C_z_all(trial_80Rwo_ind,:,:),1);mean(C_z_all(trial_50Rwo_ind,:,:),1);mean(C_z_all(trial_20Rwo_ind,:,:),1);...
    mean(C_z_all(trial_80Rwx_ind,:,:),1);mean(C_z_all(trial_50Rwx_ind,:,:),1);mean(C_z_all(trial_20Rwx_ind,:,:),1)];
zzset = permute(zzset,[2 3 1]);
save('RPE_forz.mat','zzset')
% Normalized 0-1 cell activity
C_event = C_raw(:,set_eventframe(:));
C_nor_all = rescale_mr(C_event');
C_nor_all = reshape(C_nor_all, [],12*fr,size(C_nor_all,2));
c_tmp = C_nor_all;
save('behavior_align.mat','C_nor_all','-append')
C_nor_cue=cat(3,permute(mean(c_tmp(trial_80R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_50R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_20R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80P_ind,0*fr+1:5.5*fr,:),1),[3 2 1]));
C_nor_out=cat(3,permute(mean(c_tmp(trial_80Rwo_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80Rwx_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_50Rwo_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_50Rwx_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_20R_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80Pno_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),...
    permute(mean(c_tmp(trial_80Pnx_ind,5.5*fr+1:11*fr,:),1),[3 2 1]));

st=strsplit(cd,'\');
eval(['C_nor_cue_',st{4},'=C_nor_cue; C_nor_out_',st{4},'=C_nor_out;'])
eval(['save C_nor_',st{4},'.mat C_nor_cue_',st{4},' C_nor_out_',st{4}])

zzset_nor = [mean(C_nor_all(waterReward==1,:,:),1);mean(C_nor_all(waterReward==0,:,:),1);...
    mean(C_nor_all(trial_80Rwo_ind,:,:),1);mean(C_nor_all(trial_50Rwo_ind,:,:),1);mean(C_nor_all(trial_20Rwo_ind,:,:),1);...
    mean(C_nor_all(trial_80Rwx_ind,:,:),1);mean(C_nor_all(trial_50Rwx_ind,:,:),1);mean(C_nor_all(trial_20Rwx_ind,:,:),1)];
zzset_nor = permute(zzset_nor,[2 3 1]);
save('RPE_forz.mat','zzset_nor','-append')
%% Draw figure - zscore each neuron

cmap1= [30 45 232; 120 177 255; 125 125 125; 255 30 70]./255;
f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 6]);
hold on
% 75 Rw cue
stdshade(permute(mean(c_tmp(trial_80R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),0.3,cmap1(1,:),0*fr+1:5.5*fr,'-')
stdshade(permute(mean(c_tmp(trial_80Rwo_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(1,:),5.5*fr+1:11*fr,'-')
stdshade(permute(mean(c_tmp(trial_80Rwx_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(1,:),5.5*fr+1:11*fr,'--')
% 25 Rw cue
if sum(trial_50R_ind)>0;
    stdshade(permute(mean(c_tmp(trial_50R_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),0.3,cmap1(2,:),0*fr+1:5.5*fr,'-')
    stdshade(permute(mean(c_tmp(trial_50Rwo_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(2,:),5.5*fr+1:11*fr,'-')
    stdshade(permute(mean(c_tmp(trial_50Rwx_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(2,:),5.5*fr+1:11*fr,'--')
end
% 0 Rw cue
stdshade(permute(mean(c_tmp(trial_20R_ind,0*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(3,:),0*fr+1:11*fr,'-')
% 75 Pn cue
if sum(trial_80P_ind)>0;
    stdshade(permute(mean(c_tmp(trial_80P_ind,0*fr+1:5.5*fr,:),1),[3 2 1]),0.3,cmap1(4,:),0*fr+1:5.5*fr,'-')
    stdshade(permute(mean(c_tmp(trial_80Pno_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(4,:),5.5*fr+1:11*fr,'-')
    stdshade(permute(mean(c_tmp(trial_80Pnx_ind,5.5*fr+1:11*fr,:),1),[3 2 1]),0.3,cmap1(4,:),5.5*fr+1:11*fr,'--')
end

lim = axis;
line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.5); %0.5 baseline
line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.5); %2 cue
line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.5); %3.5 delay1
line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.5); %4 motor
line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':','linewidth',1.5); %5.5 delay2
line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.','linewidth',1.5); %8 motor off
line([1 11*fr], [0 0],'color','k','linestyle','--', 'linewidth',2);
xlim([1 11*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'Cue', 'Delay1','Delay2','Outcome'})

saveas(f2,[figure_folder,'\zscore',num2str(d),'.tif']);
%% Draw figure - zscore each neuron
mkdir([figure_folder, '\C_raw_mean']);
cmap1= [30 45 232; 120 177 255; 125 125 125; 255 30 70]./255;

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
for oo = 1:size(C_z_all,3)
    hold on
    stdshade(C_z_all(trial_80R_ind,0*fr+1:5.5*fr,oo),0.3, cmap1(1,:),0*fr+1:5.5*fr,'-')
    stdshade(C_z_all(trial_80Rwo_ind,5.5*fr+1:11*fr,oo),0.3, cmap1(1,:),5.5*fr+1:11*fr,'-')
    stdshade(C_z_all(trial_80Rwx_ind,5.5*fr+1:11*fr,oo),0.3, cmap1(1,:),5.5*fr+1:11*fr,'--')
    
    if sum(trial_50R_ind)>0;
        stdshade(C_z_all(trial_50R_ind,0*fr+1:5.5*fr,oo),0.3, cmap1(2,:),0*fr+1:5.5*fr,'-')
        stdshade(C_z_all(trial_50Rwo_ind,5.5*fr+1:11*fr,oo),0.3, cmap1(2,:),5.5*fr+1:11*fr,'-')
        stdshade(C_z_all(trial_50Rwx_ind,5.5*fr+1:11*fr,oo),0.3, cmap1(2,:),5.5*fr+1:11*fr,'--')
    end
    
    stdshade(C_z_all(trial_20R_ind,0*fr+1:11*fr,oo),0.3, cmap1(3,:),0*fr+1:11*fr,'-')
    
    if sum(trial_80P_ind)>0;
        stdshade(C_z_all(trial_80P_ind,0*fr+1:5.5*fr,oo),0.3, cmap1(4,:),0*fr+1:5.5*fr,'-')
        stdshade(C_z_all(trial_80Pno_ind,5.5*fr+1:11*fr,oo),0.3, cmap1(4,:),5.5*fr+1:11*fr,'-')
        stdshade(C_z_all(trial_80Pnx_ind,5.5*fr+1:11*fr,oo),0.3, cmap1(4,:),5.5*fr+1:11*fr,'--')
    end
    lim = axis;
    line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    xlim([0*fr 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'C', 'D1','D2','O'})
    saveas(f1,[figure_folder, '\C_raw_mean\','PSTHmean_',num2str(oo),'_craw.tif']);
    clf
end

%% Draw figure - (t-1) zscore each neuron
% for t-1 signal
trial_80R_ind1 = [false;trial_80R_ind(1:end-1)];
trial_50R_ind1 = [false;trial_50R_ind(1:end-1)];
trial_20R_ind1 = [false;trial_20R_ind(1:end-1)];
trial_80P_ind1 = [false;trial_80P_ind(1:end-1)];
trial_80Pno_ind1 = [false;trial_80Pno_ind(1:end-1)];
trial_80Pnx_ind1 = [false;trial_80Pnx_ind(1:end-1)];
mkdir([figure_folder, '\C_raw_t-1_mean']);

cmap1= [30 45 232; 120 177 255; 125 125 125; 255 30 70]./255;

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
for oo = 1:size(C_z_all,3)
    hold on
    stdshade(C_z_all(trial_80R_ind1,0*fr+1:11*fr,oo),0.3, cmap1(1,:),0*fr+1:11*fr,'-')   
    if sum(trial_50R_ind)>0;
        stdshade(C_z_all(trial_50R_ind1,0*fr+1:11*fr,oo),0.3, cmap1(2,:),0*fr+1:11*fr,'-')
    end    
    stdshade(C_z_all(trial_20R_ind1,0*fr+1:11*fr,oo),0.3, cmap1(3,:),0*fr+1:11*fr,'-')    
    if sum(trial_80P_ind)>0;
        stdshade(C_z_all(trial_80Pno_ind1,0*fr+1:11*fr,oo),0.3, cmap1(4,:),0*fr+1:11*fr,'-')
        stdshade(C_z_all(trial_80Pnx_ind1,0*fr+1:11*fr,oo),0.3, cmap1(4,:),0*fr+1:11*fr,':')
    end
    lim = axis;
    line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    xlim([0*fr 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'C', 'D1','D2','O'})
    saveas(f1,[figure_folder, '\C_raw_t-1_mean\','PSTHmean_',num2str(oo),'_craw.tif']);
    clf
end
%% correlation
%correlation rwpn
if ismember(trial_kind,[0 1 3 4]);
    for ii = 1:120;
        %rwcue-lick
        c_temp = corrcoef(lick_reg_temp(:,ii),cylinderspeed(:,ii)); c_set(1,ii) = c_temp(1,2);
        %rwcue-speed
        c_temp = corrcoef(cue1_reg(:,ii), cylinderspeed(:,ii)); c_set(2,ii) = c_temp(1,2);
        %pncue-speed
        c_temp = corrcoef(cue2_reg(:,ii), cylinderspeed(:,ii)); c_set(3,ii) = c_temp(1,2);
        %rw-speed
        c_temp = corrcoef(rew_reg(:,ii), cylinderspeed(:,ii)); c_set(4,ii) = c_temp(1,2);
        %pn-speed
        c_temp = corrcoef(pun_reg(:,ii), cylinderspeed(:,ii)); c_set(5,ii) = c_temp(1,2);
    end
    %vif
    step_num = size(cue1_reg,2);
    vif = nan(6,step_num);
    % rwcue pncue rw pn speed lick
    for istep = 1:step_num
        if sum(lick_reg(:,istep))>0
            [R,~] = corrcoef([cue1_reg(:,1), cue2_reg(:,1), rew_reg(:,1), pun_reg(:,1), cylinderspeed(:,istep),lick_reg(:,istep)]);
            temp_vif = diag(inv(R))';
            vif(1:6,istep) = temp_vif;
        else
            [R,~] = corrcoef([cue1_reg(:,1), cue2_reg(:,1), rew_reg(:,1), pun_reg(:,1),cylinderspeed(:,istep)]);
            temp_vif = diag(inv(R))';
            vif(1:5,istep) = temp_vif(1:5);
        end
    end
    % correlation-rwprob
elseif ismember(trial_kind,[0 1 2 3]);
    for ii = 1:120;
        %rwcue-lick
        c_temp = corrcoef(lick_reg_temp(:,ii),cylinderspeed(:,ii)); c_set(1,ii) = c_temp(1,2);
        %value-speed
        c_temp = corrcoef(value_reg(:,ii), cylinderspeed(:,ii)); c_set(2,ii) = c_temp(1,2);
        %outcome-speed
        c_temp = corrcoef(outcome_reg(:,ii), cylinderspeed(:,ii)); c_set(3,ii) = c_temp(1,2);
    end
    %vif
    step_num = size(value_reg,2);
    vif = nan(4,step_num);
    % value outcome speed lick
    for istep = 1:step_num
        if sum(lick_reg(:,istep))>0
            [R,~] = corrcoef([value_reg(:,1), outcome_reg(:,1), cylinderspeed(:,istep),lick_reg(:,istep)]);
            temp_vif = diag(inv(R))';
            vif(1:4,istep) = temp_vif;
        else
            [R,~] = corrcoef([value_reg(:,1), outcome_reg(:,1), cylinderspeed(:,istep)]);
            temp_vif = diag(inv(R))';
            vif(1:3,istep) = temp_vif(1:3);
        end
    end
end
save('CorrVIFsetmat','c_set','vif')
%% corrlation vif figure
day = datetime('today');
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7.5 6]);
if ismember(trial_kind,[0 1 3 4]);
    cmap = [80 80 80; 129 212 250; 239 154 154; 0 0 205; 255 30 70]./255;
elseif ismember(trial_kind,[0 1 2 3]);
    cmap = [80 80 80; 199 21 133; 34 139 34]./255;
end
hold on
for ii = 1:size(c_set,1)
    plot(c_set(ii,:),'Color',cmap(ii,:))
end
lim = [0 0 -1 1];
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([0 120], [0.8 0.8],'color','k', 'linestyle', '--'); % 0.8 lim
xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'});xlim([1 120]);
ylim([-1 1]); ylabel('Correlation')
saveas(f1,[cd,'\figure\','Corr',d,'.tif'])


f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7.5 6]);
if ismember(trial_kind,[0 1 3 4]);
    cmap = [129 212 250; 239 154 154; 0 0 205; 255 30 70; 160 82 45; 0.3 0.3 0.3]./255;
elseif ismember(trial_kind,[0 1 2 3]);
    cmap = [199 21 133; 34 139 34; 160 82 45;0.3 0.3 0.3]./255;
end

for ii = 1:size(vif,1)
    hold on
    plot(vif(ii,:)','Color',cmap(ii,:))
end
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([0 120], [5 5],'color','k', 'linestyle', '--'); % 0.8 lim
xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'});xlim([1 120]);
ylabel('VIF')
saveas(f1,[cd,'\figure\','VIF',d,'.tif'])

%% multiple regression part 1- Craw
[X2,Y2] = meshgrid(-0.5*fr:12.5*fr-1,GPIO1onframe(set_rectrial,1)); %-0.25 ~ 5.25 ; 6~105 frame
set_regressionframe = X2+Y2;
% lick for one
lickinfo = timeline{1,1};
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms
lick_reg = lick_reg(set_rectrial,:);
% lick for rev
% clear licktemp lick_reg_temp
% lickinfo = [timeline{1,1};timeline{1,2}];
% licksort = sortrows(lickinfo(:,1:3), 2);
% for iii = 1:licksort(end,2);
%     licktemp = licksort((licksort(:,2)==iii),1);
%     lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
% end
% lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms
% lick_reg = lick_reg(set_rectrial,:);
% multiple regression part2- Craw
warning off;
clear rew_regTmp pun_regTmp value_regTmp
set_rectrial= 1:420;
%cue, reward, punish
outcome = waterReward;
% outcome(waterReward==0 & ismember(odorCue, find(trial_kind==1|trial_kind==2|trial_kind==4)-1)) =-1;
for ii = 1:size(odorCue,1);
    rew_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==1)& (rwProb(ii,odorCue(ii)+1)==75);
    pun_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==2)& (rwProb(ii,odorCue(ii)+1)==75);
    value_regTmp(ii,1) = rwProb(ii,odorCue(ii)+1)./100;    
end
set_crp = [odorCue, outcome.*rew_regTmp', outcome.*pun_regTmp'];
set_crp = set_crp(set_rectrial,:);
if ismember(trial_kind,[0 1 3 4]);
    % rp reg
    cue1_reg = repmat(set_crp(:,1)==(find(trial_kind==1)-1),1,120);% Rw
    cue2_reg = repmat(set_crp(:,1)==(find(trial_kind==4)-1),1,120);% Pn
    cue1_t1_reg = repmat([0; cue1_reg(1:end-1,1)],1,120);% Rw
    cue2_t1_reg = repmat([0; cue2_reg(1:end-1,1)],1,120);% Pn
    cue1_t2_reg = repmat([0;0; cue1_reg(1:end-2,1)],1,120);% Rw
    cue2_t2_reg = repmat([0;0; cue2_reg(1:end-2,1)],1,120);% Pn
    
    rew_reg = repmat(set_crp(:,2),1,120);
    pun_reg = repmat(set_crp(:,3),1,120);
    rew_1_reg =repmat([0; set_crp(1:end-1,2)],1,120);
    pun_1_reg =repmat([0; set_crp(1:end-1,3)],1,120);
    rew_2_reg =repmat([0;0; set_crp(1:end-2,2)],1,120);
    pun_2_reg =repmat([0;0; set_crp(1:end-2,3)],1,120);
elseif ismember(trial_kind,[0 1 2 3]);
    % prob reg
    outcome = waterReward;
%     outcome(waterReward==0)=-1;
%     outcome(trial_20R_ind)=0;
    set_crp = [outcome, value_regTmp];
    set_crp = set_crp(set_rectrial,:);
    outcome_reg = repmat(set_crp(:,1), 1, 120);
    outcome_1_reg = repmat([0; set_crp(1:end-1,1)],1,120);
    outcome_2_reg = repmat([0;0; set_crp(1:end-2,1)],1,120);
    value_reg = repmat(set_crp(:,2), 1, 120);
    value_1_reg = repmat([0; set_crp(1:end-1,2)],1,120);
    value_2_reg = repmat([0;0; set_crp(1:end-2,2)],1,120);
end
%cut with neuron event
C_raw_reg = C_raw0(neuron_over15,:);

% C_raw_reg = S_raw0(neuron_over15,:);
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset
%% multiple regression - C raw all RP
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset
for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1); % trial onset 0~12s
    for tt = 1:120       
         tbl = table(cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),cylinderspeed(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),...            
            mean(craw_reg(:,3*tt-2:3*tt),2),...% craw_reg(:,3*tt-1),...   
            'VariableNames',{'Cue1','Cue2','RW','PN','Lick','Speed','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
        mdl = fitlm(tbl,...
            'C_raw~Cue1+RW+Cue2+PN+Lick+Speed+Cue1t_1+Rwt_1+Cue2t_1+Pnt_1+Cue1t_2+Rwt_2+Cue2t_2+Pnt_2'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([cue1_reg(:,tt),cue2_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),cylinderspeed(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)]);
        ystd(1,tt) = std(mean(craw_reg(:,3*tt-2:3*tt),2));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        X_svm(:,oo,tt) = mean(craw_reg(:,3*tt-2:3*tt),2);
        Y_svm(:,:,tt) = [cue1_reg(:,tt)+cue2_reg(:,tt)*2, rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),cylinderspeed(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)];
    end
    srcset(:,:,oo) = src;
    betaset(:,:,oo) = beta;
    pvalueset(:,:,oo) = pvalue;
    disp(oo)
end
save ('all_reg_20220220.mat','beta', 'pvalue','xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')
%% multiple regression - C raw all Rwprob
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset
for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    for tt = 1:120       
         tbl = table(outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt), ...
             lick_reg(:,tt),cylinderspeed(:,tt),mean(craw_reg(:,3*tt-2:3*tt),2),...
            'VariableNames',{'o','ot_1','ot_2','v','vt_1','vt_2','Lick','speed','C_raw'});
        mdl = fitlm(tbl,'C_raw~o+ot_1+ot_2+v+vt_1+vt_2+Lick+speed'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([outcome_reg(:,tt), outcome_1_reg(:,tt), outcome_2_reg(:,tt), value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt),lick_reg(:,tt),cylinderspeed(:,tt)]);
        ystd(1,tt) = std(mean(craw_reg(:,3*tt-2:3*tt),2));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        X_svm(:,oo,tt) = mean(craw_reg(:,3*tt-2:3*tt),2);
        Y_svm(:,:,tt) = [outcome_reg(:,tt), outcome_1_reg(:,tt), outcome_2_reg(:,tt), value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt),lick_reg(:,tt),cylinderspeed(:,tt)];
    end
    srcset(:,:,oo) = src;
    betaset(:,:,oo) = beta;
    pvalueset(:,:,oo) = pvalue;
    disp(oo)
end
save ('all_reg_20220220.mat','beta', 'pvalue','xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')

%% multiple regression - Craw before
%before rev
[X2,Y2] = meshgrid(-0.5*fr:12.5*fr-1,GPIO1onframe(1:sum(rev_b))); %-0.25 ~ 5.25 ; 6~105 frame
set_regressionframe = X2+Y2;

%lick
lickinfo = [timeline{1,1};timeline{1,2}];
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms

%cue, reward, punish
outcome = waterReward;
outcome(waterReward==0 & ismember(odorCue, [0 1]))=-1;
for ii = 1:size(odorCue,1);
    rew_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==1)& (rwProb(ii,odorCue(ii)+1)==75);
    pun_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==2)& (rwProb(ii,odorCue(ii)+1)==75);
    value_regTmp(ii,1) = rwProb(ii,odorCue(ii)+1)./100;
end
if ismember(trial_kind,[0 1 3 4]);
    set_crp = [odorCue, outcome.*rew_regTmp', outcome.*pun_regTmp'];
    set_crp = set_crp(set_rectrial(1:sum(rev_b)),:);
    cue1_reg = repmat(set_crp(:,1)==(find(trial_kind==1)-1),1,120); % Rw
    cue2_reg = repmat(set_crp(:,1)==(find(trial_kind==4)-1),1,120); % Pn
    cue1_t1_reg = repmat([0; cue1_reg(1:end-1,1)],1,120);% Rw
    cue2_t1_reg = repmat([0; cue2_reg(1:end-1,1)],1,120);% Pn
    cue1_t2_reg = repmat([0;0; cue1_reg(1:end-2,1)],1,120);% Rw
    cue2_t2_reg = repmat([0;0; cue2_reg(1:end-2,1)],1,120);% Pn
    
    rew_reg = repmat(set_crp(:,2),1,120);
    pun_reg = repmat(set_crp(:,3),1,120);
    rew_1_reg =repmat([0; set_crp(1:end-1,2)],1,120);
    pun_1_reg =repmat([0; set_crp(1:end-1,3)],1,120);
    rew_2_reg =repmat([0;0; set_crp(1:end-2,2)],1,120);
    pun_2_reg =repmat([0;0; set_crp(1:end-2,3)],1,120);
elseif ismember(trial_kind,[0 1 2 3]);
    set_crp = [waterReward, value_regTmp];
    set_crp = set_crp(set_rectrial(1:sum(rev_b)),:);
    outcome_reg = repmat(set_crp(:,1), 1, 120);
    outcome_1_reg = repmat([0; set_crp(1:end-1,1)],1,120);
    value_reg = repmat(set_crp(:,2), 1, 120);
    value_1_reg = repmat([0; set_crp(1:end-1,2)],1,120);
    
end

lick_reg = lick_reg(set_rectrial(1:sum(rev_b)),:);


%cut with neuron event
C_raw = C_raw0;
C_raw_reg = C_raw(neuron_over15,:);
%C_raw_reg = S_raw0(neuron_over25,:);
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset

%% with lick
for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    
    for tt = 1:120
        %mdl = fitlm([cue_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), mean(lick_reg(:,2*tt-1:2*tt),2)], mean(craw_reg(:,2*tt-1:2*tt),2)); % x1:cue, x2: rw, x3:pn, x4:lick
        tbl = table(cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),mean(craw_reg(:,3*tt-2:3*tt),2),...
            'VariableNames',{'Cue1','Cue2','RW','PN','Lick','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
        mdl = fitlm(tbl,'C_raw~Cue1+Cue2+RW+PN+Lick+Cue1t_1+Cue2t_1+Rwt_1+Pnt_1+Cue1t_2+Cue2t_2+Rwt_2+Pnt_2');
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([cue1_reg(:,tt),cue2_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)]);
        ystd(1,tt) = std(mean(craw_reg(:,3*tt-2:3*tt),2));
        src(:,tt) = beta(2:14,tt).*xstd(:,tt)./ystd(1,tt);
        X_svm(:,oo,tt) = mean(craw_reg(:,3*tt-2:3*tt),2);
        Y_svm(:,:,tt) = [cue1_reg(:,tt)+cue2_reg(:,tt)*2, rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)];
    end
    srcset(:,:,oo) = src;
    betaset(:,:,oo) = beta;
    pvalueset(:,:,oo) = pvalue;
    disp(oo)
end
save ('before_re_C15_withR(t-1)P(t-1)C(t-1)&(t-2)g.mat','beta', 'pvalue','xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')


%% multiple regression - Craw after
%after rev
[X2,Y2] = meshgrid(-0.5*fr:12.5*fr-1,GPIO1onframe(sum(rev_b)+1 : end)); %-0.25 ~ 5.25 ; 6~105 frame
set_regressionframe = X2+Y2;

%lick
lickinfo = [timeline{1,1};timeline{1,2}];
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms

%cue, reward, punish
outcome = waterReward;
outcome(waterReward==0 & ismember(odorCue, [0 1]))=-1;
for ii = 1:size(odorCue,1);
    rew_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==1)& (rwProb(ii,odorCue(ii)+1)==75);
    pun_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==2)& (rwProb(ii,odorCue(ii)+1)==75);
end
set_crp = [odorCue, outcome.*rew_regTmp', outcome.*pun_regTmp'];
set_crp = set_crp(set_rectrial(sum(rev_b)+1:end),:);

cue1_reg = repmat(set_crp(:,1)==(find(trial_kind==1)-1),1,120); % Rw
cue2_reg = repmat(set_crp(:,1)==(find(trial_kind==4)-1),1,120); % Pn
cue1_t1_reg = repmat([0; cue1_reg(1:end-1,1)],1,120);% Rw
cue2_t1_reg = repmat([0; cue2_reg(1:end-1,1)],1,120);% Pn
cue1_t2_reg = repmat([0;0; cue1_reg(1:end-2,1)],1,120);% Rw
cue2_t2_reg = repmat([0;0; cue2_reg(1:end-2,1)],1,120);% Pn

rew_reg = repmat(set_crp(:,2),1,120);
pun_reg = repmat(set_crp(:,3),1,120);
rew_1_reg =repmat([0; set_crp(1:end-1,2)],1,120);
pun_1_reg =repmat([0; set_crp(1:end-1,3)],1,120);
rew_2_reg =repmat([0;0; set_crp(1:end-2,2)],1,120);
pun_2_reg =repmat([0;0; set_crp(1:end-2,3)],1,120);


lick_reg = lick_reg(set_rectrial(sum(rev_b)+1:end),:);

%cut with neuron event
C_raw = C_raw0;
C_raw_reg = C_raw(neuron_over15,:);
%C_raw_reg = S_raw0(neuron_over25,:);
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset

for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    
    for tt = 1:120
        tbl = table(cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),mean(craw_reg(:,3*tt-2:3*tt),2),...
            'VariableNames',{'Cue1','Cue2','RW','PN','Lick','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
        mdl = fitlm(tbl,'C_raw~Cue1+Cue2+RW+PN+Lick+Cue1t_1+Cue2t_1+Rwt_1+Pnt_1+Cue1t_2+Cue2t_2+Rwt_2+Pnt_2');
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([cue1_reg(:,tt),cue2_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)]);
        ystd(1,tt) = std(mean(craw_reg(:,3*tt-2:3*tt),2));
        src(:,tt) = beta(2:14,tt).*xstd(:,tt)./ystd(1,tt);
        X_svm(:,oo,tt) = mean(craw_reg(:,3*tt-2:3*tt),2);
        Y_svm(:,:,tt) = [cue1_reg(:,tt)+cue2_reg(:,tt)*2, rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)];
    end
    srcset(:,:,oo) = src;
    betaset(:,:,oo) = beta;
    pvalueset(:,:,oo) = pvalue;
    disp(oo)
end
save ('after_re_C15_withR(t-1)P(t-1)C(t-1)&(t-2)g.mat','beta', 'pvalue','xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')
%% Regression shuffle_craw
clear beta_r pvalue_r xstd_r ystd_r src_r srcset_r betaset_r pvalueset_r
for rant = 1:1000
    for oo = 1:size(C_raw_reg,1)
        regression_movmean_r = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
        craw_reg_r = regression_movmean_r(:,0.5*fr:12.5*fr-1);
        spike_rand = randperm(size(cue1_reg,1));
        for tt = 1:120
            tbl = table(cue1_reg(spike_rand,tt), cue2_reg(spike_rand,tt),rew_reg(spike_rand,tt), pun_reg(spike_rand,tt), lick_reg(spike_rand,tt),...
                rew_1_reg(spike_rand,tt),pun_1_reg(spike_rand,tt),cue1_t1_reg(spike_rand,tt),cue2_t1_reg(spike_rand,tt),...
                rew_2_reg(spike_rand,tt),pun_2_reg(spike_rand,tt),cue1_t2_reg(spike_rand,tt),cue2_t2_reg(spike_rand,tt),mean(craw_reg_r(:,3*tt-2:3*tt),2),...
                'VariableNames',{'Cue1','Cue2','RW','PN','Lick','Rwt_1','Pnt_1','Cue1t_1','Cue2t_1','Rwt_2','Pnt_2','Cue1t_2','Cue2t_2','C_raw'});
            mdl_r = fitlm(tbl,'C_raw~Cue1+Cue2+RW+PN+Lick+Rwt_1+Pnt_1+Cue1t_1+Cue2t_1+Rwt_2+Pnt_2+Cue1t_2+Cue2t_2');
            beta_r(:,tt) = mdl_r.Coefficients.Estimate;
            pvalue_r(:,tt) = mdl_r.Coefficients.pValue;
            xstd_r(:,tt) = std([cue1_reg(:,tt),cue2_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
                rew_1_reg(:,tt),pun_1_reg(:,tt),cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),...
                rew_2_reg(:,tt),pun_2_reg(:,tt),cue1_t2_reg(:,tt),cue2_t2_reg(:,tt)]);
            ystd_r(1,tt) = std(mean(craw_reg_r(:,3*tt-2:3*tt),2));
            src_r(:,tt) = beta_r(2:14,tt).*xstd_r(:,tt)./ystd_r(1,tt);
        end
        srcset_r(:,:,oo,rant) = src_r;
        betaset_r(:,:,oo,rant) = beta_r;
        pvalueset_r(:,:,oo,rant) = pvalue_r;
    end
    disp(rant)
end
save ('all_reg_shuffle_C15_with(t-1)&(t-2).mat','srcset_r','betaset_r','pvalueset_r')

pvaluesig_r = pvalueset_r(2:14,:,:,:)<0.05;
FON_r = sum(pvaluesig_r,3)./size(C_raw,1);
FON_r = permute(FON_r,[1 2 4 3]);

%% Draw FON - 1.calculate
cmap1 = [30 45 232; 120 177 255; 232 126 58; 255 30 70]/255; %blue skyblue pink red
pvaluesig = pvalueset(2:end,:,:)<0.05;
FON = sum(pvaluesig, 3)./size(pvalueset,3);
neuronum = size(pvalueset,3); % binomial test --> make function
neuron_sigbinomial = binoinv([0.05 0.95], neuronum, 0.05);
bin_percent = neuron_sigbinomial(2)./neuronum;
if exist('FON_r')
    FON_compare = FON>FON_r(:,:,:);
    FON_sig = sum(FON_compare,3)>95;
    FON_sig = double(FON_sig);
    FON_sig(FON_sig==0)=nan;
end
%% Draw FON - 2. RP
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 30 9]);
subplot(1,3,1);
hold on
p1=plot(FON(1,:), 'color', [129 212 250]./255, 'linewidth', 1);%cue1
p2=plot(FON(2,:), 'color', [239 154 154]./255, 'linewidth', 1);%cue2
p3=plot(FON(3,:), 'color', [0 0 205]./255,'linewidth', 1);%Rw
p4=plot(FON(4,:), 'color', [255 30 70]./255, 'linewidth', 1);%Pn
p5=plot(FON(5,:), 'color', [80 80 80]./255,      'linewidth', 1);%Lick
ylim([0 1])
xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
title('t')
legend([p1,p2,p3,p4,p5],'75Rw Cue', '75Pn Cue','Rw','Pn','Lick');

subplot(1,3,2)
hold on
p6=plot(FON(7,:), 'color', [129 212 250]./255,'linewidth', 1);%cue1(Rwcue t-1)
p7=plot(FON(8,:), 'color', [239 154 154]./255,'linewidth', 1);%cue2(Pncue t-1)
p8=plot(FON(9,:), 'color', [0 0 205]./255,'linewidth', 1);%Rw(t-1)
p9=plot(FON(10,:), 'color', [255 30 70]./255,'linewidth', 1);%Pn(t-1)
ylim([0 1])
xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
title('t-1')
legend([p6,p7,p8,p9],'Rw(t-1)','Pn(t-1)','Rw Cue(t-1)','Pn Cue(t-1)')

subplot(1,3,3)
hold on
p10=plot(FON(11,:), 'color', [129 212 250]./255,'linewidth', 1);%cue1(Rwcue t-1)
p11=plot(FON(12,:), 'color', [239 154 154]./255,'linewidth', 1);%cue2(Pncue t-1)
p12=plot(FON(13,:), 'color', [0 0 205]./255,'linewidth', 1);%Rw(t-1)
p13=plot(FON(14,:), 'color', [255 30 70]./255,'linewidth', 1);%Pn(t-1)
ylim([0 1])
xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
title('t-2')
legend([p10,p11,p12,p13],'Rw(t-2)','Pn(t-2)','Rw Cue(t-2)','Pn Cue(t-2)')

saveas(f1,[cd,'\figure\','FON_C15_with(t-1)(t-2)_nomean',d,'.tif'])

%% Draw FON - 2. Prob
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 9]);
subplot(1,2,1)
hold on
p1=plot(FON(1,:), 'color', [34 139 34]./255, 'linewidth', 1);% o
p3=plot(FON(3,:), 'color', [199 21 133]./255,'linewidth', 1);% v
p5=plot(FON(5,:), 'color', [80 80 80]./255,  'linewidth', 1);% Lick
p6=plot(FON(6,:), 'color', [210 105 30]./255,  'linewidth', 1); % speed
ylim([0 1])
xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
title('t')
% legend([p1,p3,p5],'Outcome(t)', 'Value(t)','Lick');

subplot(1,2,2)
hold on

p2=plot(FON(2,:), 'color', [34 139 34]./255, 'linewidth', 1);% o(t-1)
p4=plot(FON(4,:), 'color', [199 21 133]./255, 'linewidth', 1);% v(t-1)


ylim([0 1])
xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
title('t-1')
% legend([p2,p4],'Outcome(t-1)','Value(t-1)');
saveas(f1,[cd,'\figure\','FON_C15_with(t-1)withinx0-1mix',d,'.tif'])
%% Draw SRC all
for oo = 1:size(C_raw,1)
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 5]);
    hold on
    xlim([1 50]); xticks([10 20 30]); xticklabels({'Cue', 'Delay','Outcome'})
    for ii = 1:13
    plot(srcset(ii,:,oo))
    end
    lim = axis;
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([15 15], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([25 25], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    print(f1,'-dtiff','-r600',[figure_folder,'\','SRC_',num2str(oo),'_craw.tif']);
    close all
end

%% Draw SRC mean - rwpn
shift_src =permute(abs(srcset),[3 2 1]);
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 30 9]);
subplot(1,3,1);
hold on
stdshade(shift_src(:,:,1), 0.3, [129 212 250]./255)%rwcue
stdshade(shift_src(:,:,2), 0.3, [239 154 154]./255)%pncue
stdshade(shift_src(:,:,3), 0.3, [0 0 205]./255)%rw
stdshade(shift_src(:,:,4), 0.3, [255 30 70]./255)%pn
stdshade(shift_src(:,:,5), 0.3, [80 80 80]./255)%Lick
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off

subplot(1,3,2);
hold on
stdshade(shift_src(:,:,6), 0.3, [129 212 250]./255)%rwcue
stdshade(shift_src(:,:,7), 0.3, [239 154 154]./255)%pncue
stdshade(shift_src(:,:,8), 0.3, [0 0 205]./255)%rw
stdshade(shift_src(:,:,9), 0.3, [255 30 70]./255)%pn
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome

subplot(1,3,3);
hold on
stdshade(shift_src(:,:,10), 0.3, [129 212 250]./255)%rwcue
stdshade(shift_src(:,:,11), 0.3, [239 154 154]./255)%pncue
stdshade(shift_src(:,:,12), 0.3, [0 0 205]./255)%rw
stdshade(shift_src(:,:,13), 0.3, [255 30 70]./255)%pn
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome

hold off
saveas(f1,[cd,'\figure\','SRC_C15_withR(t-1)P(t-1)C(t-1)&(t-2)',d,'.tif'])
%% Draw SRC mean - rwprob
shift_src =permute(abs(srcset),[3 2 1]);
% [h1,p] = ttest2(abs(srcset(1,:,:)),mean(abs(srcset_r(1,:,:,:)),4),'Dim',3);
% [h2,p] = ttest2(abs(srcset(2,:,:)),mean(abs(srcset_r(2,:,:,:)),4),'Dim',3);
% [h3,p] = ttest2(abs(srcset(3,:,:)),mean(abs(srcset_r(3,:,:,:)),4),'Dim',3);
% [h4,p] = ttest2(abs(srcset(4,:,:)),mean(abs(srcset_r(4,:,:,:)),4),'Dim',3);
% h1(h1==0)=nan;h2(h2==0)=nan;h3(h3==0)=nan;h4(h4==0)=nan;
% srcset_mean = mean(mean(abs(srcset_r(:,:,:,:)),4),3);
% plot(h1*0.195,'color', [50 50 50]./255, 'linewidth', 1);%cue sig
% plot(h2*0.19,'color', cmap1(1,:),      'linewidth', 1);%RW sig
% plot(h3*0.185,'color', cmap1(4,:),      'linewidth', 1);%PN sig
% plot(h4*0.18,'color', [98 173 88]./255,'linewidth', 1);%lick sig
% plot(srcset_mean(1,:), 'color', [50 50 50]./255, 'linewidth', 1.5);
% plot(srcset_mean(2,:), 'color', cmap1(1,:), 'linewidth', 1.5);
% plot(srcset_mean(3,:), 'color', cmap1(4,:), 'linewidth', 1.5);
% plot(srcset_mean(4,:), 'color', [98 173 88]./255, 'linewidth', 1.5);

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 30 9]);
hxa(1) = subplot(1,2,1);
xlim(hxa, [1 120]);
ylim(hxa, [0 0.6]);
hold on
stdshade(shift_src(:,:,5), 0.3, [80 80 80]./255)%Lick
stdshade(shift_src(:,:,6), 0.3, [160 82 45]./255)%speed
stdshade(shift_src(:,:,1), 0.3, [34 139 34]./255)%o
stdshade(shift_src(:,:,3), 0.3, [199 21 133]./255)%v

line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off

hxa(2) = subplot(1,2,2);
xlim(hxa, [1 120]);
ylim(hxa, [0 0.6]);
hold on
stdshade(shift_src(:,:,2), 0.3, [34 139 34]./255)%o t-1
stdshade(shift_src(:,:,4), 0.3, [199 21 133]./255)%v t-1
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome


hold off
saveas(f1,[cd,'\figure\','SRC_C15_with(t-1)v',d,'.tif'])

%% Coefficient of patial determinant - Rwpn
clear sse_set
warning off
alist = 1:14;
% rn = randperm(sum(trial_80R_ind));
for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    
    for tt = 1:120
        %mdl = fitlm([cue_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), mean(lick_reg(:,2*tt-1:2*tt),2)], mean(craw_reg(:,2*tt-1:2*tt),2)); % x1:cue, x2: rw, x3:pn, x4:lick
        tbl = [cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),cylinderspeed(:,tt),...
             mean(craw_reg(:,3*tt-2:3*tt),2)];
%         tbl = tbl(trial_80R_ind,:);
        mdl = fitlm(tbl(:,1:end-1),tbl(:,end));
        sse_set(1,tt,oo) = mdl.SSE;
        for ii = 1:size(tbl,2)-1
            a = alist; a(ii)=[];
            mdl = fitlm(tbl(:,a),tbl(:,end));
            sse_set(1+ii,tt,oo) = mdl.SSE;
            antn{ii,:} = mdl;
        end
    end
    disp(oo)
end
% save ('CPD.mat','sse_set')

%% Coefficient of patial determinant - period
clear sse_set

alist = 1:14;
shuffletrial = randperm(420);
cyliderspeedset = [mean(cylinderspeed(:,0.5*10+1:2*10+1),2),mean(cylinderspeed(:,2*10+1:3.5*10+1),2)...
    mean(cylinderspeed(:,0.5*10+1:2*10+1),2),mean(cylinderspeed(:,2*10+1:3.5*10+1),2),...
    mean(cylinderspeed(:,5.5*10+1:8*10+1),2),mean(cylinderspeed(:,5.5*10+1:8*10+1),2),...
    mean(cylinderspeed(:,5.5*10+1:6.5*10+1),2),mean(cylinderspeed(:,5.5*10+1:6.5*10+1),2),...
    mean(cylinderspeed(:,6.5*10+1:8*10+1),2),mean(cylinderspeed(:,6.5*10+1:8*10+1),2)];
lickset = [mean(lick_reg(:,0.5*10+1:2*10+1),2),mean(lick_reg(:,2*10+1:3.5*10+1),2)...
    mean(lick_reg(:,0.5*10+1:2*10+1),2),mean(lick_reg(:,2*10+1:3.5*10+1),2),...
    mean(lick_reg(:,5.5*10+1:8*10+1),2),mean(lick_reg(:,5.5*10+1:8*10+1),2),...
    mean(lick_reg(:,5.5*10+1:6.5*10+1),2),mean(lick_reg(:,5.5*10+1:6.5*10+1),2),...
    mean(lick_reg(:,6.5*10+1:8*10+1),2),mean(lick_reg(:,6.5*10+1:8*10+1),2)];
Oshuffle = [6 8 10];
for ios = Oshuffle
rew_reg(:,ios) = rew_shu;
pun_reg(:,ios) = pun_shu;
end
for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    craw_cue = mean(craw_reg(:,0.5*fr+1:2*fr+1),2);
    craw_delay = mean(craw_reg(:,2*fr+1:3.5*fr+1),2);
    craw_outcome = mean(craw_reg(:,5.5*fr+1:8*fr+1),2);
    craw_outcome1s = mean(craw_reg(:,5.5*fr+1:8*fr+1),2);
    craw_outcome2s = mean(craw_reg(:,5.5*fr+1:8*fr+1),2);
    craw_cue_shuffle = mean(craw_reg(shuffletrial,0.5*fr+1:2*fr+1),2);
    craw_delay_shuffle = mean(craw_reg(shuffletrial,2*fr+1:3.5*fr+1),2);
    craw_set = [craw_cue,craw_delay,craw_cue_shuffle,craw_delay_shuffle,...
        craw_outcome,craw_outcome,craw_outcome1s,craw_outcome1s,craw_outcome2s,craw_outcome2s];

    
    for tt = 1:10
        %mdl = fitlm([cue_reg(:,tt), rew_reg(:,tt), pun_reg(:,tt), mean(lick_reg(:,2*tt-1:2*tt),2)], mean(craw_reg(:,2*tt-1:2*tt),2)); % x1:cue, x2: rw, x3:pn, x4:lick
        tbl = [cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lick_reg(:,tt),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),...
            cylinderspeed(:,tt),craw_set(:,tt)];   
        mdl = fitlm(tbl(:,1:end-1),tbl(:,end));
        sse_set(1,tt,oo) = mdl.SSE;
        for ii = 1:size(tbl,2)-1
            a = alist; a(ii)=[];
            mdl = fitlm(tbl(:,a),tbl(:,end));
            sse_set(1+ii,tt,oo) = mdl.SSE;
            antn{ii,:} = mdl;
        end
    end
    disp(oo)
end
save ('CPD_period.mat','sse_set')
%% Draw CPD - Rwpn
day = datetime('today');
d = datestr(day,'yymmdd');
cpd = (sse_set(2:end,:,:)-sse_set(1,:,:))./sse_set(2:end,:,:);
cpd = permute(cpd,[3 2 1]);
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 25 6]);
hxa(1) = subplot(1,3,1);

ylimset = [0 0.06];

hold on
% stdshade(cpd(:,:,14), 0.3, [160 82 45]./255)%speed
stdshade(cpd(:,:,5), 0.3, [80 80 80]./255)%Lick
stdshade(cpd(:,:,4), 0.3, [255 30 70]./255)%pn
stdshade(cpd(:,:,3), 0.3, [0 0 205]./255)%rw
stdshade(cpd(:,:,1), 0.3, [129 212 250]./255)%rwcue
stdshade(cpd(:,:,2), 0.3, [239 154 154]./255)%pncue

ylim(ylimset);lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off

hxa(2) = subplot(1,3,2);
hold on
stdshade(cpd(:,:,6), 0.3, [129 212 250]./255)%rwcue
stdshade(cpd(:,:,7), 0.3, [239 154 154]./255)%pncue
stdshade(cpd(:,:,8), 0.3, [0 0 205]./255)%rw
stdshade(cpd(:,:,9), 0.3, [255 30 70]./255)%pn
ylim(ylimset);lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off

hxa(3) = subplot(1,3,3);
hold on
stdshade(cpd(:,:,10), 0.3, [129 212 250]./255)%rwcue
stdshade(cpd(:,:,11), 0.3, [239 154 154]./255)%pncue
stdshade(cpd(:,:,12), 0.3, [0 0 205]./255)%rw
stdshade(cpd(:,:,13), 0.3, [255 30 70]./255)%pn
xlim(hxa,[1 120]);
ylim(ylimset); lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %-0.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off
xticks(hxa,[12 27 47 67 100]); xticklabels(hxa, {'C', 'D1','D2','O','ITI'})
saveas(f1,[cd,'\figure\','CPD_withv',d,'.tif'])
%% Coefficient of patial determinant - period
clear sse_set
warning off
shuffletrial = randperm(420);
alist = 1:8;
value_2 = value_reg(:,1); value_2(value_reg(:,1)==0)=0.25; value_2(value_reg(:,1)==0.25)=0; %change 0-25 sort
value_3 = value_reg(:,1); value_3(value_reg(:,1)==0.25)=0.75; value_3(value_reg(:,1)==0.75)=0.25; %change 0-25 sort
value_set = [value_reg(:,1),value_reg(:,1),value_reg(:,1),value_reg(:,1),value_2,value_3,value_2,value_3,...
    value_reg(:,1),value_reg(:,1),value_reg(:,1),value_reg(:,1),value_reg(:,1),value_reg(:,1)];
value_1_set = [0 0 0 0 0 0 0 0 0 0 0 0 0 0; value_set(1:end-1,:)];
value_2_set = [0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0; value_set(1:end-2,:)];
cyliderspeedset = [mean(cylinderspeed(:,0.5*10+1:2*10+1),2),mean(cylinderspeed(:,2*10+1:3.5*10+1),2),...
    mean(cylinderspeed(:,0.5*10+1:2*10+1),2),mean(cylinderspeed(:,2*10+1:3.5*10+1),2),...
    mean(cylinderspeed(:,0.5*10+1:2*10+1),2),mean(cylinderspeed(:,2*10+1:3.5*10+1),2),...
    mean(cylinderspeed(:,0.5*10+1:2*10+1),2),mean(cylinderspeed(:,2*10+1:3.5*10+1),2),...
    mean(cylinderspeed(:,5.5*10+1:8*10+1),2),mean(cylinderspeed(:,5.5*10+1:8*10+1),2),...
    mean(cylinderspeed(:,5.5*10+1:6.5*10+1),2),mean(cylinderspeed(:,5.5*10+1:6.5*10+1),2),...
    mean(cylinderspeed(:,6.5*10+1:8*10+1),2),mean(cylinderspeed(:,6.5*10+1:8*10+1),2)];
lickset = [mean(lick_reg(:,0.5*10+1:2*10+1),2),mean(lick_reg(:,2*10+1:3.5*10+1),2),...
    mean(lick_reg(:,0.5*10+1:2*10+1),2),mean(lick_reg(:,2*10+1:3.5*10+1),2),...
    mean(lick_reg(:,0.5*10+1:2*10+1),2),mean(lick_reg(:,2*10+1:3.5*10+1),2),...
    mean(lick_reg(:,0.5*10+1:2*10+1),2),mean(lick_reg(:,2*10+1:3.5*10+1),2),...
    mean(lick_reg(:,5.5*10+1:8*10+1),2),mean(lick_reg(:,5.5*10+1:8*10+1),2),...
    mean(lick_reg(:,5.5*10+1:6.5*10+1),2),mean(lick_reg(:,5.5*10+1:6.5*10+1),2),...
    mean(lick_reg(:,6.5*10+1:8*10+1),2),mean(lick_reg(:,6.5*10+1:8*10+1),2)];
for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    craw_cue = mean(craw_reg(:,0.5*fr+1:2*fr+1),2);
    craw_delay = mean(craw_reg(:,2*fr+1:3.5*fr+1),2);
    craw_outcome = mean(craw_reg(:,5.5*fr+1:8*fr+1),2);
    craw_outcome1 = mean(craw_reg(:,5.5*fr+1:6.5*fr+1),2);
    craw_outcome2 = mean(craw_reg(:,6.5*fr+1:8*fr+1),2);
    craw_cue_shuffle = mean(craw_reg(shuffletrial,0.5*fr+1:2*fr+1),2);
    craw_delay_shuffle = mean(craw_reg(shuffletrial,2*fr+1:3.5*fr+1),2);
    craw_set = [craw_cue,craw_delay,craw_cue_shuffle,craw_delay_shuffle,craw_cue,craw_delay,craw_cue,craw_delay,...
        craw_outcome,craw_outcome,craw_outcome1,craw_outcome,craw_outcome2,craw_outcome];
    
    for tt =13:14
        tbl = [outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_set(:,tt),value_1_set(:,tt),value_2_set(:,tt),cylinderspeed(:,tt),...
             lickset(:,tt),craw_set(:,tt)];
        mdl = fitlm(tbl(:,1:end-1),tbl(:,end));
        sse_set(1,tt,oo) = mdl.SSE;
        for ii = 1: size(tbl,2)-1
            a = alist; a(ii)=[];
            mdl = fitlm(tbl(:,a),tbl(:,end));
            sse_set(1+ii,tt,oo) = mdl.SSE;
        end        
    end
    disp(oo)
end
save ('CPD_period_with rw-101.mat','sse_set')
%% Coefficient of patial determinant - Rwprob
clear sse_set
warning off
alist = 1:8;

for oo = 1:size(C_raw_reg,1)
    regression_movmean = movmean(reshape(C_raw_reg(oo, set_regressionframe),[],13*fr), 0.5*fr,2); %movmean for 500ms
    craw_reg = regression_movmean(:,0.5*fr:12.5*fr-1);
    
    for tt = 1:120
        tbl = [outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt), cylinderspeed(:,tt),...
             lick_reg(:,tt),mean(craw_reg(:,3*tt-2:3*tt),2)];
        mdl = fitlm(tbl(:,1:end-1),tbl(:,end));
        sse_set(1,tt,oo) = mdl.SSE;
        for ii = 1: size(tbl,2)-1
            a = alist; a(ii)=[];
            mdl = fitlm(tbl(:,a),tbl(:,end));
            sse_set(1+ii,tt,oo) = mdl.SSE;
        end        
    end
    disp(oo)
end
% save ('CPD.mat','sse_set')
%% Draw CPD - Rwprob
day = datetime('today');
d = datestr(day,'yymmdd');
cpd = (sse_set(2:end,:,:)-sse_set(1,:,:))./sse_set(2:end,:,:);
cpd = permute(cpd,[3 2 1]);
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 17 6]);
hxa(1) = subplot(1,2,1);
xlim(hxa, [1 120]);
ylim(hxa, [0 0.6]);
hold on
stdshade(cpd(:,:,5), 0.3, [80 80 80]./255)%Lick
% stdshade(cpd(:,:,6), 0.3, [160 82 45]./255)%speed
stdshade(cpd(:,:,1), 0.3, [34 139 34]./255)% o
stdshade(cpd(:,:,3), 0.3, [199 21 133]./255)% v
lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off

hxa(2) = subplot(1,2,2);
hold on
stdshade(cpd(:,:,2), 0.3, [34 139 34]./255)%o t-1
stdshade(cpd(:,:,4), 0.3, [199 21 133]./255)%v t-1
xlim(hxa, [1 120]);
ylim(hxa, [0 0.06]); lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
hold off


xticks(hxa,[12 27 47 67 100]); xticklabels(hxa, {'C', 'D1','D2','O','ITI'})
% saveas(f1,[cd,'\figure\','CPD_withv',d,'.tif'])
%% Reversal setting
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,0.5*fr+1:3.5*fr,:),2); %cue +delay1
C_c = mean(C_event(:,0.5*fr+1:2*fr,:),2); %cue +delay1
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %cue +delay1
if ismember(trial_kind,[0 1 3 4]); %rwpn
    for itrial = rev_index(1,5):nTrial-100 % find after trial
        [p,~,stats] = anovan(lickNum(itrial:itrial+99),odorCue(itrial:itrial+99),'display','off');
        if p<0.05;
            c = multcompare(stats,'display','off');
            rwcue = find(trial_kind==4)-1;
            rwcueindex = find(strcmp(stats.grpnames{1,1},num2str(rwcue)));
            highindex = ismember(c(:,1),rwcueindex) | ismember(c(:,2),rwcueindex);
            if c(highindex,6)<0.05 & [c(ismember(c(:,1),rwcueindex),3)>0, c(ismember(c(:,2),rwcueindex),3)<0];
                aftertrial = itrial;
                break;
            end
        end
    end

elseif ismember(trial_kind, [0 1 2 3]); %rwprob
    for itrial = rev_index(1,5):nTrial-100 % find after trial
        [p,~,stats] = anovan(lickNum(itrial:itrial+99),odorCue(itrial:itrial+99),'display','off');
        if p<0.05;
            c = multcompare(stats,'display','off');
            rw80cue = find(trial_kind==2)-1;
            rw50cue = find(trial_kind==1)-1;
            rw20cue = find(trial_kind==3)-1;
            rw80cueindex = find(strcmp(stats.grpnames{1,1},num2str(rw80cue)));
            rw50cueindex = find(strcmp(stats.grpnames{1,1},num2str(rw50cue)));
            rw20cueindex = find(strcmp(stats.grpnames{1,1},num2str(rw20cue)));
            if c(:,6)<0.05 ...
                    & [c(ismember(c(:,1),rw80cueindex),3)>0; c(ismember(c(:,2),rw80cueindex),3)<0;...
                     any([c(ismember(c(:,1),rw50cueindex)&ismember(c(:,2),rw20cueindex),3)>0, c(ismember(c(:,1),rw20cueindex)&ismember(c(:,2),rw50cueindex),3)<0])];
                aftertrial = itrial;
                break;
            end
        end
    end

end
beforetrial = rev_index(1,5)-100;
analysistrial= [beforetrial:beforetrial+99,aftertrial:aftertrial+99];
save('behavior_align.mat','beforetrial','aftertrial','-append')
%% Reversal individual decoding
for oo = 1:size(C_event,3)
    ran_before = randperm(100,80); ran_after = randperm(100,60)+100; 
    Y_svm = mean(C_event(:,0.5*fr+1:3.5*fr,oo),2);
    X_svm = odorCue;
    t = templateSVM('Standardize',true,'KernelFunction','gaussian');
    Mdl1 = fitcecoc(Y_svm(analysistrial(ran_before(1:20))),X_svm(analysistrial(ran_before(1:20))),...
        'Learners',t,'FitPosterior',true,...
    'ClassNames',[0,1,3]);
    decoding_before = predict(Mdl1, Y_svm(analysistrial(ran_before(21:80))));
    decoding_after = predict(Mdl1, Y_svm(analysistrial(ran_after(1:60))));
    decoding_set(:,oo) = [mean(decoding_before==X_svm(analysistrial(ran_before(21:80)))),...
        mean(decoding_after==X_svm(analysistrial(ran_after(1:60))))];
    disp(oo)
end

save('decoding_ind','decoding_set')
scatter(decoding_set(1,:),decoding_set(2,:))
line([1/3 1/3], [0 1],'color','k', 'linestyle', ':');
line([0 1],[1/3 1/3], 'color','k', 'linestyle', ':');
xlim([0 1]); ylim([0 1]);
%% Individual decoding -time
for oo = 1:size(C_event,3)
    ran=randperm(100,20);
    Y_svm = mean(C_event(:,0.5*fr+1:3.5*fr,oo),2);
    X_svm = odorCue;
    t = templateSVM('Standardize',true,'KernelFunction','gaussian');
    Mdl1 = fitcecoc(Y_svm(analysistrial(ran(1:20)+100)),X_svm(analysistrial(ran(1:20)+100)),...
        'Learners',t,'FitPosterior',true,...
        'ClassNames',[0,1,3]);
    for itrial = 1:nTrial
        decoding_time = predict(Mdl1, Y_svm(itrial));
        decoding_set_time(itrial,oo) = decoding_time==X_svm(itrial);
    end
    disp(oo)
end
%% Cue + Delay regression -Rwpn : multi session
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,2*fr+1:3.5*fr,:),2); %cue +delay1
C_ITI = mean(C_event(:,8*fr+1:11*fr,:),2); %cue +delay1
cylinderspeed_cd = mean(cylinderspeed(:,2*10:3.5*10),2);
cylinderspeed_ITI = mean(cylinderspeed(:,8*10:11*10),2);
for ii = 1:2
    if ii ==1;     C_reg = C_cd;  speed_reg = cylinderspeed_cd;  savename = 'cuereg.mat';
    elseif ii== 2;  C_reg = C_ITI; speed_reg = cylinderspeed_ITI; savename = 'ITIreg.mat';
    end
    for oo = 1:size(C_raw)
        %before
        tbl = table(cue1_reg(:,1), cue2_reg(:,1),rew_reg(:,1), pun_reg(:,1), ...
            cue1_t1_reg(:,1),cue2_t1_reg(:,1),rew_1_reg(:,1),pun_1_reg(:,1),...
            cue1_t2_reg(:,1),cue2_t2_reg(:,1),rew_2_reg(:,1),pun_2_reg(:,1),speed_reg,C_reg(:,:,oo),...
            'VariableNames',{'Cue1','Cue2','RW','PN','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','Speed','C_raw'});
        mdl = fitlm(tbl,'C_raw~Cue1+Cue2+RW+PN+Cue1t_1+Cue2t_1+Rwt_1+Pnt_1+Cue1t_2+Cue2t_2+Rwt_2+Pnt_2+Speed');
        beta(:,oo,1) = mdl.Coefficients.Estimate;
        pvalue(:,oo,1) = mdl.Coefficients.pValue;
        xstd(:,oo,1) = std([cue1_reg(:,1), cue2_reg(:,1),rew_reg(:,1), pun_reg(:,1), ...
            cue1_t1_reg(:,1),cue2_t1_reg(:,1),rew_1_reg(:,1),pun_1_reg(:,1),...
            cue1_t2_reg(:,1),cue2_t2_reg(:,1),rew_2_reg(:,1),pun_2_reg(:,1),speed_reg]);
        ystd(:,oo,1) = std(C_reg(:,:,oo));
        src(:,oo,1) = beta(2:end,oo,1).*xstd(:,oo,1)./ystd(1,oo,1);
        
        
    end
    save(savename,'beta','pvalue','xstd','ystd','src')
end


%% Cue + Delay regression -Rwpn : reversal
beforet = beforetrial:beforetrial+99;
aftert = aftertrial:aftertrial+99;
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,0.5*fr+1:3.5*fr,:),2); %cue +delay1
C_o = mean(C_event(:,5.5*fr+1:11*fr,:),2); %cue +delay1
C_reg = C_cd;
for oo = 1:size(C_raw)
    %before
    tbl = table(cue1_reg(beforet,1), cue2_reg(beforet,1),rew_reg(beforet,1), pun_reg(beforet,1), ...
        cue1_t1_reg(beforet,1),cue2_t1_reg(beforet,1),rew_1_reg(beforet,1),pun_1_reg(beforet,1),...
        cue1_t2_reg(beforet,1),cue2_t2_reg(beforet,1),rew_2_reg(beforet,1),pun_2_reg(beforet,1),C_reg(beforet,:,oo),...
        'VariableNames',{'Cue1','Cue2','RW','PN','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
    mdl = fitlm(tbl,'C_raw~Cue1+Cue2+RW+PN+Cue1t_1+Cue2t_1+Rwt_1+Pnt_1+Cue1t_2+Cue2t_2+Rwt_2+Pnt_2');
    beta(:,oo,1) = mdl.Coefficients.Estimate;
    pvalue(:,oo,1) = mdl.Coefficients.pValue;
    xstd(:,oo,1) = std([cue1_reg(beforet,1), cue2_reg(beforet,1),rew_reg(beforet,1), pun_reg(beforet,1), ...
        cue1_t1_reg(beforet,1),cue2_t1_reg(beforet,1),rew_1_reg(beforet,1),pun_1_reg(beforet,1),...
        cue1_t2_reg(beforet,1),cue2_t2_reg(beforet,1),rew_2_reg(beforet,1),pun_2_reg(beforet,1)]);
    ystd(:,oo,1) = std(C_reg(beforet,:,oo));
    src(:,oo,1) = beta(2:end,oo,1).*xstd(:,oo,1)./ystd(1,oo,1);
    %after
    tbl = table(cue1_reg(aftert,1), cue2_reg(aftert,1),rew_reg(aftert,1), pun_reg(aftert,1), ...
        cue1_t1_reg(aftert,1),cue2_t1_reg(aftert,1),rew_1_reg(aftert,1),pun_1_reg(aftert,1),...
        cue1_t2_reg(aftert,1),cue2_t2_reg(aftert,1),rew_2_reg(aftert,1),pun_2_reg(aftert,1),C_reg(aftert,:,oo),...
        'VariableNames',{'Cue1','Cue2','RW','PN','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
    mdl = fitlm(tbl,'C_raw~Cue1+Cue2+RW+PN+Cue1t_1+Cue2t_1+Rwt_1+Pnt_1+Cue1t_2+Cue2t_2+Rwt_2+Pnt_2');
    beta(:,oo,2) = mdl.Coefficients.Estimate;
    pvalue(:,oo,2) = mdl.Coefficients.pValue;
    xstd(:,oo,2) = std([cue1_reg(aftert,1), cue2_reg(aftert,1),rew_reg(aftert,1), pun_reg(aftert,1), ...
        cue1_t1_reg(aftert,1),cue2_t1_reg(aftert,1),rew_1_reg(aftert,1),pun_1_reg(aftert,1),...
        cue1_t2_reg(aftert,1),cue2_t2_reg(aftert,1),rew_2_reg(aftert,1),pun_2_reg(aftert,1)]);
    ystd(:,oo,2) = std(C_reg(aftert,:,oo));
    src(:,oo,2) = beta(2:end,oo,2).*xstd(:,oo,2)./ystd(1,oo,2);
    disp(oo)
    
end
save('cue_reg.mat','beta','pvalue','xstd','ystd','src')
%% beta plot
sig_before = pvalue(3,:,1)<0.05;
sig_after = pvalue(3,:,2)<0.05;
sig_non = pvalue(3,:,1)>0.05&pvalue(3,:,2)>0.05;
figure
hold on
scatter(src(2,sig_non,1),src(2,sig_non,2),10,[0.3 0.3 0.3],'o','filled')
scatter(src(2,sig_before,1),src(2,sig_before,2),20,'b','o','filled')
scatter(src(2,sig_after,1),src(2,sig_after,2),30,'r','o','linewidth',1.5)
%% beta plot
sig_before = pvalue(5,:,1)<0.05;
sig_after = pvalue(5,:,2)<0.05;
sig_non = pvalue(5,:,1)>0.05&pvalue(5,:,2)>0.05;
figure
hold on
scatter(src(4,sig_non,1),src(4,sig_non,2),10,[0.3 0.3 0.3],'o','filled')
scatter(src(4,sig_before,1),src(4,sig_before,2),20,'b','o','filled')
scatter(src(4,sig_after,1),src(4,sig_after,2),30,'r','o','linewidth',1.5)

%% Cue + Delay regression -Rwprob : Reversal
beforet = beforetrial:beforetrial+99;
aftert = aftertrial:aftertrial+99;
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,0.5*fr+1:3.5*fr,:),2); %cue +delay1
C_o = mean(C_event(:,5.5*fr+1:11*fr,:),2); %cue +delay1
C_reg = C_cd;
beforet = beforetrial:beforetrial+99;
aftert = aftertrial:aftertrial+99;
for oo = 1:size(C_raw)
    %before
    tbl = table(outcome_reg(beforet,1), outcome_1_reg(beforet,1),value_reg(beforet,1), value_1_reg(beforet,1), C_reg(beforet,:,oo),...
        'VariableNames',{'o','ot_1','v','vt_1','C_raw'});
    mdl = fitlm(tbl,'C_raw~v+vt_1+o+ot_1');
    beta(:,oo,1) = mdl.Coefficients.Estimate;
    pvalue(:,oo,1) = mdl.Coefficients.pValue;
    xstd(:,oo,1) = std([outcome_reg(beforet,1), outcome_1_reg(beforet,1),value_reg(beforet,1), value_1_reg(beforet,1)]);
    ystd(:,oo,1) = std(C_reg(beforet,:,oo));
    src(:,oo,1) = beta(2:end,oo,1).*xstd(:,oo,1)./ystd(1,oo,1);
    %after
    tbl = table(outcome_reg(aftert,1), outcome_1_reg(aftert,1),value_reg(aftert,1), value_1_reg(aftert,1),C_reg(aftert,:,oo),...
        'VariableNames',{'o','ot_1','v','vt_1','C_raw'});
    mdl = fitlm(tbl,'C_raw~v+vt_1+o+ot_1');
    beta(:,oo,2) = mdl.Coefficients.Estimate;
    pvalue(:,oo,2) = mdl.Coefficients.pValue;
    xstd(:,oo,2) = std([outcome_reg(aftert,1), outcome_1_reg(aftert,1),value_reg(aftert,1), value_1_reg(aftert,1)]);
    ystd(:,oo,2) = std(C_reg(aftert,:,oo));
    src(:,oo,2) = beta(2:end,oo,2).*xstd(:,oo,2)./ystd(1,oo,2);
    disp(oo)
    
end
save('cue_reg.mat','beta','pvalue','xstd','ystd','src')
%% beta plot
sig_before = pvalue(4,:,1)<0.05;
sig_after = pvalue(4,:,2)<0.05;
sig_non = pvalue(4,:,1)>0.05&pvalue(4,:,2)>0.05;
hold on
scatter(src(3,sig_non,1),src(3,sig_non,2),10,[0.3 0.3 0.3],'o','filled')
scatter(src(3,sig_before,1),src(3,sig_before,2),20,'b','o','filled')
scatter(src(3,sig_after,1),src(3,sig_after,2),30,'r','o','linewidth',1.5)

%% population decoding
Y_svm_cue_temp=Y_svm(1:sum(rev_b),1,1);
pn0_trial = find(Y_svm_cue_temp==0 | Y_svm_cue_temp==2);
Y_svm_cue = Y_svm_cue_temp(pn0_trial);
X_svms = X_svm(pn0_trial,:,:);

for trials = 1:size(X_svms,1)
    samplex = X_svms(trials,:,:);
    X_svm_train = X_svms;
    X_svm_train(trials,:,:)=[];
    sampley = Y_svm_cue(trials);
    Y_svm_train = Y_svm_cue;
    Y_svm_train(trials)=[];
    for tt = 1:50
        Mdl1 = fitcsvm(X_svm_train(:,:,tt),Y_svm_train); %cue
        label = predict(Mdl1,samplex(:,:,tt));
        Prob_svm(trials,tt,1) = (label== sampley);
        
    end
    disp(trials)
end

%% Decoding shuffle
%cue
Y_svm_cue_temp=Y_svm(1:sum(rev_b),1,1);
pn0_trial = find(Y_svm_cue_temp==0 | Y_svm_cue_temp==2);
Y_svm_cue = Y_svm_cue_temp(pn0_trial);
X_svms = X_svm(pn0_trial,:,:);

%100 permutation
for rant = 1:100
    %cue
    Y_svms_shu = Y_svm_cue(randperm(size(Y_svm_cue,1)),:,:);
    for trials = 1:size(X_svms,1)
        samplex = X_svms(trials,:,:);
        X_svm_train = X_svms;
        X_svm_train(trials,:,:)=[];
        sampley = Y_svms_shu(trials);
        Y_svm_train = Y_svms_shu;
        Y_svm_train(trials,:,:)=[];
        for tt = 1:50
            Mdl1 = fitcsvm(X_svm_train(:,:,tt),Y_svm_train); %cue
            label = predict(Mdl1,samplex(:,:,tt));
            Prob_svm_cue_shu(trials,tt,rant) = (label== sampley);
        end
    end
    
    disp(rant)
end


%% Draw decoding
Prob_cue_sig = mean(Prob_svm(:,:,1))>mean(Prob_svm_cue_shu);
Prob_cue_sig2 = sum(Prob_cue_sig,3)>95;

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 9]);
hold on
xlim([0 50]); xticks([10 20 30]); xticklabels({'Cue', 'Delay','Outcome'})
ylim([0.3 1])
plot(mean(Prob_svm(:,:,1),1),'color', [50 50 50]./255,'linewidth',2)
plot(Prob_cue_sig2*0.99,'color', [50 50 50]./255,'marker','*', 'linestyle','none')
plot(mean(mean(Prob_svm_cue_shu(:,:,:),3)),'color', [50 50 50]./255,'linewidth',1)

lim = axis;
line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([15 15], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([25 25], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
%line([lim(1) lim(2)], [0.5 0.5],'color','k', 'linestyle', '-.'); %2.5 outcome

saveas(f1,[figure_folder,'\','decoding','_craw.tif'])
%% GLM - Rwpn
%parameter
term = [3 5 3];
winstep  = {-2:1:35;...
    -2:1:35;...
    -2:1:35;...
    -2:1:57;...
    -2:1:57;...
    -2:1:57;...
    -2:1:57;...
    -2:1:57;...
    -6:1:13;...
    -6:1:78;...
    -6:1:78;};
nWinBin = cellfun(@length, winstep);
varName = {'Rw Cue', 'Pn Cue', 'Cue C', ...
    'Reward', 'Non-reward', 'Punishment', 'Non-punishment', 'No-outcome',...
    'Lick onset', 'Mid-bout lick','Lick offset'};
framelist = {'cuerw_frame','cuepn_frame','cue0_frame','o_rw','o_rwx','o_pn','o_pnx','o_no','lick_frame'};
xlist = {'cuerw_frame','cuepn_frame','cue0_frame','o_rw','o_rwx','o_pn','o_pnx','o_no','lickonset','lickmid','lickoffset'};

%z-scored f/f
[X3,Y3] = meshgrid(0:3:11*fr-1,GPIO1onframe(1:sum(rev_b))); %0 ~ 11 s
set_glmframe = X3+Y3;
% normalization
C_mean = mean(C_raw(:,GPIO1onframe(1):GPIO1onframe(end)),2);
C_std = std(C_raw(:,GPIO1onframe(1):GPIO1onframe(end)),0,2);
C_n = (C_raw-C_mean)./C_std;
%get frame
C_glmf = movmean(C_n,[7 7],2); %filter
C_glm100 = movmean(C_glmf,[1 1],2); %100ms
C_glm100 = C_glm100(:,2:3:end,:);
%trial type
odorCue_revb = odorCue(1:sum(rev_b));
ind_R = odorCue_revb==(find(trial_kind==1)-1);
ind_P = odorCue_revb==(find(trial_kind==4)-1);
ind_0 = odorCue_revb==(find(trial_kind==3)-1);
%period
cue_period = [0.5]; outcome_period = [5.5];
%cue frame
cue_frame = zeros(sum(rev_b),11*10); cue_frame(:,cue_period(1)*10+1)=1; %only onset==1
cuerw_frame = ceil(set_glmframe.*cue_frame.*ind_R./3); 
cuepn_frame = ceil(set_glmframe.*cue_frame.*ind_P./3);
cue0_frame=ceil(set_glmframe.*cue_frame.*ind_0./3);
%outcomeframe
outcome_revb = waterReward(1:sum(rev_b));
outcome_frame = zeros(sum(rev_b),11*10);outcome_frame(:,outcome_period*10+1)=1;
o_rw = ceil(set_glmframe.*outcome_frame.*ind_R.*outcome_revb./3); 
o_rwx = ceil(set_glmframe.*outcome_frame.*ind_R.*(1-outcome_revb)./3);
o_pn = ceil(set_glmframe.*outcome_frame.*ind_P.*outcome_revb./3); 
o_pnx = ceil(set_glmframe.*outcome_frame.*ind_P.*(1-outcome_revb)./3);
o_no = ceil(set_glmframe.*outcome_frame.*ind_0./3);
%lick frame
lickinfo = [timeline{1,1}];%timeline{1,2}
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_frame_temp(iii,:) = histcounts(licktemp, 0:100000:11000000);
end
lick_frame_temp(lick_frame_temp>1)=1;
lick_frame = ceil(set_glmframe.*lick_frame_temp(1:sum(rev_b),:)./3);

%% GLM - Rwprob
%parameter
winstep  = {-2:1:35;...
    -2:1:35;...
%     -2:1:35;...
    -2:1:57;...
    -2:1:57;...
    -6:1:13;...
    -6:1:78;...
    -6:1:78;};
nWinBin = cellfun(@length, winstep);
varName = {'Value',...%'Cue75', 'Cue25', 'Cue0', ...
    'Reward', 'Non-reward', ...
    'Lick onset', 'Mid-bout lick','Lick offset'};
framelist = {'cue75_frame','cue25_frame','cue0_frame','cue_frame','o_rw','o_rwx','lick_frame'};
% xlist = {'cue75_frame','cue25_frame','cue0_frame','o_rw','o_rwx','lickonset','lickmid','lickoffset'};
% framelist = {'Value','o_rw','o_rwx','lick_frame'};
xlist = {'cue_frame','Value','o_rw','o_rwx','lickonset','lickmid','lickoffset'};

%z-scored f/f
[X3,Y3] = meshgrid(0:3:11*fr-1,GPIO1onframe(1:sum(rev_b))); %0 ~ 11 s
set_glmframe = X3+Y3;
% normalization
C_mean = mean(C_raw(:,GPIO1onframe(1):GPIO1onframe(end)),2);
C_std = std(C_raw(:,GPIO1onframe(1):GPIO1onframe(end)),0,2);
C_n = (C_raw-C_mean)./C_std;
%get frame
C_glmf = movmean(C_n,[7 7],2); %filter
C_glm100 = movmean(C_glmf,[1 1],2); %100ms
C_glm100 = C_glm100(:,2:3:end,:);
%trial type
odorCue_revb = odorCue(1:sum(rev_b));
ind_75 = odorCue_revb==(find(trial_kind==1)-1);
ind_25 = odorCue_revb==(find(trial_kind==2)-1);
ind_0 = odorCue_revb==(find(trial_kind==3)-1);
%period
cue_period = [0.5]; outcome_period = [5.5];
%cue frame
cue_frame = zeros(sum(rev_b),11*10); cue_frame(:,cue_period(1)*10+1)=1; %only onset==1
cue75_frame = ceil(set_glmframe.*cue_frame.*ind_75./3); 
cue25_frame = ceil(set_glmframe.*cue_frame.*ind_25./3);
cue0_frame=ceil(set_glmframe.*cue_frame.*ind_0./3);
%outcomeframe
outcome_revb = waterReward(1:sum(rev_b));
outcome_frame = zeros(sum(rev_b),11*10);outcome_frame(:,outcome_period*10+1)=1;
o_rw = ceil(set_glmframe.*outcome_frame.*outcome_revb./3); 
o_rwx = ceil(set_glmframe.*outcome_frame.*(1-outcome_revb)./3);

%lick frame
lickinfo = [timeline{1,1}];%timeline{1,2}
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_frame_temp(iii,:) = histcounts(licktemp, 0:100000:11000000);
end
lick_frame_temp(lick_frame_temp>1)=1;
lick_frame = ceil(set_glmframe.*lick_frame_temp(1:sum(rev_b),:)./3);

%% make linear & separate lick term
for ii = 1:length(framelist)
    frametmp = zeros(1,size(C_glm100,2));
    eval(['frametmp(',framelist{ii},'(',framelist{ii},'~=0)) = 1;'])
    eval([framelist{ii},'=','frametmp;']);
end
if exist('cue75_frame')
Value = cue75_frame*0.75+cue25_frame*0.25+cue0_frame*0;
end
%find lick onset, offset
licksum = movsum(lick_frame,[19 0], 2);
lickonset = (licksum==1)& ([0 licksum(1:end-1)]==0);
licksum = movsum(lick_frame,[0 19], 2);
lickoffset = (licksum==1)& ([licksum(2:end) 0]==0);
lickmid = lick_frame;
framelength = size(C_glm100,2);
%%
Xbintmp = zeros(framelength,sum(nWinBin));
for ivar = 1:length(winstep)
    eval(['xtmp=',xlist{ivar},';']);
    for istep = 1:length(winstep{ivar,1})
        nmove = winstep{ivar,1}(1,istep);
        inrange = [1:length(xtmp)]-nmove ;
        inrange(inrange<=0)=[];inrange(inrange>length(xtmp))=[];
        Xbintmp(nmove*(nmove>0)+1:end+nmove*(nmove<0),sum(nWinBin(1:ivar-1))+istep) = xtmp(inrange);
    end
end

Xbin = Xbintmp(ceil(min(set_glmframe(:))./3):ceil(max(set_glmframe(:))./3),:);
C_glmset = C_glm100(:,ceil(min(set_glmframe(:))./3):ceil(max(set_glmframe(:))./3));
%% GLM fit
glmresult=cell(1,size(C_raw,1));
for oo = 1:size(C_raw,1)
    C_event = C_glmset(oo,:);
    [b, dev, stats] = glmfit(Xbin, C_event);
    glmresult{oo}=stats;
    disp(oo);
end
save('glm.mat','glmresult','xlist','winstep');

%% draw all glm coefficient
clear glmbetaset betamean
for oo = 1:size(glmresult,2)
    glmbetaset(:,oo) = glmresult{1,oo}.beta;
    glmtvalueset(:,oo) = glmresult{1,oo}.t; 
end

betamean = mean(glmbetaset(2:end,:),2);
tmean = mean(glmtvalueset(2:end,:),2);
%%
nst=1;
for ii = 1:size(winstep,1)
subplot(6,2,ii)
plot(betamean(nst:nst+size(winstep{ii,1},2)-1,1))

nst = nst+size(winstep{ii,1},2);
end
figure;
nst=1;
for ii = 1:size(winstep,1)
subplot(6,2,ii)
plot(tmean(nst:nst+size(winstep{ii,1},2)-1,1))
nst = nst+size(winstep{ii,1},2);
end
%% draw glm t value-rwpn
nst=1;
f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_outcome = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);

f_licks = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_lickm = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_licke = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
for ii = 1:size(winstep,1)
    if contains(xlist{1,ii},'cue')
        cset = [0 0 205; 255 30 70; 80 80 80]./255;
        colortmp = cset(contains(xlist{1,ii},'rw')+2*contains(xlist{1,ii},'pn')+3*contains(xlist{1,ii},'0'),:);
        figure(f_cue); hold on;

    elseif contains(xlist{1,ii},'o_rw')
        cset = [0 0 205; 129 212 250]./255;
        colortmp = cset(contains(xlist{1,ii},'x')+1,:);
        figure(f_outcome);hold on;
    elseif contains(xlist{1,ii},'o_pn')
        cset = [255 30 70; 239 154 154]./255;
        colortmp = cset(contains(xlist{1,ii},'x')+1,:);
        figure(f_outcome);hold on;
    elseif contains(xlist{1,ii},'o_no')
        colortmp = [80 80 80]./255;
        figure(f_outcome);hold on;
    elseif contains(xlist{1,ii},'lick')
        colortmp = [0 0 0]./255;
        if contains(xlist{1,ii},'onset'); figure(f_licks);
        elseif contains(xlist{1,ii},'mid'); figure(f_lickm);
        elseif contains(xlist{1,ii},'offset'); figure(f_licke);
        end
        hold on;
    end
    stdshade(glmbetaset(nst:nst+size(winstep{ii,1},2)-1,:)', 0.3, colortmp,winstep{ii,1})%pncue
    xlim([winstep{ii,1}(1,1) winstep{ii,1}(1,end)])
    nst = nst+size(winstep{ii,1},2);
end
figure(f_cue); title('Cue'); ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_cue,[figure_folder,'\glm_cue.tif'])
    
    figure(f_outcome); title('Outcome');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([25 25],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_outcome,[figure_folder,'\glm_outcome.tif'])
    
    figure(f_licks); title('Lick onset');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_licks,[figure_folder,'\glm_licks.tif'])
    figure(f_lickm); title('Lick');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_lickm,[figure_folder,'\glm_lickm.tif'])
    figure(f_licke); title('Lick offset');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_licke,[figure_folder,'\glm_licke.tif'])
    close all
%% draw glm t value - rwprob
nst=1;
ylimset = [0 0.25; 0 0.5; 0 0.2; 0 0.1; 0 0.35;... %beta
    0 3.5;  0 5;   0 3;   0 2.5; 0 5];
f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_value = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_outcome = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_licks = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_lickm = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
f_licke = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
for ii = 1:size(winstep,1)
    if contains(xlist{1,ii},'cue')
        colortmp = [0 0 205]./255;
%         cset = [0 0 205; 0 128 255; 80 80 80]./255;
%         colortmp = cset(contains(xlist{1,ii},'75')+2*contains(xlist{1,ii},'25')+3*contains(xlist{1,ii},'0'),:);
        figure(f_cue); hold on; xticks([0 15 30 35]);
    elseif contains(xlist{1,ii},'Value')
        colortmp = [199 21 133]./255;
        figure(f_value); hold on;
    elseif contains(xlist{1,ii},'o_rw')
        cset = [0 0 205; 80 80 80]./255;
        colortmp = cset(contains(xlist{1,ii},'x')+1,:);
        figure(f_outcome);hold on; xticks([0 25]);
    elseif contains(xlist{1,ii},'lick')
        colortmp = [0 0 0]./255;
        if contains(xlist{1,ii},'onset'); figure(f_licks);
        elseif contains(xlist{1,ii},'mid'); figure(f_lickm);
        elseif contains(xlist{1,ii},'offset'); figure(f_licke);
        end
        hold on;
    end
    stdshade(glmbetaset(nst:nst+size(winstep{ii,1},2)-1,:)', 0.3, colortmp,winstep{ii,1})%pncue
    xlim([winstep{ii,1}(1,1) winstep{ii,1}(1,end)])
    nst = nst+size(winstep{ii,1},2);
    
end
figure(f_cue); title('Cue'); ax = axis;
line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
saveas(f_cue,[figure_folder,'\glm_cue.tif'])
figure(f_value); title('Value'); ax = axis;
line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
saveas(f_value,[figure_folder,'\glm_value.tif'])
figure(f_outcome); title('Reward');ax = axis;
line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([25 25],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
saveas(f_outcome,[figure_folder,'\glm_rw.tif'])
figure(f_licks); title('Lick onset');ax = axis;
line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
saveas(f_licks,[figure_folder,'\glm_licks.tif'])
figure(f_lickm); title('Lick');ax = axis;
line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
saveas(f_lickm,[figure_folder,'\glm_lickm.tif'])
figure(f_licke); title('Lick offset');ax = axis;
line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
saveas(f_licke,[figure_folder,'\glm_licke.tif'])
close all


%% time correlation
[X,Y] = meshgrid(0:12*fr-1,GPIO1onframe); % 0~12
set_eventframe = X+Y; %meshgrid for event frame
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
%C_event = reshape(permute(C_corr2,[2 1 3]), [],1,size(C_raw,1));
C_corr = movmean(C_z,[0 29],2); %window 1s
C_corr2 = permute(C_corr(:,1:6:end,:),[3 2 1]); %moving 200ms
clear Rset Pset
for ii = 1:nTrial
    [R,P] = corrcoef(C_corr2(:,:,ii));
    Rset(:,:,ii) = R;
    Pset(:,:,ii) = P;
end
clear Rset_rwo Pset_rwo
[R,P] = corrcoef(C_corr2(:,:,1:10));

for ii = 1:nTrial
    
    [R,P] = corrcoef(C_corr2(unique(cro),:,ii));
    Rset_rwo(:,:,ii) = R;
    Pset_rwo(:,:,ii) = P;
end
colormap(hot)
imagesc(mean(Rset,3))
%imagesc(mean(Rset(:,:,trial_80Rwo_ind),3))
linelist = [0.5*fr/3 0.5*fr/3 0 120;...
            2*fr/3     2*fr/3 0 120;...
            3.5*fr/3 3.5*fr/3 0 120;...
            4*fr/3     4*fr/3 0 120;...
            5.5*fr/3 5.5*fr/3 0 120;...
            8*fr/3     8*fr/3 0 120];
for Li = 1:size(linelist,1)*2
    line([linelist(Li,1) linelist(Li,2)],[linelist(Li,3) linelist(Li,4)],'color','w', 'linestyle', ':','linewidth',1.5); %2 cue
    line([linelist(Li,3) linelist(Li,4)],[linelist(Li,1) linelist(Li,2)],'color','w', 'linestyle', ':','linewidth',1.5); %2 cue
end
%% active time
C_raw_filter = movmean(C_raw, [7 7],2); % movmean window 500ms
[X2,Y2] = meshgrid(-0.5*fr:12.5*fr-1,GPIO1onframe); %-0.5 ~ 12.5 ; 6~105 frame
set_regressionframe = X2+Y2;
test_region= [0.5 2 3.5 4 5.5 8 11]+0.5;
for oo = 1:size(C_raw,1)
    C_event = reshape(C_raw_filter(oo,set_regressionframe(:)),[],13*fr);
    %     eval(['zind=', zlist{fix((ii+1)/2)},';']);
    %     Ctemp = mean(C_event(zind,1:0.5*fr),2);
    %     CM = mean(Ctemp); CSTD=std(Ctemp);
    for rr = 1:length(test_region)-1
        tempcomp = mean(C_event(:,1+fr*(test_region(rr)-1):fr*test_region(rr)),2);
        testcomp = movmean(C_event(:,1+fr*test_region(rr):fr*test_region(rr+1)),[0 5],2,'Endpoints','discard');
        window100 = 1:3:size(testcomp,2);
        for ww=1:length(window100);
            [p,h]=signrank(tempcomp,testcomp(:,ww));
            pset(ww) = p;
        end
        psetco= mafdr(pset,'BHFDR',true);
        pset_noco{rr,oo}= pset';
        pset_co{rr,oo} = psetco';
        pset =[];
    end
end



pset_noco_s = cellfun(@(x) sum(x<0.001), pset_noco);
pset_co_s = cellfun(@(x) sum(x<0.001), pset_co);

sig_cue = sum(pset_co_s(1,:)~=0)./size(pset_co_s,2);
sig_d1 = sum(pset_co_s(2,:)~=0)./size(pset_co_s,2);
sig_m = sum(pset_co_s(3,:)~=0)./size(pset_co_s,2);
sig_d2 = sum(pset_co_s(4,:)~=0)./size(pset_co_s,2);
sig_o = sum(pset_co_s(5,:)~=0)./size(pset_co_s,2);
sig_iti = sum(pset_co_s(6,:)~=0)./size(pset_co_s,2);

%% Co-activity
% [X,Y] = meshgrid(0:12*fr-1,GPIO1onframe); % 0~12
% set_eventframe = X+Y; %meshgrid for event frame
% S_event = S_raw(:,set_eventframe(:));
% S_event = reshape(S_event', [],12*fr,size(S_event,1));
% S_event01 = S_event>0;
mouseid = {'dHP03','dHP04','dHP06','dHP07','dHP08'...
    'vHP06','vHP07','vHP08','vHP11'};
rwpng = {'rwpn01','rwpn02','rwpn03',...
    'rwpnrev01','rwpnrev02','rwpnrev03','rwpnrev04','rwpnrev05',...
    'rwpnrev06','rwpnrev07','rwpnrev08','rwpnrev09','rwpnrev10',...
    'rwpnrev11','rwpnrev12','rwpnrev13'};
CoaNum10_300 = NaN(9,16);
for imouse = 1:9
    for ii = 5
        try
            cd(['E:\data\vHPC\',mouseid{1,imouse},'\',rwpng{1,ii}])
            load('behavior_align.mat','fr','GPIO1onframe','S_raw','set_eventframe')
            S_all = S_raw(:,set_eventframe(1,1):set_eventframe(size(set_eventframe,1),size(set_eventframe,2))+30)>0;
            coaThr = 0.1*size(S_all,1); %coactivity threshold
            clear S_coa
            %500ms:15 300ms:9
            for iwin = 1:size(S_all,2)-9
                S_coa(iwin)=sum(sum(S_all(:,iwin:iwin+8),2)>0,1);
            end
            CoaNum10_300(imouse,ii) = sum(diff(S_coa>coaThr)==1);
            S_coaset{imouse,ii} = S_coa;
%             if strcmp(rwpng{1,ii},'rwpnrev01');
%                 load('behavior_align.mat','rev_index')
                rev_index=[0 0 0 0 142];
                set_eventframe = set_eventframe-set_eventframe(1,1)+1;
                CoaNum_rev2(imouse,:) = [sum(diff(S_coa(1,set_eventframe(rev_index(5)-100,1):set_eventframe(rev_index(5)-51,360))>coaThr)==1),...
                    sum(diff(S_coa(1,set_eventframe(rev_index(5)-50,1):set_eventframe(rev_index(5)-1,360))>coaThr)==1),...
                    sum(diff(S_coa(1,set_eventframe(rev_index(5),1):set_eventframe(rev_index(5)+49,360))>coaThr)==1),...
                    sum(diff(S_coa(1,set_eventframe(rev_index(5)+50,1):set_eventframe(rev_index(5)+99,360))>coaThr)==1),...
                    sum(diff(S_coa(1,set_eventframe(rev_index(5)+100,1):set_eventframe(rev_index(5)+149,360))>coaThr)==1),...
                    sum(diff(S_coa(1,set_eventframe(rev_index(5)+150,1):set_eventframe(rev_index(5)+199,360))>coaThr)==1),...
                    sum(diff(S_coa(1,set_eventframe(rev_index(5)+200,1):set_eventframe(rev_index(5)+249,360))>coaThr)==1),];
            
%             end
                
        end
    end
end

%% All coactivity
cmap = [128 0 128; 255 0 255;  238 130 238; 153 50 204; 218 112 214; ...
    34 139 34; 50 205 50; 144 238 144; 0 128 0]./255; 
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
for ii = 1:9
    if rev(1,ii)==0; lines= '-o'; elseif rev(1,ii)==1; lines= '-+'; end
    hold on
    subplot(2,1,dvHP(1,ii)+1)
    plot(CoaNum10_300(ii,:), lines,'Color',cmap(ii,:))
end
subplot(2,1,1)
lim = axis;
line([3.5 3.5],[lim(3) lim(4)],'color','k', 'linestyle', ':')
xticks([1 2 3 4 8 13 ]); xlim([0 16]); ylim([0 1000])
xticklabels({'L1','L2','L3','Rev1','Rev5','Rev10'})
subplot(2,1,2)
lim = axis;
line([3.5 3.5],[lim(3) lim(4)],'color','k', 'linestyle', ':')
xticks([1 2 3 4 8 13 ]); xlim([0 16]);ylim([0 1000])
xticklabels({'L1','L2','L3','Rev1','Rev5','Rev10'})
saveas(f1,['E:\data\vHPC\all\CoaNum10_300.tif'])
%%
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 15]);

dvHP = ~logical(dvHP);rev = logical(rev);
subplot(2,1,1)
hold on
errorbar(1:7, mean(CoaNum_l1(dvHP,:),1),sem(CoaNum_l1(dvHP,:)),'Color',[0.9 0 0.9]);
errorbar(9:15,mean(CoaNum_l2(dvHP,:),1),sem(CoaNum_l2(dvHP,:)),'Color',[0.8 0 0.8]);
errorbar(17:23,mean(CoaNum_l3(dvHP,:),1),sem(CoaNum_l3(dvHP,:)),'Color',[0.7 0 0.7]);
errorbar(25:31,mean(CoaNum_rev1(dvHP,:),1),sem(CoaNum_rev1(dvHP,:)),'Color',[0.3 0 0.3]);
errorbar(1:7,mean(CoaNum_l1(~dvHP,:),1),sem(CoaNum_l1(~dvHP,:)),'Color',[0 0.9 0]);
errorbar(9:15,mean(CoaNum_l2(~dvHP,:),1),sem(CoaNum_l2(~dvHP,:)),'Color',[0 0.8 0]);
errorbar(17:23,mean(CoaNum_l3(~dvHP,:),1),sem(CoaNum_l3(~dvHP,:)),'Color',[0 0.7 0]);
errorbar(25:31,mean(CoaNum_rev1(~dvHP,:),1),sem(CoaNum_rev1(~dvHP,:)),'Color',[0 0.3 0]);

subplot(2,1,2)
hold on
errorbar(1:7, mean(CoaNum_l1(dvHP&~rev,:),1),sem(CoaNum_l1(dvHP&~rev,:)),'Color',[0.9 0 0.9]);
errorbar(9:15,mean(CoaNum_l2(dvHP&~rev,:),1),sem(CoaNum_l2(dvHP&~rev,:)),'Color',[0.8 0 0.8]);
errorbar(17:23,mean(CoaNum_l3(dvHP&~rev,:),1),sem(CoaNum_l3(dvHP&~rev,:)),'Color',[0.7 0 0.7]);
errorbar(25:31,mean(CoaNum_rev1(dvHP&~rev,:),1),sem(CoaNum_rev1(dvHP&~rev,:)),'Color',[0.3 0 0.3]);
errorbar(1:7,mean(CoaNum_l1(~dvHP&~rev,:),1),sem(CoaNum_l1(~dvHP&~rev,:)),'Color',[0 0.9 0]);
errorbar(9:15,mean(CoaNum_l2(~dvHP&~rev,:),1),sem(CoaNum_l2(~dvHP&~rev,:)),'Color',[0 0.8 0]);
errorbar(17:23,mean(CoaNum_l3(~dvHP&~rev,:),1),sem(CoaNum_l3(~dvHP&~rev,:)),'Color',[0 0.7 0]);
errorbar(25:31,mean(CoaNum_rev1(~dvHP&~rev,:),1),sem(CoaNum_rev1(~dvHP&~rev,:)),'Color',[0 0.3 0]);
saveas(f1,['E:\data\vHPC\all\CoaNum_trial.tif'])

%%
figure
mkdir([figure_folder, '\coa\'])
for ii = 1:size(S_event01,1)
imagesc(permute(S_event01(ii,:,:),[3 2 1]))
colormap([1 1 1; 0 0 0])
lim = axis;
line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
xlim([0 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'C', 'D1','D2','O'})
hold off
saveas(gcf,[figure_folder, '\coa\',num2str(ii),'.tif']);
end
%%
S_all = S_raw(:,set_eventframe(1,1):set_eventframe(size(set_eventframe,1),size(set_eventframe,2))+30)>0;
coaThr = 0.05*size(S_all,1); %coactivity threshold
clear S_coa
%500ms:15 300ms:9
for iwin = 1:size(S_all,2)-9
    S_coa(iwin)=sum(sum(S_all(:,iwin:iwin+8),2)>0,1);
end
CoaNum(imouse,ii) = sum(diff(S_coa>size(S_all,1)*0.05)==1);
[cotrial, cotime] = find(diff(coa,1,1)==1);



%% ITI coactivated cell

S_ITI=S_event01(:,8*fr+1:11*fr,:);
clear S_co ITI_cell
for ii =1:size(S_ITI,2)-15;
    S_co(:,ii,:)=sum(S_ITI(:,ii:ii+14,:),2)>0;
end
[S_coind,S_cotime] = max(sum(S_co,3),[],2);
S_coa = S_event01(S_coind>=coaThr,:,:);
for ii = 1:size(S_event01,1)
    ITI_cell(ii,:,:) = sum(S_ITI(ii,S_cotime(ii):S_cotime(ii)+14,:),2)>0; % max time fire cell numbering
end

ITI_trial.rwo = find((S_coind>=coaThr) & trial_80Rwo_ind); [rro, cro] = find(ITI_cell(ITI_trial.rwo,:)>0); %colume = cell number
ITI_trial.rwx = find((S_coind>=coaThr) & trial_80Rwx_ind); [rrx, crx] = find(ITI_cell(ITI_trial.rwx,:)>0);
ITI_trial.pno = find((S_coind>=coaThr) & trial_80Pno_ind); [rpo, cpo] = find(ITI_cell(ITI_trial.pno,:)>0);
ITI_trial.pnx = find((S_coind>=coaThr) & trial_80Pnx_ind); [rpx, cpx] = find(ITI_cell(ITI_trial.pnx,:)>0);
ITI_trial.nto = find((S_coind>=coaThr) & trial_20R_ind);   [roo, coo] = find(ITI_cell(ITI_trial.nto,:)>0);


itilist = 1*((S_coind>=coaThr) & trial_80Rwo_ind)+2*((S_coind>=coaThr) & trial_80Rwx_ind)+...
    3*((S_coind>=coaThr) & trial_80Pno_ind)+4*((S_coind>=coaThr) & trial_80Pnx_ind)+...
    5*((S_coind>=coaThr) & trial_20R_ind);
itilist(itilist==0)=[];


mean(sum(ITI_cell,2))/size(ITI_cell,2)
size(S_coa,1)
mean(sum(ITI_cell(S_coind>=coaThr,:),2))./size(S_coa,3)
[sum(itilist==1), sum(itilist==2),sum(itilist==3), sum(itilist==4), sum(itilist==5)]
[sum(itilist==1)/sum(trial_80Rwo_ind),sum(itilist==2)/sum(trial_80Rwx_ind),...
    sum(itilist==3)/sum(trial_80Pno_ind), sum(itilist==4)/sum(trial_80Pnx_ind),...
    sum(itilist==5)/sum(trial_20R_ind)]

%% total number of S
S_event = S_raw(:,set_eventframe(:));
S_event = reshape(S_event', [],12*fr,size(S_event,1));
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
s = mean(mean(S_event>0,1),3);
hold on
plot(s);
lim = axis;
line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
xlim([0 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'C', 'D1','D2','O'})
hold off
saveas(f1,[figure_folder, '\totalS.tif']);

%total number of S according to trial
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
triallist= {'trial_20R_ind','trial_80Rwo_ind','trial_80Pno_ind','trial_80Rwx_ind','trial_80Pnx_ind'};
colorset = [ 100 100 100; 30 45 232; 255 30 70; 120 177 255; 232 126 58]/255; %blue red skyblue pink grey
for ii=1:5
    eval(['trialind=',triallist{ii},';',]);
    s = mean(mean(S_event01(trialind&rev_b',:,:),1),3);
    hold on
    plot(s,'color',colorset(ii,:),'linewidth',2)
end
lim = axis;
line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
xlim([0 12*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'C', 'D1','D2','O'})
hold off
saveas(f1,[figure_folder, '\totalS_trial_brev.tif']);

%histogram mean activated cell
task_region= [0.5 2 3.5 4 5.5 8 12];
triallist= {'trial_80Rwo_ind','trial_80Pno_ind','trial_20R_ind','trial_80Rwx_ind','trial_80Pnx_ind'};
titlelist = {'Rw', 'Pn', 'NT', 'No Rw','NoPn','Mean activated cells'};
S_ITI=sum(S_event(:,8*fr+1:11.5*fr,:),2)>0;
S_C=sum(S_event(:,0.5*fr+1:2*fr,:),2)>0;
S_D=sum(S_event(:,4*fr+1:5.5*fr,:),2)>0;
S_O=sum(S_event(:,5.5*fr+1:8*fr,:),2)>0;

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 25 10]);
S_g = S_ITI;
S_80ro = sum(S_g(trial_80Rwo_ind,:,:),3); S_80rx=sum(S_g(trial_80Rwx_ind,:,:),3);
S_80po = sum(S_g(trial_80Pno_ind,:,:),3); S_80px=sum(S_g(trial_80Pnx_ind,:,:),3);
S_00xx = sum(S_g(trial_20R_ind,:,:),3);
S_xg = [ones(1,length(S_80ro)),2*ones(1,length(S_80po)),3*ones(1,length(S_00xx)),4*ones(1,length(S_80rx)),5*ones(1,length(S_80px))];
[p,tbl,stats]=anovan([S_80ro;S_80po;S_00xx;S_80rx;S_80px],{S_xg});
[c,m,h,nms] = multcompare(stats);

for ii = 1:6
    subplot(2,3,ii)
    if ii == 6;
        hold on
        bar([1,2,3,4,5], [mean(S_80ro),mean(S_80po),mean(S_00xx),mean(S_80rx),mean(S_80px)])
        errorbar([mean(S_80ro),mean(S_80po),mean(S_00xx),mean(S_80rx),mean(S_80px)],...
            [sem(S_80ro),sem(S_80po),sem(S_00xx),sem(S_80rx),sem(S_80px)],'linestyle','none')
        hold off
        xlim([0 6]); xticks(1:5); xticklabels({'Rw','No Rw','Pn','No Pn','NT'});
    else
        eval(['trialind=',triallist{ii},';',]);
        a= sum(S_g(trialind,:,:));
        histogram(a,0:3:60);
        xlim([0 60])
    end
    title(titlelist{ii})
end
 saveas(f1,[figure_folder,'\act d2.tif']);
%%
clear S_shuffle_coa
sh_N =20;
S_shuffle_coa= zeros([size(S_event01),sh_N]);
for oo = 1:size(S_event,3)
    for ii = 1:sh_N
        S_shuffle_coa(:,:,oo,ii) = S_event01(:,randperm(360),oo);
    end
    disp(oo)
end
S_shu_cri = mean(S_shuffle_coa,4)+3*std(S_shuffle_coa,0,4);
S_shu_cri= movmean(sum(S_shu_cri,3),[0 15],2,'endpoints','discard');
S_eve_cri = movmean(sum(S_event01,3),[0 15],2,'endpoints','discard');

sum(S_eve_cri(:)>S_shu_cri(:));
%% Anova
task_region= [0.5 2 3.5 4 5.5 8 11];
[X,Y] = meshgrid(0:12*fr-1,GPIO1onframe); % 0~12
set_eventframe = X+Y; %meshgrid for event frame
S_event = S_raw(:,set_eventframe(:));
S_event = reshape(S_event', [],12*fr,size(S_event,1));
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
S_event01 = S_event>0;
S_B=sum(S_event(:,0*fr+1:0.5*fr,:),2);
S_C=sum(S_event(:,0.5*fr+1:2*fr,:),2);
S_D1=sum(S_event(:,2*fr+1:3.5*fr,:),2);
S_D2=sum(S_event(:,4*fr+1:5.5*fr,:),2);
S_O=sum(S_event(:,5.5*fr+1:8*fr,:),2);
S_ITI_a=sum(S_event(:,8*fr+1:11*fr,:),2);
S_set = [S_C,S_D1,S_D2, S_O,S_ITI_a];
clear p_set c_set
c_set = cell(size(S_ITI_a,3),5);
pset = ones(size(S_ITI_a,3),5);
for jj = 1:5
    S_g= S_set(:,jj,:);
    for ii = 1:size(S_ITI_a,3)
        if jj>=4
        S_80ro = sum(S_g(trial_80Rwo_ind,:,ii),3); S_80rx=sum(S_g(trial_80Rwx_ind,:,ii),3);
        S_80po = sum(S_g(trial_80Pno_ind,:,ii),3); S_80px=sum(S_g(trial_80Pnx_ind,:,ii),3);
        S_00xx = sum(S_g(trial_20R_ind,:,ii),3);        
        S_xg = [ones(1,length(S_80ro)),2*ones(1,length(S_80po)),3*ones(1,length(S_00xx)),4*ones(1,length(S_80rx)),5*ones(1,length(S_80px))];
        [p,~,stats]=anovan([S_80ro;S_80po;S_00xx;S_80rx;S_80px],{S_xg},'display','off','alpha', 0.01);
        [c,~,~,~] = multcompare(stats,'display','off','alpha', 0.01);
        else
        S_80r = sum(S_g(trial_80R_ind,:,ii),3);
        S_80p = sum(S_g(trial_80P_ind,:,ii),3); 
        S_00x = sum(S_g(trial_20R_ind,:,ii),3);        
        S_xg = [ones(1,length(S_80r)),2*ones(1,length(S_80p)),3*ones(1,length(S_00x))];
        [p,~,stats]=anovan([S_80r;S_80p;S_00x;],{S_xg},'display','off','alpha', 0.01);
        [c,~,~,~] = multcompare(stats,'display','off','alpha', 0.01);
        end
        p_set(ii,jj) = p;
        c_set{ii,jj} = c;
    end
            disp(jj);
end
%% temporal shuffle
clear S_shuffle
for oo = 1:size(S_event,3)
for ii = 1:100
S_shuffle_temp = S_event(:,randperm(360),oo);
S_shuffle(:,:,oo,ii) = [sum(S_shuffle_temp(:,0*fr+1:0.5*fr),2),sum(S_shuffle_temp(:,0.5*fr+1:2*fr),2),...
                       sum(S_shuffle_temp(:,2*fr+1:3.5*fr),2),sum(S_shuffle_temp(:,4*fr+1:5.5*fr),2),...
                       sum(S_shuffle_temp(:,5.5*fr+1:8*fr),2),sum(S_shuffle_temp(:,8*fr+1:11*fr),2)];
end
disp(oo)
end

%% temporal shuffle test
S_test =  [sum(mean(S_C)>mean(S_shuffle(:,2,:,:)),4),sum(mean(S_D1)>mean(S_shuffle(:,3,:,:)),4),...
    sum(mean(S_D2)>mean(S_shuffle(:,4,:,:)),4),sum(mean(S_O)>mean(S_shuffle(:,5,:,:)),4),...
    sum(mean(S_ITI_a)>mean(S_shuffle(:,6,:,:)),4)];

S_test = permute(S_test,[3 2 1]);
save('S_shuffle.mat','S_shuffle','p_set','c_set','S_test');
%%
a = S_test>95;
actN1 = sum(a(sum(a,2)==1,:));
actN2 = sum(a(sum(a,2)==2,:));sum(a(sum(a,2)==3,:));

%% multiple regression part 1- Craw
[X4,Y4] = meshgrid(0*fr:12*fr-1,GPIO1onframe); %-0.25 ~ 5.25 ; 6~105 frame
set_regressionframe = X4+Y4;
%% lick for one
period = [0.5 2 3.5 4 5.5 8 11];
clear licktemp lick_reg_temp
lickinfo = timeline{1,1};
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
for tt = 1:length(period)-1
    lick_reg(:,tt) = sum(lick_reg_temp(:,period(tt)*10+1:period(tt+1)*10),2);
end

%% lick for rev
period = [0.5 2 3.5 4 5.5 8 11];
clear licktemp lick_reg_temp
lickinfo = [timeline{1,1};timeline{1,2}];
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
for tt = 1:length(period)-1
    lick_reg(:,tt) = sum(lick_reg_temp(:,period(tt)*10+1:period(tt+1)*10),2);
end
%% multiple regression part2- Craw
warning off;
period = [0.5 2 3.5 4 5.5 8 11];

%cue, reward, punish
clear rew_regTmp pun_regTmp
outcome = waterReward;
outcome(waterReward==0 & ismember(odorCue, [0 1]))=-1;
for ii = 1:size(odorCue,1);
    rew_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==1)& (rwProb(ii,odorCue(ii)+1)==75);
    pun_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==2)& (rwProb(ii,odorCue(ii)+1)==75);
end
set_crp = [odorCue, outcome.*rew_regTmp', outcome.*pun_regTmp'];
set_crp = set_crp(set_rectrial,:);

cue1_reg = repmat(set_crp(:,1)==(find(trial_kind==1)-1),1,120);% Rw
cue2_reg = repmat(set_crp(:,1)==(find(trial_kind==4)-1),1,120);% Pn
cue1_t1_reg = repmat([0; cue1_reg(1:end-1,1)],1,120);% Rw
cue2_t1_reg = repmat([0; cue2_reg(1:end-1,1)],1,120);% Pn
cue1_t2_reg = repmat([0;0; cue1_reg(1:end-2,1)],1,120);% Rw
cue2_t2_reg = repmat([0;0; cue2_reg(1:end-2,1)],1,120);% Pn

rew_reg = repmat(set_crp(:,2),1,120);
pun_reg = repmat(set_crp(:,3),1,120);
rew_1_reg =repmat([0; set_crp(1:end-1,2)],1,120);
pun_1_reg =repmat([0; set_crp(1:end-1,3)],1,120);
rew_2_reg =repmat([0;0; set_crp(1:end-2,2)],1,120);
pun_2_reg =repmat([0;0; set_crp(1:end-2,3)],1,120);
%% Period src beta - rwpn
%cut with neuron event
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_c = mean(C_event(:,0.5*fr+1:2*fr,:),2); %cue
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %delay1
C_O = mean(C_event(:,5.5*fr+1:8*fr,:),2); %Outcome
C_ITI = mean(C_event(:,5.5*fr+1:6.5*fr,:),2); %outcome 1s

if exist('cylinderspeed','var');
    speedset = [mean(cylinderspeed(set_rectrial,21:35),2),mean(cylinderspeed(set_rectrial,56:80),2),mean(cylinderspeed(set_rectrial,56:65),2)];
else
    speedset = zeros(size(set_rectrial,2),3); end

lickset = [mean(lick_reg(set_rectrial,21:35),2),mean(lick_reg(set_rectrial,56:80),2),mean(lick_reg(set_rectrial,56:65),2)];

clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset
for iper=1:3
    if iper==1; craw_reg = C_d(set_rectrial,:,:);
    elseif iper==2; craw_reg = C_O(set_rectrial,:,:);
    elseif iper==3; craw_reg = C_ITI(set_rectrial,:,:);end;
    for oo = 1:size(craw_reg,3)
        tt = 1;
        tbl = table(cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lickset(:,iper),speedset(:,iper),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),...            
            craw_reg(:,:,oo),...
            'VariableNames',{'Cue1','Cue2','RW','PN','Lick','Speed','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
        mdl = fitlm(tbl,...
            'C_raw~Cue1+RW+Cue2+PN+Lick+Speed+Cue1t_1+Rwt_1+Cue2t_1+Pnt_1+Cue1t_2+Rwt_2+Cue2t_2+Pnt_2'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lickset(:,iper),speedset(:,iper),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)]);
        ystd(1,tt) = std(craw_reg(:,:,oo));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        
        srcset(:,oo,iper) = src;
        betaset(:,oo,iper) = beta;
        pvalueset(:,oo,iper) = pvalue;
        
    end
    disp(iper)
end
save ('period_reg_20211118.mat','srcset','betaset', 'pvalueset')%,'xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')

%% Period src beta - rwprob
%cut with neuron event
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_c = mean(C_event(:,0.5*fr+1:2*fr,:),2); %cue
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %delay1
C_O = mean(C_event(:,5.5*fr+1:8*fr,:),2); %Outcome
C_ITI = mean(C_event(:,8*fr+1:11*fr,:),2); %ITI

if exist('cylinderspeed','var');
    speedset = [mean(cylinderspeed(set_rectrial,21:35),2),mean(cylinderspeed(set_rectrial,56:80),2),mean(cylinderspeed(set_rectrial,81:110),2)];
else
    speedset = zeros(size(set_rectrial,2),2); end

lickset = [mean(lick_reg(set_rectrial,21:35),2),mean(lick_reg(set_rectrial,56:80),2),mean(lick_reg(set_rectrial,81:110),2)];

clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset
for iper=1:3
    if iper==1; craw_reg = C_d(set_rectrial,:,:);
    elseif iper==2; craw_reg = C_O(set_rectrial,:,:);
    elseif iper==3; craw_reg = C_ITI(set_rectrial,:,:);end;
    for oo = 1:size(craw_reg,3)
        tt = 1;
        tbl = table(outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt), ...
            lickset(:,iper),speedset(:,iper),craw_reg(:,:,oo),...
            'VariableNames',{'o','ot_1','ot_2','v','vt_1','vt_2','Lick','speed','C_raw'});
        mdl = fitlm(tbl,'C_raw~o+ot_1+ot_2+v+vt_1+vt_2+Lick+speed'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt),lickset(:,iper), speedset(:,iper)]);
        ystd(1,tt) = std(craw_reg(:,:,oo));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        
        srcset(:,oo,iper) = src;
        betaset(:,oo,iper) = beta;
        pvalueset(:,oo,iper) = pvalue;
        
    end
    disp(iper)
end
save ('period_reg_20211118.mat','srcset','betaset', 'pvalueset')%,'xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')

%% RPE - rwprob
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_ITI = mean(C_event(:,8*fr+1:11*fr,:),2); %ITI
C_O1 = mean(C_event(:,5.5*fr+1:7*fr,:),2); %Outcome part
C_Oall = mean(C_event(:,5.5*fr+1:8*fr,:),2); %Outcome all

clear beta pvalue xstd ystd src srcset_RPE betaset_RPE pvalueset_RPE

if exist('cylinderspeed','var');
    speedset = [mean(cylinderspeed(set_rectrial,81:110),2),mean(cylinderspeed(set_rectrial,56:70),2),mean(cylinderspeed(set_rectrial,56:80),2) ];
else
    speedset = zeros(size(set_rectrial,2),2); end
lickset = [zeros(length(set_rectrial),1),mean(lick_reg(set_rectrial,56:70),2),mean(lick_reg(set_rectrial,56:80),2),];

for iper=1:2
    if iper==1; craw_reg = C_ITI(set_rectrial,:,:);
    elseif iper==2; craw_reg = C_O1(set_rectrial,:,:);end;
    for oo = 1:size(craw_reg,3)
        tt = 1;
        tbl = table(outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt), ...
            lickset(:,iper),speedset(:,iper),craw_reg(:,:,oo),...
            'VariableNames',{'o','ot_1','ot_2','v','vt_1','vt_2','Lick','speed','C_raw'});
        mdl = fitlm(tbl,'C_raw~o+ot_1+ot_2+v+vt_1+vt_2+Lick+speed'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt), ...
            lickset(:,iper),speedset(:,iper)]);
        ystd(1,tt) = std(craw_reg(:,:,oo));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        
        srcset_RPE(:,oo,iper,1) = src;
        betaset_RPE(:,oo,iper,1) = beta;
        pvalueset_RPE(:,oo,iper,1) = pvalue;
        
    end
    disp(iper)
end
clear beta pvalue xstd ystd src
sel_trial = [waterReward(:,:)==1,waterReward(:,:)==0];
% sel_trial = [trial_80Rwo_ind | trial_50Rwo_ind,trial_80Rwx_ind | trial_50Rwx_ind];
for iper=1:2
    if iper==1; craw_reg = C_ITI(set_rectrial,:,:);
    elseif iper==2; craw_reg = C_O1(set_rectrial,:,:);end;
    for irw = 1:2        
        for oo = 1:size(craw_reg,3)
            tt = 1;
            tbl = table(outcome_reg(sel_trial(:,irw),tt),outcome_1_reg(sel_trial(:,irw),tt),outcome_2_reg(sel_trial(:,irw),tt),value_reg(sel_trial(:,irw),tt),value_1_reg(sel_trial(:,irw),tt),value_2_reg(sel_trial(:,irw),tt), ...
                lickset(sel_trial(:,irw),iper),speedset(sel_trial(:,irw),iper),craw_reg(sel_trial(:,irw),:,oo),...
                'VariableNames',{'o','ot_1','ot_2','v','vt_1','vt_2','Lick','speed','C_raw'});
            mdl = fitlm(tbl,'C_raw~ot_1+ot_2+v+vt_1+vt_2+Lick+speed'); %
            beta(:,tt) = mdl.Coefficients.Estimate;
            pvalue(:,tt) = mdl.Coefficients.pValue;
            xstd(:,tt) = std([outcome_1_reg(sel_trial(:,irw),tt),outcome_2_reg(sel_trial(:,irw),tt),value_reg(sel_trial(:,irw),tt),value_1_reg(sel_trial(:,irw),tt),value_2_reg(sel_trial(:,irw),tt), ...
                lickset(sel_trial(:,irw),iper),speedset(sel_trial(:,irw),iper)]);
            ystd(1,tt) = std(craw_reg(sel_trial(:,irw),:,oo));
            src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
            
            srcset_RPE(:,oo,iper,1+irw) = [0; src];
            betaset_RPE(:,oo,iper,1+irw) = [0; beta];
            pvalueset_RPE(:,oo,iper,1+irw) = [1; pvalue];
            
        end
    end
    disp(iper)
end
save('RPEset_20220309.mat','srcset_RPE','betaset_RPE','pvalueset_RPE')
%% RPE - rwpn
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_ITI = mean(C_event(:,8*fr+1:11*fr,:),2); %ITI
C_O1 = mean(C_event(:,5.5*fr+1:7*fr,:),2); %Outcome 1.5s
C_O2 = mean(C_event(:,5.5*fr+1:6.5*fr,:),2); %Outcome 1s

clear beta pvalue xstd ystd src srcset_RPE betaset_RPE pvalueset_RPE

if exist('cylinderspeed','var');
    speedset = [mean(cylinderspeed(set_rectrial,81:110),2),mean(cylinderspeed(set_rectrial,56:70),2),mean(cylinderspeed(set_rectrial,56:80),2) ];
else
    speedset = zeros(size(set_rectrial,2),2); end
lickset = [mean(lick_reg(set_rectrial,56:70),2),mean(lick_reg(set_rectrial,56:65),2)];

for iper=1:2
    if iper==1; craw_reg = C_O1(set_rectrial,:,:);
    elseif iper==2; craw_reg = C_O2(set_rectrial,:,:); end;
    for oo = 1:size(craw_reg,3)
        tt = 1;
        tbl = table(cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lickset(:,iper),speedset(:,iper),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt),...
            craw_reg(:,:,oo),...
            'VariableNames',{'Cue1','Cue2','RW','PN','Lick','Speed','Cue1t_1','Cue2t_1','Rwt_1','Pnt_1','Cue1t_2','Cue2t_2','Rwt_2','Pnt_2','C_raw'});
        mdl = fitlm(tbl,...
            'C_raw~Cue1+RW+Cue2+PN+Lick+Speed+Cue1t_1+Rwt_1+Cue2t_1+Pnt_1+Cue1t_2+Rwt_2+Cue2t_2+Pnt_2'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([cue1_reg(:,tt), cue2_reg(:,tt),rew_reg(:,tt), pun_reg(:,tt), lickset(:,iper),speedset(:,iper),...
            cue1_t1_reg(:,tt),cue2_t1_reg(:,tt),rew_1_reg(:,tt),pun_1_reg(:,tt),...
            cue1_t2_reg(:,tt),cue2_t2_reg(:,tt),rew_2_reg(:,tt),pun_2_reg(:,tt)]);
        ystd(1,tt) = std(craw_reg(:,:,oo));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        
        srcset_RPE(:,oo,iper,1) = src;
        betaset_RPE(:,oo,iper,1) = beta;
        pvalueset_RPE(:,oo,iper,1) = pvalue;
        
    end
    disp(iper)
end

save('RPEset.mat','srcset_RPE','betaset_RPE','pvalueset_RPE')
%% Reversal - cue period
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,0.5*fr+1:3.5*fr,:),2); %cue +delay1
C_c = mean(C_event(:,0.5*fr+1:3.5*fr,:),2); %cue +delay1
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %cue +delay1
beforetrial = rev_index(1,5)-100;
if ismember(trial_kind,[0 1 3 4]); %rwpn
    if ~exist('aftertrial','var')
        for itrial = rev_index(1,5):nTrial-100 % find after trial
            [p,~,stats] = anovan(lickNum(itrial:itrial+100),odorCue(itrial:itrial+100),'display','off');
            if p<0.05;
                c = multcompare(stats);
                rwcue = find(trial_kind==4)-1;
                rwcueindex = find(strcmp(stats.grpnames{1,1},num2str(rwcue)));
                highindex = ismember(c(:,1),rwcueindex) | ismember(c(:,2),rwcueindex);
                if c(highindex,6)<0.05 & [c(ismember(c(:,1),rwcueindex),3)>0, c(ismember(c(:,2),rwcueindex),3)<0];
                    aftertrial = itrial;
                    break;
                end
            end
        end
        save('behavior_align','aftertrial','beforetrial','-append')
    end
    cueanova = trial_80R_ind+-1*trial_80P_ind;
    outanova = (odorCue==find(trial_kind==1)-1).*(outcomeContingency(:,find(trial_kind==1)))+...
        (odorCue==find(trial_kind==4)-1).*(outcomeContingency(:,find(trial_kind==4)));
    %cueanova(cueanova==0)=[];outanova(outanova==0)=[];
elseif ismember(trial_kind, [0 1 2 3]); %rwprob
    if ~exist('aftertrial','var')
        for itrial = rev_index(1,5):nTrial-100 % find after trial
            [p,~,stats] = anovan(lickNum(itrial:itrial+100),odorCue(itrial:itrial+100),'display','off');
            if p<0.05;
                c = multcompare(stats);
                rw80cue = find(trial_kind==2)-1;
                rw50cue = find(trial_kind==1)-1;
                rw20cue = find(trial_kind==3)-1;
                rw80cueindex = find(strcmp(stats.grpnames{1,1},num2str(rw80cue)));
                rw50cueindex = find(strcmp(stats.grpnames{1,1},num2str(rw50cue)));
                rw20cueindex = find(strcmp(stats.grpnames{1,1},num2str(rw20cue)));
                if c(:,6)<0.05 ...
                        & [c(ismember(c(:,1),rw80cueindex),3)>0; c(ismember(c(:,2),rw80cueindex),3)<0;...
                        any([c(ismember(c(:,1),rw50cueindex)&ismember(c(:,2),rw20cueindex),3)>0, c(ismember(c(:,1),rw20cueindex)&ismember(c(:,2),rw50cueindex),3)<0])];
                    aftertrial = itrial;
                    break;
                end
            end
        end
        save('behavior_align','aftertrial','beforetrial','-append')
    end
    
    cueanova = 75*trial_80R_ind+25*trial_50R_ind;
    outanova = (odorCue==find(trial_kind==1)-1).*(rwProb(:,find(trial_kind==1)))+...
        (odorCue==find(trial_kind==2)-1).*(rwProb(:,find(trial_kind==2)));
end


analysistrial= [beforetrial:beforetrial+99,aftertrial:aftertrial+99];

clear panovaset
for ii = 1:size(C_raw,1)
    [p,tbl,stats] = anovan(C_d(analysistrial,:,ii),{cueanova(analysistrial),outanova(analysistrial)},'model','full','varnames',{'Identity','Contingency'},...
        'display','off');
    panovaset(ii,:) = p;
%     [p,tbl,stats] = anovan(C_c(analysistrial,:,ii),{cueanova(analysistrial),outanova(analysistrial)},'varnames',{'Identity','Contingency'},...
%         'display','off');
%     panovaset(ii,:) = p;
%     [p,tbl,stats] = anovan(C_d(analysistrial,:,ii),{cueanova(analysistrial),outanova(analysistrial)},'varnames',{'Identity','Contingency'},...
%         'display','off');
%     panovaset(ii,:) = p;
end

endpien = [sum(panovaset(:,2)>=0.05&panovaset(:,1)<0.05), sum(panovaset(:,2)<0.05&panovaset(:,1)>=0.05),...
    sum(panovaset(:,2)<0.05&panovaset(:,1)<0.05), sum(panovaset(:,2)>=0.05&panovaset(:,1)>=0.05)];
save('endpien.mat','endpien','panovaset');
pie(endpien)
colormap([65 105 225; 255 99 71;138 43 226;  130 130 130;]./255)
legend('Identity','Contingency','Both','None')
%% draw plot
for oo = 1:size(C_raw,1)
    C_tmp = reshape(C_raw(oo,set_eventframe(:)),[],12*fr);
    C_eventset(:,:,oo) = C_tmp;
end
%% cue dependency
load('behavior_align.mat','C_z','odorCue','trial_kind','trial_80R_ind','trial_80P_ind','trial_50R_ind','trial_20R_ind','rev_b','fr')
C_cue = mean(C_z(:,0.5*fr+1:3.5*fr,:),2);
C_cuerw = C_cue.*(odorCue==find(trial_kind==1)-1); C_cuerw(odorCue==find(trial_kind==3)-1,:,:)=[];
C_cuepn = C_cue.*(odorCue==find(trial_kind==4|trial_kind==2)-1); C_cuepn(odorCue==find(trial_kind==3)-1,:,:)=[];
trial_rw = trial_80R_ind; trial_rw(trial_20R_ind) = [];
trial_pn = trial_50R_ind; trial_pn(trial_20R_ind) = [];
clear C_cuedep
for ii= 1: length(trial_rw)-25
    C_cuedep(ii,:,:) = sum(C_cuerw(ii:ii+24,:,:))./sum(trial_rw(ii:ii+24))-sum(C_cuepn(ii:ii+24,:,:))./sum(trial_pn(ii:ii+24));
end
ctmp = C_cuedep';
% ctmp = ctmp(panovaset(:,1)<0.05,:);
ctmp(:,length(trial_rw)-24) = mean(ctmp(:,1:sum(rev_b)-sum(trial_20R_ind(1:sum(rev_b)))),2);
figure
imagesc(sortrows(ctmp,length(trial_rw)-24))
line([sum(rev_b)-sum(trial_20R_ind(1:sum(rev_b))) sum(rev_b)-sum(trial_20R_ind(1:sum(rev_b)))],[0 800],'color','k','linewidth',5)
 set(gca,'clim', [-1 1]); 
%% correlation ITI & whole
[X,Y] = meshgrid(0:12*fr-1,GPIO1onframe); % 0~12
set_eventframe = X+Y; %meshgrid for event frame
S_event = S_raw(:,set_eventframe(:));
S_event = reshape(S_event', [],12*fr,size(S_event,1));
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
S_event01 = S_event>0;
% ITI coactivated cell
coaThr = 0.05*size(S_raw,1);
S_ITI=S_event01(:,8*fr+1:11*fr,:);
clear S_co ITI_cell
for ii =1:size(S_ITI,2)-15;
    S_co(:,ii,:)=sum(S_ITI(:,ii:ii+14,:),2)>0;
end
[S_coind,S_cotime] = max(sum(S_co,3),[],2);
S_coa = S_event01(S_coind>=coaThr,:,:);
for ii = 1:size(S_event01,1)
    ITI_cell(ii,:,:) = sum(S_ITI(ii,S_cotime(ii):S_cotime(ii)+14,:),2)>0; % max time fire cell numbering
end

ITI_trial.rwo = find((S_coind>=coaThr) & trial_80Rwo_ind); [rro, cro] = find(ITI_cell(ITI_trial.rwo,:)>0); %colume = cell number
ITI_trial.rwx = find((S_coind>=coaThr) & trial_80Rwx_ind); [rrx, crx] = find(ITI_cell(ITI_trial.rwx,:)>0);
ITI_trial.pno = find((S_coind>=coaThr) & trial_80Pno_ind); [rpo, cpo] = find(ITI_cell(ITI_trial.pno,:)>0);
ITI_trial.pnx = find((S_coind>=coaThr) & trial_80Pnx_ind); [rpx, cpx] = find(ITI_cell(ITI_trial.pnx,:)>0);
ITI_trial.nto = find((S_coind>=coaThr) & trial_20R_ind);   [roo, coo] = find(ITI_cell(ITI_trial.nto,:)>0);
%%
S_ITI_cor = movsum(S_ITI,[0 14],2,'Endpoints','discard');
S_ITI_cor = S_ITI_cor(1:sum(rev_b),1:15:size(S_ITI_cor,2),:);
S_cor = movsum(S_raw,[0 14],2,'Endpoints','discard');
S_cor= S_cor(:,1:15:size(S_cor,2));
clear Rset_all Pset_all
for icell = 1:size(S_cor,1);
    for jcell = 1:size(S_cor,1);
        [R,P]  = corrcoef(S_cor(icell,:),S_cor(jcell,:));
        Rset_all(icell,jcell) = R(1,2); Pset_all(icell,jcell)=P(1,2); 
        [R,P]  = corrcoef(reshape(S_ITI_cor(:,:,icell)',[],1),reshape(S_ITI_cor(:,:,jcell)',[],1));
        Rset_ITI(icell,jcell) =R(1,2); Pset_ITI(icell,jcell)=P(1,2);
    end
end
%%
triallist = {'ro','rx','po','px','oo'};
for itrial = 1:5
    eval(['celllist=unique(c',triallist{itrial},');']);
    cellcomb = nchoosek(unique(celllist),2);
    clear cortemp
    for icell = 1:size(cellcomb,1)
        cortemp(icell,1) = Rset_all(cellcomb(icell,1),cellcomb(icell,2));
    end
    disp(length(unique(celllist)))
    eval(['cor_',triallist{itrial},'= cortemp;'])
end
[mean(cor_ro),mean(cor_rx),mean(cor_po),mean(cor_px),mean(cor_oo);]
%% cue activity compare (r,p,0)
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,2*fr+1:3.5*fr,:),2); %cue +delay1
% C_o = mean(C_event(:,5.5*fr+1:11*fr,:),2); %cue +delay1
cueorder = zeros(size(C_cd,3),3);
for oo = 1:size(C_cd,3)
    if ismember(trial_kind,[0 1 3 4]); %rwpn
        [p,~,stats] = anovan(C_cd(:,:,oo), trial_80R_ind+2*trial_80P_ind+3*trial_20R_ind,'display','off');
        cuesetname = {'R','P','N','R>P>0','R>0>P','P>R>0','P>0>R','0>P>R','0>R>P'};
    elseif ismember(trial_kind, [0 1 2 3]); %rwprob
        [p,~,stats] = anovan(C_cd(:,:,oo), trial_80R_ind+2*trial_50R_ind+3*trial_20R_ind,'display','off');
        cuesetname = {'75','25','0','75>25>0','75>0>25','25>75>0','25>0>75','0>25>75','0>75>25'};
    end
    cueorder(oo,1) = p;
    if p<0.05
        c = multcompare(stats,'display','off');
        if sum(c(:,6)<0.05)==3
            cueorder(oo,2) = 1;
        elseif c([1 2],6)<0.05
            cueorder(oo,3) = 1;
        end
        %         if sum(c(:,6)<0.05)==2 | sum(c(:,6)<0.05)==1
        %             if c([1,2],4)>0
        %                 cueorder(1,1) = cueorder(1,1)+1;
        %             elseif c(1,4)<0 & c(3,4)>0
        %                 cueorder(1,2) = cueorder(1,2)+1;
        %             elseif c([2,3],4)<0
        %                 cueorder(1,3) = cueorder(1,3)+1;
        %             end
        %         elseif sum(c(:,6)<0.05)==3
        %             if c(:,4)>0
        %                 cueorder(1,4) = cueorder(1,4)+1;
        %             elseif c([1,2],4)>0 & c([3],4)<0
        %                 cueorder(1,5) = cueorder(1,5)+1;
        %             elseif c([2,3],4)>0 & c([1],4)<0
        %                 cueorder(1,6) = cueorder(1,6)+1;
        %             elseif c([3],4)>0 & c([1,2],4)<0
        %                 cueorder(1,7) = cueorder(1,7)+1;
        %             elseif c(:,4)<0
        %                 cueorder(1,8) = cueorder(1,8)+1;
        %             elseif c([1],:)>0 & c([2,3],4)<0
        %                 cueorder(1,9) = cueorder(1,9)+1;
        %             end
        %         end
        
    end
end

 save('cue_order.mat','cueorder','cuesetname')
% save('cue_order.mat','cueorder','cuesetname')
% dc = cat(1,cueorder_dHP03,cueorder_dHP04,cueorder_dHP06,cueorder_dHP08);
% vc = cat(1,cueorder_vHP06,cueorder_vHP07,cueorder_vHP08,cueorder_vHP10,cueorder_vHP11);


%% z-score
C_mean = mean(C_raw(:,GPIO1onframe(1):GPIO1onframe(end)),2);
C_std = std(C_raw(:,GPIO1onframe(1):GPIO1onframe(end)),0,2);
C_n = (C_raw-C_mean)./C_std;
%eventframe
C_event = C_n(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_O = mean(C_event(:,5.5*fr+1:8*fr,:),2);

%%
rwDepAct = mean(C_O(trial_80Rwo_ind|trial_50Rwo_ind,:,:))-mean(C_O(trial_80Rwx_ind|trial_50Rwx_ind,:,:));
pnDepAct = mean(C_O(trial_80Pno_ind,:,:))-mean(C_O(trial_80Pnx_ind,:,:));
for ii = 1:100
    ranii = randperm(size(C_O,1));
    C_rand = C_O(ranii,:,:);
    rwDepAct_rand = mean(C_rand(trial_80Rwo_ind|trial_50Rwo_ind,:,:))-mean(C_rand(trial_80Rwx_ind|trial_50Rwx_ind,:,:));
    pertest(ii,:)= (rwDepAct>rwDepAct_rand)-(rwDepAct<rwDepAct_rand);        
end
sig = abs(sum(pertest))>97.5;
figure
hold on
xrange = 1:length(rwDepAct);
scatter(xrange,rwDepAct)
scatter(xrange(sig),rwDepAct(sig),'r')
%  clearvars -except rwDepAct
%%
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_cd = mean(C_event(:,0.5*fr+1:3.5*fr,:),2); %cue +delay1
C_c = mean(C_event(:,0.5*fr+1:2*fr,:),2); %cue +delay1
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %cue +delay1
[p,~,stats] = anovan(C_d([beforetrial:beforetrial+99, beforetrial+100:beforetrial+199,aftertrial:aftertrial+99]),...
    {[ones(1,100),2*ones(1,100),3*ones(1,100)]});
[p,stats] = anova_rm(C_d([beforetrial:beforetrial+99; beforetrial+100:beforetrial+199;aftertrial:aftertrial+99])');
figure
hold on
bar(1,mean(C_d(beforetrial:beforetrial+99)));
bar(2,mean(C_d(beforetrial+100:beforetrial+199)));
bar(3,mean(C_d(aftertrial:aftertrial+99)));
errorbar(1,mean(C_d(beforetrial:beforetrial+99)),sem(C_d(beforetrial:beforetrial+99)));
errorbar(2,mean(C_d(beforetrial+100:beforetrial+199)),sem(C_d(beforetrial+100:beforetrial+199)));
errorbar(3,mean(C_d(aftertrial:aftertrial+99)),sem(C_d(aftertrial:aftertrial+99)));

%% Reversal delay outcome ITI change
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %delay1
C_o = mean(C_event(:,5.5*fr+1:8*fr,:),2); %outcome
C_ITI = mean(C_event(:,8*fr+1:11*fr,:),2); %ITI
trialL= [1,2,3,4];
cmap = [0 0 205; 255 30 70; 80 80 80]./255;
C_rev = NaN(420,4,size(C_d,3));

for ii= 1:4
xtmp = find(odorCue==find(trial_kind==ii)-1);
C_rev(xtmp,ii,:,1) = C_d(odorCue==find(trial_kind==ii)-1,:,:);
C_rev(xtmp,ii,:,2) = C_o(odorCue==find(trial_kind==ii)-1,:,:);
C_rev(xtmp,ii,:,3) = C_ITI(odorCue==find(trial_kind==ii)-1,:,:);
end
C_rev_onset = C_rev(rev_index(1,5)-100:rev_index(1,5)+150);
C_rev_learned = C_rev(aftertrial-100:aftertrial+150);

save('reversal_z.mat','C_rev','C_rev_onset','C_rev_learned');

%%
% trialL= [1,2,3];
% cmap = [0 0 205; 129 212 250; 80 80 80]./255;
mkdir([cd,'\figure\C_trace'])
f1 = figure;
for oo= 1:size(C_d,3)
    for ii = 1:3
        xtmp = find(odorCue==find(trial_kind==trialL(1, ii))-1);
        ytmp_d = NaN(1,420);
        ytmp_d(xtmp) = C_d(odorCue==find(trial_kind==trialL(1, ii))-1,:,oo);
        ytmp_d = movmean(ytmp_d, 50, 'omitnan','Endpoints','discard');
        ytmp_ITI = NaN(1,420);
        ytmp_ITI(xtmp) = C_ITI(odorCue==find(trial_kind==trialL(1, ii))-1,:,oo);
        ytmp_ITI = movmean(ytmp_ITI, 50, 'omitnan','Endpoints','discard');
        ax1 = subplot(1,2,1);
        hold on
        plot(ytmp_d, 'Color',cmap(ii,:))
        
        ax2 = subplot(1,2,2);
        hold on
        plot(ytmp_ITI, 'Color',cmap(ii,:))
    end
    
    linkaxes([ax1 ax2],'y')
    if sum(rev_b)<420; 
        line(ax1,[beforetrial+100 beforetrial+100],[0 10], 'Color','k','lineStyle',':');
        line(ax2,[beforetrial+100 beforetrial+100],[0 10], 'Color','k','lineStyle',':');
        if exist('aftertrial')
        line(ax1,[aftertrial aftertrial],[0 10], 'Color','k','lineStyle','-.'); 
        line(ax2,[aftertrial aftertrial],[0 10], 'Color','k','lineStyle','-.');
        end
    end;
%     saveas(f1,[cd,'\figure\C_trace\',num2str(oo),'.tif'])
    clf
end

%%
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_d = mean(C_event(:,2*fr+1:3.5*fr,:),2); %cue +delay1
C_ITI = mean(C_event(:,8*fr+1:11*fr,:),2); %cue +delay1
figure; plot(mean(permute(C_d,[1 3 2]),2))
figure; plot(mean(permute(C_ITI,[1 3 2]),2))
figure; plot(mean(mean(C_event,2),3))
%%
trialL= [1,4,3];
cmap = [0 0 205; 255 30 70; 80 80 80]./255;
% trialL= [1,2,3];
% cmap = [0 0 205; 129 212 250; 80 80 80]./255;
mkdir([cd,'\figure\C_trace'])
f1 = figure;
for oo= 1:size(C_d,3)
    for ii = 1:3
        xtmp = find(odorCue==find(trial_kind==trialL(1, ii))-1);
        ytmp_all = NaN(1,420);
        ytmp_all(xtmp) = mean(C_event(odorCue==find(trial_kind==trialL(1, ii))-1,1:11*30,oo),2);
        ytmp_all = movmean(ytmp_all, 50, 'omitnan','Endpoints','discard');
        ytmp_d = NaN(1,420);
        ytmp_d(xtmp) = C_d(odorCue==find(trial_kind==trialL(1, ii))-1,:,oo);
        ytmp_d = movmean(ytmp_d, 50, 'omitnan','Endpoints','discard');
        ytmp_ITI = NaN(1,420);
        ytmp_ITI(xtmp) = C_ITI(odorCue==find(trial_kind==trialL(1, ii))-1,:,oo);
        ytmp_ITI = movmean(ytmp_ITI, 50, 'omitnan','Endpoints','discard');
        
        ax1 = subplot(1,3,1);
        hold on
        plot(ytmp_d, 'Color',cmap(ii,:))
        title('Delay')
        
        ax2 = subplot(1,3,2);
        hold on
        plot(ytmp_ITI, 'Color',cmap(ii,:))
        title('ITI')
        
        ax3 = subplot(1,3,3);
        hold on
        plot(ytmp_all, 'Color',cmap(ii,:))
        title('All')
        
    end
    
    linkaxes([ax1 ax2 ax3],'y')
    if sum(rev_b)<420; 
        line(ax1,[beforetrial+100 beforetrial+100],[0 10], 'Color','k','lineStyle',':');
        line(ax2,[beforetrial+100 beforetrial+100],[0 10], 'Color','k','lineStyle',':');
        line(ax3,[beforetrial+100 beforetrial+100],[0 10], 'Color','k','lineStyle',':');
        if exist('aftertrial')
        line(ax1,[aftertrial aftertrial],[0 10], 'Color','k','lineStyle','-.'); 
        line(ax2,[aftertrial aftertrial],[0 10], 'Color','k','lineStyle','-.');
        line(ax3,[aftertrial aftertrial],[0 10], 'Color','k','lineStyle','-.');
        end
    end;
    saveas(f1,[cd,'\figure\C_trace\',num2str(oo),'.tif'])
    clf
end


%% Delay v(t-1) v(t) sig
clear
load('behavior_align.mat','C_z_all','fr','odorCue','trial_kind','waterReward',...
    'trial_80R_ind','trial_50R_ind','trial_80P_ind','trial_20R_ind')
d_c = mean(C_z_all(:,0.5*fr+1:2*fr,:),2); d_c = permute(d_c,[1 3 2]);
d_m= mean(C_z_all(:,2*fr+1:3.5*fr,:),2); d_m = permute(d_m,[1 3 2]);

for oo = 1:size(d_m,2)
    [p,~,stats] = anovan(d_m(:,oo), odorCue,'display','off');
    siglist(oo,1)=p;
    s = multcompare(stats);
    if p<0.05 & s(:,6)<0.05; siglist(oo,2)=0; end
    [p,~,stats] = anovan(d_c(:,oo), odorCue,'display','off');
    siglist_cue(oo,1)=p;
    s = multcompare(stats);
    if p<0.05 & s(:,6)<0.05; siglist_cue(oo,2)=0; end
end
sigdelayneuron=siglist<0.05;
sigcueneuron= siglist_cue<0.05;
if ismember(trial_kind, [0 1 2 3])
    d_t(:,:,1) = [mean(d_m(trial_80R_ind,:),1); mean(d_m(trial_50R_ind,:),1); mean(d_m(trial_20R_ind,:),1)];
    d_t(:,:,2) = [mean(d_c(trial_80R_ind,:),1); mean(d_c(trial_50R_ind,:),1); mean(d_c(trial_20R_ind,:),1)];

elseif ismember(trial_kind,[0 1 3 4])
    d_t(:,:,1) = [mean(d_m(trial_80R_ind,:),1); mean(d_m(trial_80P_ind,:),1); mean(d_m(trial_20R_ind,:),1)];
    d_t(:,:,2) = [mean(d_c(trial_80R_ind,:),1); mean(d_c(trial_80P_ind,:),1); mean(d_c(trial_20R_ind,:),1)];
end;

save('CorRorV.mat','d_t','sigdelayneuron','sigcueneuron')

%% history effect
load('behavior_align.mat','waterReward','C_z_all','fr',...
    'trial_80R_ind','trial_80Rwo_ind','trial_80Rwx_ind',...
    'trial_50R_ind','trial_50Rwo_ind','trial_50Rwx_ind',...
    'trial_20R_ind','trial_20Rwo_ind','trial_20Rwx_ind',... 
    'trial_80P_ind','trial_80Pno_ind','trial_80Pnx_ind')
% reward
% all trial, cue+delay, cue, delay, outcome, ITI
region = [0 11; 0.5 3.5; 0.5 2; 2 3.5; 5.5 8; 8 11];
rwhistory{1,1} = C_z_all([false;waterReward(2:420)==1&waterReward(1:419)==1],:,:); %r-->r
rwhistory{1,2} = C_z_all([false;waterReward(2:420)==1&waterReward(1:419)==0],:,:); %u-->r
rwhistory{2,1} = C_z_all([false;waterReward(2:420)==0&waterReward(1:419)==1],:,:); %r-->u
rwhistory{2,2} = C_z_all([false;waterReward(2:420)==0&waterReward(1:419)==0],:,:); %u-->u

trial_all = [trial_80R_ind, trial_80Rwo_ind,trial_80Rwx_ind,...
    trial_50R_ind,trial_50Rwo_ind,trial_50Rwx_ind,...
    trial_20R_ind,trial_20Rwo_ind,trial_20Rwx_ind,...
    trial_80P_ind,trial_80Pno_ind,trial_80Pnx_ind];
save('historyeffect.mat','C_z_all','rwhistory','waterReward','trial_all')

%% Residual
load('behavior_align.mat','C_z_all','fr','cylinderspeed','timeline',...
    'trial_80R_ind','trial_50R_ind','trial_20R_ind','trial_80P_ind',...
    'trial_80Rwo_ind','trial_80Rwx_ind','trial_50Rwo_ind','trial_50Rwx_ind','trial_80Pno_ind','trial_80Pnx_ind');
Residual_speed_set = NaN(420,120,size(C_z_all,3),2);
lickinfo = timeline{1,1};
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms
if ~exist('cylinderspeed','var'); cylinderspeed = zeros(420,120); end;
for oo = 1:size(C_z_all,3);
    C_z_reg = movmean(C_z_all(:,:,oo),0.5*fr,2);
    for tt = 1:120
        tbl = table(cylinderspeed(:,tt), lick_reg(:,tt),mean(C_z_reg(:,3*tt-2:3*tt),2),...
            'VariableNames',{'Speed','Lick','C_raw'});

        mdl = fitlm(tbl,'C_raw~Speed+Lick'); %
        Residual_speed_set(:,tt,oo) = table2array(mdl.Residuals(:,1));
    end
    disp(oo)
end
for ii = 1:3
    c_tmp=Residual_speed_set(:,:,:,ii);
    C_residual_cue(:,:,:,ii)=cat(3,permute(mean(c_tmp(trial_80R_ind,0*10+1:5.5*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_50R_ind,0*10+1:5.5*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_20R_ind,0*10+1:5.5*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_80P_ind,0*10+1:5.5*10,:),1),[3 2 1]));
    C_residual_out(:,:,:,ii)=cat(3,permute(mean(c_tmp(trial_80Rwo_ind,5.5*10+1:11*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_80Rwx_ind,5.5*10+1:11*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_50Rwo_ind,5.5*10+1:11*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_50Rwx_ind,5.5*10+1:11*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_20R_ind,5.5*10+1:11*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_80Pno_ind,5.5*10+1:11*10,:),1),[3 2 1]),...
        permute(mean(c_tmp(trial_80Pnx_ind,5.5*10+1:11*10,:),1),[3 2 1]));
end
save('Residual.mat','Residual_speed_set','C_residual_cue','C_residual_out')

%% Separation lick cell
% clear 
% load('behavior_align.mat')
lick_bout_onset = [0; find(diff(lickTime(:,1))>1*10^6)];
lick_bout_offset = [find(diff(lickTime(:,1))>1*10^6);size(lickTime(:,1),1)];
lick_bout_offset([false;diff(lick_bout_onset)<2])=[];
lick_bout_onset([false;diff(lick_bout_onset)<2])=[];

lick_trial_index = sortrows(lick_trial,[2 1]);
lick_trial_index(lick_trial_index(:,1)<0,:)=[];
lick_onset_time = [lick_trial_index(lick_bout_onset+1,2),floor(lick_trial_index(lick_bout_onset+1,1)*30/10^6)];
lick_offset_time = [lick_trial_index(lick_bout_offset,2),floor(lick_trial_index(lick_bout_offset,1)*30/10^6)];
lick_mean_activity = NaN(size(lick_onset_time,1),2,size(C_z_all,3));
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
for oo = 1:size(C_z_all,3)
    for ii = 1:size(lick_onset_time,1)
        if lick_onset_time(ii,2)<166 % lick --> outcome
            if outcomeContingency(lick_onset_time(ii,1),odorCue(lick_onset_time(ii,1))+1) ==1 && ...
                    waterReward(lick_onset_time(ii,1))==1; %rewarded
                lick_mean_activity(ii,1,oo) = mean(C_event(lick_onset_time(ii,1),[-30:-1]+lick_onset_time(ii,2),oo),2);
                lick_mean_activity(ii,2,oo) = mean(C_event(lick_onset_time(ii,1),lick_onset_time(ii,2):min(165-3,lick_onset_time(ii,2)+29),oo),2);
            else %unrewarded, punishment, no outcome
                lick_mean_activity(ii,1,oo) = mean(C_event(lick_onset_time(ii,1),[-30:-1]+lick_onset_time(ii,2),oo),2);
                lick_mean_activity(ii,2,oo) = mean(C_event(lick_onset_time(ii,1),[0:29]+lick_onset_time(ii,2),oo),2);
            end
        elseif lick_onset_time(ii,2)>=166 % outcome --> lick
            if outcomeContingency(lick_onset_time(ii,1),odorCue(lick_onset_time(ii,1))+1) ==1 && ...
                    waterReward(lick_onset_time(ii,1))==1; %rewarded
            else %unrewarded, punishment, no outcome
                lick_mean_activity(ii,1,oo) = mean(C_event(lick_onset_time(ii,1),[-30:-1]+lick_onset_time(ii,2),oo),2);
                lick_mean_activity(ii,2,oo) = mean(C_event(lick_onset_time(ii,1),[0:29]+lick_onset_time(ii,2),oo),2);
            end
        end
        neuronlick(ii,:,oo) = C_z_all(lick_onset_time(ii,1),[-30:29]+lick_onset_time(ii,2),oo);
    end
    [~,p] = ttest(lick_mean_activity(:,1,oo),lick_mean_activity(:,2,oo));
    pset_lick(oo,1) = p;
%     pset_lick(oo,:)=[p,...
%         mean(lick_mean_activity(~isnan(lick_mean_activity(:,2,oo)),1,oo)<lick_mean_activity(~isnan(lick_mean_activity(:,2,oo)),2,oo))>0.7,...
%         mean(lick_mean_activity(~isnan(lick_mean_activity(:,2,oo)),1,oo)<lick_mean_activity(~isnan(lick_mean_activity(:,2,oo)),2,oo))>0.75];


end

save('lick_neuron.mat','pset_lick');

%% outcome activity with select trials
load('behavior_align.mat','trial_80Rwo_ind','trial_80Rwx_ind','trial_50Rwo_ind','trial_50Rwx_ind','trial_20R_ind','trial_80Pno_ind','trial_80Pnx_ind',...
    'C_z_nor_Outcome','trial_kind','fr')
load('lick_neuron.mat')
% outcome : 4~11 s

if ismember(trial_kind,[0 1 3 4]) %rwpn
    seltrialnum = 27;
    trialname = {'trial_80Rwo_ind','trial_80Rwx_ind','trial_80Pno_ind','trial_80Pnx_ind'};
elseif ismember(trial_kind,[0 1 2 3]) %rwprob
    seltrialnum = 24;
    trialname = {'trial_80Rwo_ind','trial_80Rwx_ind','trial_50Rwo_ind','trial_50Rwx_ind'};
end

for ii = 1:4
    eval(['indtmp = find(',trialname{1,ii},');']);
    indselect = indtmp(randperm(length(indtmp),seltrialnum));
    C_z_selected(:,:,:,ii) = C_z_nor_Outcome(indselect,4*fr+1:11*fr,pset_lick(:,1)>=0.05);
    indset{1,ii} = indselect;
end
save('C_z_selected_Nor.mat','C_z_selected','seltrialnum','indset')

%% Outcome period src beta - rwprob - lick control
%cut with neuron event
load('lick_neuron.mat'); load('C_z_selected.mat','indset');
C_event = C_raw(pset_lick(:,1)>=0.05,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_O1 = mean(C_event(:,5.5*fr+1:8*fr,:),2); %Outcome
C_O2 = mean(C_event(:,5.5*fr+1:6.5*fr,:),2); %Outcome
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset
seltrial = cell2mat(indset');

for oo = 1:size(C_event,3)
    craw_reg = movmean(C_event, 0.5*fr,2); %movmean for 500ms
    for tt = 1:120       
         tbl = table(outcome_reg(:,tt),outcome_1_reg(:,tt),outcome_2_reg(:,tt),value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt), ...
             lick_reg(:,tt),cylinderspeed(:,tt),mean(craw_reg(:,3*tt-2:3*tt),2),...
            'VariableNames',{'o','ot_1','ot_2','v','vt_1','vt_2','Lick','speed','C_raw'});
        mdl = fitlm(tbl,'C_raw~o+ot_1+ot_2+v+vt_1+vt_2+Lick+speed'); %
        beta(:,tt) = mdl.Coefficients.Estimate;
        pvalue(:,tt) = mdl.Coefficients.pValue;
        xstd(:,tt) = std([outcome_reg(:,tt), outcome_1_reg(:,tt), outcome_2_reg(:,tt), value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt),lick_reg(:,tt),cylinderspeed(:,tt)]);
        ystd(1,tt) = std(mean(craw_reg(:,3*tt-2:3*tt),2));
        src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
        X_svm(:,oo,tt) = mean(craw_reg(:,3*tt-2:3*tt),2);
        Y_svm(:,:,tt) = [outcome_reg(:,tt), outcome_1_reg(:,tt), outcome_2_reg(:,tt), value_reg(:,tt),value_1_reg(:,tt),value_2_reg(:,tt),lick_reg(:,tt),cylinderspeed(:,tt)];
    end
    srcset(:,:,oo) = src;
    betaset(:,:,oo) = beta;
    pvalueset(:,:,oo) = pvalue;
    disp(oo)

end

save ('Outcome_reg_20211118.mat','srcset','betaset', 'pvalueset')%,'xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')
%% Outcome period src beta - rwpn - lick control
%cut with neuron event
load('lick_neuron.mat'); load('C_z_selected.mat','indset');
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_O1 = mean(C_event(:,5.5*fr+1:8*fr,:),2); %Outcome
C_O2 = mean(C_event(:,5.5*fr+1:6.5*fr,:),2); %Outcome
clear beta pvalue xstd ystd src X_svm Y_svm srcset betaset pvalueset


for itrial = 1:2
    if itrial ==1; sel_trial= 1:420';
    elseif itrial==2; sel_trial = cell2mat(indset');sel_trial = sel_trial';end;
    if exist('cylinderspeed','var');
        speedset = [mean(cylinderspeed(sel_trial,56:80),2),mean(cylinderspeed(sel_trial,56:65),2)];
    else
        speedset = zeros(size(sel_trial,2),2); end
    lickset = [mean(lick_reg(sel_trial,56:80),2),mean(lick_reg(sel_trial,56:65),2)];
    for iper=1:2
        if iper==1; craw_reg = C_O1(sel_trial,:,pset_lick(:,1)>=0.05);
        elseif iper==2; craw_reg = C_O2(sel_trial,:,pset_lick(:,1)>=0.05);end;
        for oo = 1:size(craw_reg,3)
            tt = 1;
            tbl = table(rew_reg(sel_trial,tt), pun_reg(sel_trial,tt),lickset(:,iper),speedset(:,iper),craw_reg(:,:,oo),...
                'VariableNames',{'RW','PN','Lick','speed','C_raw'});
            mdl = fitlm(tbl,'C_raw~RW+PN+Lick+speed'); %
            beta(:,tt) = mdl.Coefficients.Estimate;
            pvalue(:,tt) = mdl.Coefficients.pValue;
            xstd(:,tt) = std([rew_reg(sel_trial,tt), pun_reg(sel_trial,tt),lickset(:,iper), speedset(:,iper)]);
            ystd(1,tt) = std(craw_reg(:,:,oo));
            src(:,tt) = beta(2:end,tt).*xstd(:,tt)./ystd(1,tt);
            
            srcset(:,oo,iper,itrial) = src;
            betaset(:,oo,iper,itrial) = beta;
            pvalueset(:,oo,iper,itrial) = pvalue;
            
        end
        disp(iper)
    end
end
save ('Outcome_reg_20211118.mat','srcset','betaset', 'pvalueset')%,'xstd','ystd','src','X_svm','Y_svm','srcset','betaset','pvalueset')

%% k-means clustering
fr=30;
cmap1= [30 45 232; 120 177 255; 125 125 125; 255 30 70]./255;
load('E:\data\vHPC\all\rwpn\C_z_selected.mat')
C_z_selected_dHP = cat(3,C_z_selected_dHP03,C_z_selected_dHP04,C_z_selected_dHP06,C_z_selected_dHP07,C_z_selected_dHP08,C_z_selected_dHP10);%
C_z_selected_vHP = cat(3,C_z_selected_vHP06,C_z_selected_vHP07,C_z_selected_vHP08,C_z_selected_vHP11,C_z_selected_vHP12,C_z_selected_vHP14);%
% pn onset = 1.5 fr
C_z_pn = mean(C_z_selected_vHP(:,0.5*fr+1:2.5*fr,:,3));- mean(mean(C_z_selected_vHP(:,0.5*fr+1:2.5*fr,:,4),1),2);
C_z_pnfigure = mean(C_z_selected_vHP(:,:,:,3));% - C_z_selected_dHP(:,:,:,4),1);
X = permute(C_z_pn,[3 2 1]);
% [pcacoeff,pcascore,~,~,explained,~] = pca(permute(C_z_pn,[3 2 1]));
% X = pcascore(:,1:3);
% [pcacoeff,pcascore,~,~,explained,~] = pca(permute(C_z_pn,[2 3 1]));
% X = pcacoeff(:,1:3);
rng(1);
%% find best group number
distset = NaN(10,2);
for ii = 1:10
[idx,~,~,D] = kmeans(X,ii,'Replicates',10);
dist = min(D,[],2);
[silh] = silhouette(X,idx);
distset(ii,:) = [sum(dist),nanmean(silh)];
end

%% flgure 
% load('E:\data\vHPC\all\rwpn\kclustering.mat','idx','C_z_pn');
for ii = 2:8
    [idx] = kmeans(X,ii,'Replicates',10);
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 15]);    
    imagesort = sortrows([permute(C_z_pnfigure,[3 2 1]),idx],size(C_z_pnfigure,2)+1);
    hold on
    imagesc(imagesort(:,1:size(C_z_pnfigure,2)))
    caxis([-1 1]);
    yticklist = [];
    for iline = 1:ii
        line([1 210],[sum(idx<=iline) sum(idx<=iline)],'linestyle',':','color','k','linewidth',2)
        yticklist(iline)=sum(idx<=iline);
    end
    ax = gca;
    ax.YDir = 'reverse';
    yticks(yticklist)
    line([1.5*fr 1.5*fr], [0 1000],'linestyle',':','color','k','linewidth',2)
    line([4*fr 4*fr], [0 1000],'linestyle',':','color','k','linewidth',2)
    xticks([1 1.5*fr 4*fr]); xticklabels([-1.5 0 2.5]);
    xlim([0 210]); ylim([1 size(C_z_pn,3)]);
%     saveas(f1,['E:\data\vHPC\all\figure\cluster\dtrace_clusterbypn',num2str(ii),'.tif']);
    
    f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 15]);
    C_z_pn_sort = permute(mean(C_z_selected_dHP(:,:,:,3),1),[3 2 1]);
    C_z_nopn_sort = permute(mean(C_z_selected_dHP(:,:,:,4),1),[3 2 1]);
    for iplot = 1:ii
        subplot(ii,1,iplot)
        hold on
        stdshade(C_z_pn_sort(idx==iplot,:),0.3,[255 30 70]./255)
        stdshade(C_z_nopn_sort(idx==iplot,:),0.3,[255 30 70]./255,[],':')
        line([1.5*fr 1.5*fr], [-1 1],'linestyle',':','color','k','linewidth',2)
        line([4*fr 4*fr], [-1 1],'linestyle',':','color','k','linewidth',2)
        xticks([1 1.5*fr 4*fr]); xticklabels([-1.5 0 2.5]);
        xlim([0 210]); ylim([-1 1]);
    end
%         saveas(f2,['E:\data\vHPC\all\figure\cluster\dmean_clusterbypn',num2str(ii),'.tif']);
% close all
    
end

%% collecting pn trials
load('behavior_align.mat','C_raw','set_eventframe','fr','trial_80P_ind','trial_80Pno_ind','trial_80Pnx_ind')
C_event = C_raw(:,set_eventframe(:));
C_event = reshape(C_event', [],12*fr,size(C_event,1));
C_pncueset = C_event(trial_80P_ind,:,:);
C_pnset = C_event(trial_80Pno_ind,:,:);
C_nopnset = C_event(trial_80Pnx_ind,:,:);
C_pnbasem = mean(mean(C_pncueset(:,4.5*fr+1:5.5*fr,:),2));
C_pnbasestd = std(mean(C_pncueset(:,4.5*fr+1:5.5*fr,:),2));
C_pnz = (C_pnset-C_pnbasem)./C_pnbasestd;
C_nopnz = (C_nopnset-C_pnbasem)./C_pnbasestd;
for oo = 1:size(C_pnz,3)
[~,p1] = ttest(mean(C_pnz(:,4.5*fr+1:5.5*fr,oo),2),mean(C_pnz(:,5.5*fr+1:6.5*fr,oo),2));
[~,p2] = ttest2(mean(C_pnz(:,5.5*fr+1:6.5*fr,oo),2),mean(C_nopnz(:,5.5*fr+1:6.5*fr,oo),2));
pnsigcell(oo,:) = [p1,p2];
end
save('pnsigcell.mat','C_pncueset','C_pnset','C_nopnset','C_pnz','C_nopnz','pnsigcell')


%% decoding period - plot
dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
    'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
fileID = 'E:\data\vHPC\all\decoding_period\';
fileID_base = 'E:\data\vHPC\all\decoding_period\';

cmap1=([0 56 66; 57 197 187; 55 0 102; 179 143 177;]./255);
for itask = 1:2
    if itask==1; load([fileID,'decoding_rwprob.mat']); load([fileID_base,'decoding_base_rwprob.mat']); tasktitle = 'rwprob';
    elseif itask==2; load([fileID,'decoding_rwpn.mat']); load([fileID_base,'decoding_base_rwpn.mat']); tasktitle = 'rwpn';
    end        
    for ii = 1:12 %mean decoding performance
        eval(['dec_tmp = decodingresult_period_',dirn{ii},';']);
        eval(['dec_tmp_base = decodingresult_base_period_',dirn{ii},';']);
        cue_performance= mean(mean([dec_tmp(1,1,:,:)==1;dec_tmp(2,1,:,:)==2;dec_tmp(3,1,:,:)==3],1),3);
        delay_performance= mean(mean([dec_tmp(1,2,:,:)==1;dec_tmp(2,2,:,:)==2;dec_tmp(3,2,:,:)==3],1),3);
        cue_base_performance= mean(mean([dec_tmp_base(1,1,:,:)==1;dec_tmp_base(2,1,:,:)==2;dec_tmp_base(3,1,:,:)==3],1),3);
        delay_base_performance= mean(mean([dec_tmp_base(1,2,:,:)==1;dec_tmp_base(2,2,:,:)==2;dec_tmp_base(3,2,:,:)==3],1),3);
        cue_performance_set(:,ii) = squeeze(cue_performance);
        delay_performance_set(:,ii) = squeeze(delay_performance);
        cue_base_performance_set(:,ii) = squeeze(cue_base_performance);
        delay_base_performance_set(:,ii) = squeeze(delay_base_performance);
        % cue, period, leave one out,mouse, task
        lbefore(:,:,:,ii,itask) = [dec_tmp_base(1,:,:,:)==1;dec_tmp_base(2,:,:,:)==2;dec_tmp_base(3,:,:,:)==3];
        lafter(:,:,:,ii,itask) = [dec_tmp(1,:,:,:)==1;dec_tmp(2,:,:,:)==2;dec_tmp(3,:,:,:)==3];
        pcuebefore(:,:,:,ii,itask) = [dec_tmp_base];
        pcueafter(:,:,:,ii,itask) = [dec_tmp];
        
    end
    
    [~,p1,~,stat1] = ttest(cue_performance_set(:,[1:6]),1/3,'tail','right');
    [~,p2,~,stat2] = ttest(cue_performance_set(:,[7:12]),1/3,'tail','right');
    [~,p3,~,stat3] = ttest(delay_performance_set(:,[1:6]),1/3,'tail','right');
    [~,p4,~,stat4] = ttest(delay_performance_set(:,[7:12]),1/3,'tail','right');
    [~,p5,~,stat5] = ttest(cue_base_performance_set(:,[1:6]),1/3,'tail','right');
    [~,p6,~,stat6] = ttest(cue_base_performance_set(:,[7:12]),1/3,'tail','right');
    [~,p7,~,stat7] = ttest(delay_base_performance_set(:,[1:6]),1/3,'tail','right');
    [~,p8,~,stat8] = ttest(delay_base_performance_set(:,[7:12]),1/3,'tail','right');
    pset(:,itask)= [p1,p2,p3,p4,p5,p6,p7,p8];
    statset(:,itask) = [stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8];
    for iperiod = 1:2
        if iperiod ==1;
            perf_tmp= cue_performance_set; base_tmp = cue_base_performance_set; tasktitle2 = '_cue';
        elseif iperiod ==2;
            perf_tmp= delay_performance_set; base_tmp = delay_base_performance_set; tasktitle2 = '_delay';end
        f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 2.5]);
        hold on
        for ibar = 1:2
            bar(0.8*(ibar-1)+0.1, mean(base_tmp(:,6*(ibar-1)+[1:6]),2),'barwidth',0.7,'FaceColor','none','Edgecolor',cmap1(ibar+2*(itask-1),:))
            errorbar(0.8*(ibar-1)+0.1, mean(base_tmp(:,6*(ibar-1)+[1:6]),2),sem(base_tmp(:,6*(ibar-1)+[1:6])),'color','k','capsize',3)
            %         plot([1,1]*2*(ibar-1)+0.1, [-0.5 0.5]*sem(mean(cue_performance_set(:,6*(ibar-1)+[1:6]),2))+mean(mean(cue_performance_set(:,6*(ibar-1)+[1:6]),2),1),'color','k')
            bar(0.8*(ibar-1)+2.1, mean(perf_tmp(:,6*(ibar-1)+[1:6]),2),'barwidth',0.7,'FaceColor',cmap1(ibar+2*(itask-1),:),'Edgecolor','none')
            errorbar(0.8*(ibar-1)+2.1, mean(perf_tmp(:,6*(ibar-1)+[1:6]),2),sem(perf_tmp(:,6*(ibar-1)+[1:6])),'color','k','capsize',3)
            %         plot([1,1]*2*(ibar-1)+0.9, [-0.5 0.5]*sem(mean(delay_performance_set(:,6*(ibar-1)+[1:6]),2))+mean(mean(delay_performance_set(:,6*(ibar-1)+[1:6]),2),1),'color','k')
        end
        line([-1 4], [1/3 1/3], 'LineStyle',':','color','k')
        xlim([-0.4 3.4]); xticks([0.1 0.9 2.1 2.9]); xticklabels({});
        ylim([0.25 0.75]); yticks([0.25 0.5 0.75]); yticklabels({});
%         print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\decoding\',tasktitle,tasktitle2,'.ai']);
        
%         anova_t= [base_tmp(:,1:6),perf_tmp(:,1:6),base_tmp(:,7:12),perf_tmp(:,7:12)];
%         mouse_anova = [ones(1,12),2*ones(1,12)];
%         learning_anova = [ones(1,6),2*ones(1,6),ones(1,6),2*ones(1,6)];
%         [c]=multcompare(rm,'Drug','by','Mouse');
%         [~,p,stat] = anovan(anova_t,{mouse_anova, learning_anova},'model','interaction','varnames',{'mouse','learning'});
%         pset{:,2*(itask-1)+iperiod} = p; stats{2*(itask-1)+iperiod} = stat;
%         eval(['cue_',tasktitle,'=cue_performance_set;'])
%         eval(['delay_',tasktitle,'=delay_performance_set;'])

        decoding_anova= [base_tmp',perf_tmp'];
        mouseG = [ones(1,6),2*ones(1,6)]';
        [tbl,rm] = simple_mixed_anova(decoding_anova, mouseG, {'learning'},{'Mouse'});
        tblset{:,2*(itask-1)+iperiod} = tbl;
    end
end

%% each cue
cmap1=([255 255 255; 255 255 255; 255 255 255; 255 255 255;...
    0 56 66; 57 197 187; 55 0 102; 179 143 177; ]./255);
totall =cat(6,lbefore,lafter);
for itask= 1:2
    figure
    for iperiod = 1:2 % cue or delay1
        for ilearning = 1:2 % before or after learning
            subplot(1,4,2*(iperiod-1) + ilearning)
            hold on
            for icue = 1:3 % cue 1,2,3
                settmpd = mean(totall(icue,iperiod,:,1:6,itask,ilearning));
                bar(icue, mean(settmpd),'facecolor',cmap1(2*(itask-1)+4*(ilearning-1)+1,:),'edgecolor',cmap1(2*(itask-1)+5,:))
                errorbar(icue, mean(settmpd),sem(settmpd),'color','k')
                settmpv = mean(totall(icue,iperiod,:,7:12,itask,ilearning));
                bar(icue+4, mean(settmpv),'facecolor',cmap1(2*(itask-1)+4*(ilearning-1)+2,:),'edgecolor',cmap1(2*(itask-1)+6,:))
                errorbar(icue+4, mean(settmpv),sem(settmpv),'color','k')
            end
            
            xlim([0 8]); xticks([]); ylim([0 1]);
        end
    end
end
%% confusion matrix
% cue, period, leave one out,mouse, task
truelabel = repmat( [1 ; 2 ; 3],1,1,10,6);
retrue = reshape(truelabel, 1,[]);
groupname = {'dHP1','vHP1','dHP2','vHP2'};
for itask= 1:2
    for iperiod = 1:2 % cue or delay1
        dHP1=pcuebefore(:,iperiod, :, 1:6,itask);
        vHP1=pcuebefore(:,iperiod, :, 7:12,itask);
        dHP2=pcueafter(:,iperiod, :, 1:6,itask);        
        vHP2=pcueafter(:,iperiod, :, 7:12,itask);
        for ii = 1:4
            eval(['delick = ', groupname{1,ii},';'])
            redelick = reshape(delick,1,[]);
            [m,order] = confusionmat(retrue, redelick);
            m = flip(m,2);
            mset(:,:,ii,iperiod,itask) = m;
        end
    end
end

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 18 18]);
% cmap_pink = [linspace(255, 0,300)', linspace(255, 0,300)', linspace(255, 0,300)']./255;
% cmap_pink = [linspace(255, 153,300)', linspace(255, 12,300)', linspace(255, 88,300)']./255; %(255 20 147)
% cmap_pink = [linspace(255, 0,300)', linspace(255, 56,300)', linspace(255, 66,300)']./255; % rwprob
cmap_pink = [linspace(255, 55,300)', linspace(255, 0,300)', linspace(255, 102,300)']./255; %rwpn
for itask= 1:2
    for iperiod = 1:2 % cue or delay1
         for ii = 1:4
             hold on
             subplot(4,4,ii+4*(iperiod-1)+8*(itask-1))
             imagesc(mset(:,:,ii,iperiod,itask)./60)
             caxis([0.05 0.6])
             colormap(cmap_pink)
             caxis([0 0.7])
%              colorbar
%              xlabel('Predicted cue')
%              ylabel('Actual cue')
         end
    end
end
saveas(f1,['E:\data\vHPC\all\figure\decoding\confusionmatrix2.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\decoding\confusionmatrix.ai']);
% cmsize = size(decodingresult_lick);
% truelabel = cat(1,ones([1,cmsize(2:4)]),2*ones([1,cmsize(2:4)]),3*ones([1,cmsize(2:4)]));
% retrue = reshape(truelabel, 1,[]);
% redelick = reshape(decodingresult_lick, 1,[]);
% [m,order] = confusionmat(retrue, redelick);


%% decoding anova
timer= [ones(200,1);2*ones(200,1)];
HPG= [ones(100,1);2*ones(100,1);ones(100,1);2*ones(100,1)];
cue_rwpn_d = mean(cue_rwpn(:,1:6),2);
cue_rwpn_v = mean(cue_rwpn(:,7:12),2);
delay_rwpn_d = mean(delay_rwpn(:,1:6),2);
delay_rwpn_v = mean(delay_rwpn(:,7:12),2);
[p1,~,stat1] = anovan([cue_rwpn_d;cue_rwpn_v;delay_rwpn_d;delay_rwpn_v],{timer,HPG},'full')
%%
cue_rwprob_d = mean(cue_rwprob(:,1:6),2);
cue_rwprob_v = mean(cue_rwprob(:,7:12),2);
delay_rwprob_d = mean(delay_rwprob(:,1:6),2);
delay_rwprob_v = mean(delay_rwprob(:,7:12),2);
[p1,~,stat1] = anovan([cue_rwprob_d;cue_rwprob_v;delay_rwprob_d;delay_rwprob_v],{timer,HPG},'full')

%% decoding_0 across session
clear
dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
    'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
fileID = 'E:\data\vHPC\all\decoding_0\';
cmap1=([0 56 66; 57 197 187; 55 0 102; 179 143 177;]./255);
cmap2 = ([255 255 255; 55 0 102; 179 143 177;]./255);
cd(fileID)
% cue(3), period(cue,delay,2),mintrial, ntrain(4), ntest(4)

for ii = 1:12 %mean decoding performance
    load([dirn{ii},'_rwpn_0_decoding_period_samecell.mat'])
    decoding_rawset{ii} = decodingresult_period; 
    dec_tmp = decodingresult_period;
    cue_performance(:,:,ii)= squeeze(mean(mean([dec_tmp(1,1,:,:,:)==1;dec_tmp(2,1,:,:,:)==2;dec_tmp(3,1,:,:,:)==3],1),3));
    delay_performance(:,:,ii)= squeeze(mean(mean([dec_tmp(1,2,:,:,:)==1;dec_tmp(2,2,:,:,:)==2;dec_tmp(3,2,:,:,:)==3],1),3));
end

for iperiod = 1:2
    if iperiod ==1;
        perf_tmp= cue_performance;  tasktitle2 = '_cue';
    elseif iperiod ==2;
        perf_tmp= delay_performance; tasktitle2 = '_delay';end
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 8]);
    hold on
    for itrain = 1:4
        subplot(2,4,itrain)
        hold on
        for itest = 1:4
            bar(itest,mean(perf_tmp(itrain,itest,1:6),3),    'FaceColor',cmap2(1+(itrain==itest)*1,:),'Edgecolor',cmap1(3,:))
            errorbar(itest,mean(perf_tmp(itrain,itest,1:6),3), sem(perf_tmp(itrain,itest,1:6)),'color','k')
            [~,p]=ttest(perf_tmp(itrain,itest,1:6),1/3,'tail','right');
            pset(itest, itrain, iperiod) = p;
        end
        line([0 5], [0.3 0.3],'color',[0.5 0.5 0.5])
        xticks([]); xlim([0 5]); ylim([0 1]); yticks([0 1])
        
        subplot(2,4,itrain+4)
        hold on
        for itest = 1:4
            bar(itest,mean(perf_tmp(itrain,itest,7:12),3) ,'FaceColor',cmap2(1+(itrain==itest)*2,:),'Edgecolor',cmap1(4,:))
            errorbar(itest,mean(perf_tmp(itrain,itest,7:12),3) ,sem(perf_tmp(itrain,itest,7:12)),'color','k')
            [~,p]=ttest(perf_tmp(itrain,itest,7:12),1/3,'tail','right');
            pset(itest, itrain+4, iperiod) = p;
        end
        line([0 5], [0.3 0.3],'color',[0.5 0.5 0.5])
        xticks([]); xlim([0 5]); ylim([0 1]); yticks([0 1])
        anova_tmp = squeeze(perf_tmp(itrain,:,:))'; %anova_tmp = reshape(anova_tmp,1,[]);
        mouse_tmp = [zeros(1,6), ones(1,6)]'; %mouse_tmp = reshape(mouse_tmp,1,[]);
        [tbl,rm] = simple_mixed_anova(anova_tmp,mouse_tmp, {'Test'},{'Mouse'});
        [c_test]=multcompare(rm,'Test','ComparisonType','tukey-kramer');%cue
        tblset{iperiod,itrain} = tbl;
        c_testset{iperiod,itrain} = c_test;

    end
%     print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\decoding\across_session_decoding',tasktitle2,'.ai']);
    
%     for ibar = 1:2
%         bar(0.8*(ibar-1)+0.1, mean(base_tmp(:,6*(ibar-1)+[1:6]),2),'barwidth',0.7,'FaceColor','none','Edgecolor',cmap1(ibar+2*(itask-1),:))
%         errorbar(0.8*(ibar-1)+0.1, mean(base_tmp(:,6*(ibar-1)+[1:6]),2),sem(base_tmp(:,6*(ibar-1)+[1:6])),'color','k','capsize',3)
%         %         plot([1,1]*2*(ibar-1)+0.1, [-0.5 0.5]*sem(mean(cue_performance_set(:,6*(ibar-1)+[1:6]),2))+mean(mean(cue_performance_set(:,6*(ibar-1)+[1:6]),2),1),'color','k')
%         bar(0.8*(ibar-1)+2.1, mean(perf_tmp(:,6*(ibar-1)+[1:6]),2),'barwidth',0.7,'FaceColor',cmap1(ibar+2*(itask-1),:),'Edgecolor','none')
%         errorbar(0.8*(ibar-1)+2.1, mean(perf_tmp(:,6*(ibar-1)+[1:6]),2),sem(perf_tmp(:,6*(ibar-1)+[1:6])),'color','k','capsize',3)
%         %         plot([1,1]*2*(ibar-1)+0.9, [-0.5 0.5]*sem(mean(delay_performance_set(:,6*(ibar-1)+[1:6]),2))+mean(mean(delay_performance_set(:,6*(ibar-1)+[1:6]),2),1),'color','k')
%     end
%     line([-1 4], [1/3 1/3], 'LineStyle',':','color','k')
%     xlim([-0.4 3.4]); xticks([0.1 0.9 2.1 2.9]); xticklabels({});
%     ylim([0.25 0.75]); yticks([0.25 0.5 0.75]); yticklabels({});
%     %         print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\decoding\',tasktitle,tasktitle2,'.ai']);
%     
%     %         anova_t= [base_tmp(:,1:6),perf_tmp(:,1:6),base_tmp(:,7:12),perf_tmp(:,7:12)];
%     %         mouse_anova = [ones(1,12),2*ones(1,12)];
%     %         learning_anova = [ones(1,6),2*ones(1,6),ones(1,6),2*ones(1,6)];
%     %         [c]=multcompare(rm,'Drug','by','Mouse');
%     %         [~,p,stat] = anovan(anova_t,{mouse_anova, learning_anova},'model','interaction','varnames',{'mouse','learning'});
%     %         pset{:,2*(itask-1)+iperiod} = p; stats{2*(itask-1)+iperiod} = stat;
%     %         eval(['cue_',tasktitle,'=cue_performance_set;'])
%     %         eval(['delay_',tasktitle,'=delay_performance_set;'])
%     
%     decoding_anova= [base_tmp',perf_tmp'];
%     mouseG = [ones(1,6),2*ones(1,6)]';
%     [tbl,rm] = simple_mixed_anova(decoding_anova, mouseG, {'learning'},{'Mouse'});
%     tblset{:,2*(itask-1)+iperiod} = tbl;
end

%% number of cell of intersect
nx=cellfun(@(x) sum(x~=0), cellindexset','Uniformoutput',false);
nsumx_dHP = sum(cell2mat(nx(1:6)));
nsumx_vHP = sum(cell2mat(nx(7:12)));
ny = cellfun(@(y) sum(y(:,1)~=0 & y(:,2)~=0),cellindexset','Uniformoutput',false);
nsumy_dHP = sum(cell2mat(ny(1:6)));
nsumy_vHP = sum(cell2mat(ny(7:12)));
%% speed cell
load('behavior_align.mat','C_z_all','cylinderspeed','timeline','fr')

lickinfo = timeline{1,1};
licksort = sortrows(lickinfo(:,1:3), 2);
for iii = 1:licksort(end,2);
    licktemp = licksort((licksort(:,2)==iii),1);
    lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
end
lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms
lick= reshape(lick_reg(:,41:80)',420*40,[]);
cs = reshape(cylinderspeed(:,1:120)',420*120,[]);

ctmp = movmean(C_z_all,3,2);
craw_lick = permute(ctmp(:,4*fr+1:3:8*fr,:),[2 1 3]);
craw_lick = reshape(craw_lick,size(craw_lick,3),[420*40]);
craw_lick = craw_lick';

craw_speed = permute(ctmp(:,1:3:end,:),[2 1 3]);
craw_speed = reshape(craw_speed,size(craw_speed,3),[420*120]);
craw_speed = craw_speed';

clear rset pset
for oo = 1:size(craw_speed,2)
    [r,p] = corrcoef(craw_speed(:,oo),cs);
    rset(oo,1) = r(1,2);
    pset(oo,1) = p(1,2);
    [r,p] = corrcoef(craw_lick(:,oo),lick);
    rset(oo,2) = r(1,2);
    pset(oo,2) = p(1,2);
end
eval(['rset_',mousename,'=rset;'])
eval(['pset_',mousename,'=pset;'])


rset_dHP = [rset_dHP03;rset_dHP04;rset_dHP06;rset_dHP07;rset_dHP08;rset_dHP10];
rset_vHP = [rset_vHP06;rset_vHP07; rset_vHP08; rset_vHP11; rset_vHP12; rset_vHP14];
xbin = {[-0.02 : 0.005 : 0.04], [-0.6 : 0.05 : 0.6]};
tl = {'Speed','Lick'};
for ii = 1:2
subplot(1,2,ii)
hold on
histogram(rset_dHP(:,ii),xbin{ii})
histogram(rset_vHP(:,ii),xbin{ii})
xlabel('r'); ylabel('Neuron')
title(tl{ii})
end
legend('dHP','vHP')
%% check 25% activity

load('behavior_align.mat','C_z_all','odorCue','fr', 'trial_80R_ind', 'trial_50R_ind', 'trial_20R_ind')
load('fon_period01_220228.mat','pvalueset')
iterm=4+1; %value
iperiod=7; % delay period
signeuron = find(pvalueset(iterm, iperiod, :)<0.05);
C_z_value = squeeze(mean(C_z_all(:,2*fr+1:3.5*fr,signeuron),2));
trialInd = 1*trial_80R_ind +2*trial_50R_ind+3*trial_20R_ind;


f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);%
for ineuron = 1:length(signeuron)
    clf
    hold on
    for ibar = 1:3
        scatter(ibar+0.5*rand(sum(trialInd==ibar),1), C_z_value(trialInd==ibar,ineuron))
        line([ibar ibar+0.5], mean(C_z_value(trialInd==ibar,ineuron))*[1 1] ,'color','k','linestyle','-','lineWidth',5)
    end   
    saveas(f1,[cd,'\figure\25 distribution\',num2str(signeuron(ineuron)),'.tif'])
end







