%% anticipatoty lick
clear
figure_folder = 'E:\data\vHPC\all\figure\behavior\';

dirlist= dir('E:\data\vHPC\all\behavior\rwprob');
cd('E:\data\vHPC\all\behavior\rwprob');
% dirlist= dir('E:\data\vHPC\all\behavior\rwpn');
% cd('E:\data\vHPC\all\behavior\rwpn');
% dirlist= dir('E:\data\vHPC\all\behavior\rwprobrev');
% cd('E:\data\vHPC\all\behavior\rwprobrev');
% dirlist= dir('E:\data\vHPC\all\behavior\rwpnrev');
% cd('E:\data\vHPC\all\behavior\rwpnrev');
dirname={dirlist(:).name};
dHPlist = dirname(startsWith(dirname,'dHP'));
vHPlist = dirname(startsWith(dirname,'vHP'));
%% anticipatory lick 
indexHP=[ones(1,size(dHPlist,2)),zeros(1,size(vHPlist,2))];
HP_trace = NaN(4,42,size(indexHP,2));
for ifile = 1:size(indexHP,2)
    if indexHP(1,ifile)==1; load(dHPlist{1,ifile});
    elseif indexHP(1,ifile)==0; load(vHPlist{1,ifile-sum(indexHP)}); end
    %     ntrial(ifile,:) = [sum(odorCue==0) sum(odorCue==1) sum(odorCue==2) sum(odorCue==3)];
    for ii = 1:4
        HP_trace(ii,:,ifile) = histcounts(timeline{1,1}(timeline{1,1}(:,3)==ii,1),[4000-100:100:8000+100]*1000);
    end
    odorCueset(:,ifile) = odorCue;
    lickNumset(:,ifile) = lickNum;
    cueKindset(:,ifile) = lickEachCue.cue_kind;
    if ismember(lickEachCue.cue_kind,[0 1 3 4])
        meanr= [mean(lickEachCue.r80), mean(lickEachCue.p80),mean(lickEachCue.r20)];
        meanz = mean([lickEachCue.r80;lickEachCue.p80;lickEachCue.r20]);
    stdz = std([lickEachCue.r80;lickEachCue.p80;lickEachCue.r20]);
        colormap = [0 25 145; 230 0 1; 128 128 128]./255;
        lickz = ([mean(lickEachCue.r80-meanz),mean(lickEachCue.p80-meanz),mean(lickEachCue.r20-meanz)]./stdz);
        behaviortitle = 'behavior_rwpn';
    elseif ismember(lickEachCue.cue_kind,[0 1 2 3])
        meanr= [mean(lickEachCue.r80), mean(lickEachCue.r50),mean(lickEachCue.r20)];
        meanz = mean([lickEachCue.r80;lickEachCue.r50;lickEachCue.r20]);
        stdz = std([lickEachCue.r80;lickEachCue.r50;lickEachCue.r20]);
        colormap = [0 25 145; 0 153 255; 128 128 128]./255;
        lickz = ([mean(lickEachCue.r80-meanz),mean(lickEachCue.r50-meanz),mean(lickEachCue.r20-meanz)]./stdz);
        behaviortitle = 'behavior_rwprob';
    end
    HP_lick(ifile,:,1)=meanr./1.5; HP_lick(ifile,:,2)=lickz;
end
HP_trace = movmean(HP_trace,5,2)/10;
[tbl,rm] = simple_mixed_anova(HP_lick(:,:,1), indexHP', {'Cue'},{'Mouse'});
[c]=multcompare(rm,'Cue','ComparisonType','tukey-kramer');
%% Anticipatory trace lick plot
day = datetime('today');
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
cmap1= [0 25 145; 0 153 255; 128 128 128;230 0 1]./255;

hold on
stdshade(permute(HP_trace(3,:,:),[3 2 1]),0.3, cmap1(3,:))
stdshade(permute(HP_trace(1,:,:),[3 2 1]),0.3, cmap1(1,:))
if ismember(lickEachCue.cue_kind,[0 1 3 4])
stdshade(permute(HP_trace(4,:,:),[3 2 1]),0.3, cmap1(4,:))
elseif ismember(lickEachCue.cue_kind,[0 1 2 3])
stdshade(permute(HP_trace(2,:,:),[3 2 1]),0.3, cmap1(2,:))
end
xlim([1 42]); xticks([2 1+15 1+40]); xticklabels([])
line([2 2], [0 10], 'linestyle',':','color','k')
line([16 16], [0 10], 'linestyle',':','color','k')
line([41 41], [0 10], 'linestyle',':','color','k')
ylim([0 9]); yticks([0 9]); yticklabels([]);

% saveas(f1,['E:\data\vHPC\all\figure\behavior\',behaviortitle,d,'.tif'])
print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\behavior\trace_',behaviortitle,d,'.ai']);

% cue=[ones(1,size(indexHP,2)),2*ones(1,size(indexHP,2)),zeros(1,size(indexHP,2))];
% [p,~,stats] = anovan(reshape(HP_lick(:,:,2),1,[]),{cue,repmat(indexHP,1,3)},...
%     'model','linear','varnames',{'Cue','HP'});
% s = multcompare(stats,'Dimension',1);
%% Anticipatory mean lick plot
day = datetime('today');
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.5 3]);
ii=1; % for raw lick
hold on
% errorbar(1,mean(HP_lick(:,1,ii)),sem(HP_lick(:,1,ii)),'k','CapSize',3)
% errorbar(2,mean(HP_lick(:,2,ii)),sem(HP_lick(:,2,ii)),'k','CapSize',3)
% errorbar(3,mean(HP_lick(:,3,ii)),sem(HP_lick(:,3,ii)),'k','CapSize',3)
% scatter(1,mean(HP_lick(:,1,ii)),15,colormap(1,:),'s','filled')
% scatter(2,mean(HP_lick(:,2,ii)),15,colormap(2,:),'s','filled')
% scatter(3,mean(HP_lick(:,3,ii)),15,colormap(3,:),'s','filled')
bar(1,mean(HP_lick(:,1,ii)),0.6,'Facecolor',colormap(1,:),'Edgecolor','none')
bar(2,mean(HP_lick(:,2,ii)),0.6,'Facecolor',colormap(2,:),'Edgecolor','none')
bar(3,mean(HP_lick(:,3,ii)),0.6,'Facecolor',colormap(3,:),'Edgecolor','none')
errorbar(1,mean(HP_lick(:,1,ii)),sem(HP_lick(:,1,ii)),'k')
errorbar(2,mean(HP_lick(:,2,ii)),sem(HP_lick(:,2,ii)),'k')
errorbar(3,mean(HP_lick(:,3,ii)),sem(HP_lick(:,3,ii)),'k')
scatter(0.8+0.1*rand(size(dHPlist,2),1),HP_lick(indexHP==1,1,ii),3,'k','o','MarkerFaceColor',[0 0 0])
scatter(1.8+0.1*rand(size(dHPlist,2),1),HP_lick(indexHP==1,2,ii),3,'k','o','MarkerFaceColor',[0 0 0])
scatter(2.8+0.1*rand(size(dHPlist,2),1),HP_lick(indexHP==1,3,ii),3,'k','o','MarkerFaceColor',[0 0 0])
scatter(1.1+0.1*rand(size(vHPlist,2),1),HP_lick(indexHP==0,1,ii),3,'k','o','MarkerFaceColor',[1 1 1])
scatter(2.1+0.1*rand(size(vHPlist,2),1),HP_lick(indexHP==0,2,ii),3,'k','o','MarkerFaceColor',[1 1 1])
scatter(3.1+0.1*rand(size(vHPlist,2),1),HP_lick(indexHP==0,3,ii),3,'k','o','MarkerFaceColor',[1 1 1])
xlim([0.5 3.5]); xticks([1 2 3]); xticklabels([])
ylim([0 9]); yticks([0 9]); yticklabels([]);


% saveas(f1,['E:\data\vHPC\all\figure\behavior\',behaviortitle,d,'.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\behavior\',behaviortitle,d,'.ai']);

% cue=[ones(1,size(indexHP,2)),2*ones(1,size(indexHP,2)),zeros(1,size(indexHP,2))];
% [p,~,stats] = anovan(reshape(HP_lick(:,:,1),1,[]),{cue,repmat(indexHP,1,3)},...
%     'model','full','varnames',{'Cue','HP'});
% s = multcompare(stats,'Dimension',1);
[tbl,rm] = simple_mixed_anova(HP_lick(:,:,1),dHPlist, {'Cue'},{'Mouse'});
[c_cue]=multcompare(rm,'Cue','ComparisonType','tukey-kramer');%cue

%% Reversal lick plot
day = datetime('today');
d = datestr(day,'yymmdd');
cmap1= [101 110 236; 150 200 255; 200 200 200; 255 130 150]./255;
cmap2= [30 45 232; 120 177 255; 125 125 125; 255 30 70]./255;
rev_end = 200;
indexHP=[ones(1,size(dHPlist,2)),zeros(1,size(vHPlist,2))];
HP_revlickmean = NaN(size(indexHP,2),6,4);


for ifile = 1:size(indexHP,2)
    if indexHP(1,ifile)==1; load(dHPlist{1,ifile});
    elseif indexHP(1,ifile)==0; load(vHPlist{1,ifile-sum(indexHP)}); end
    revstart = rev_index(1,5);
    rev_range = revstart-100:50:revstart+rev_end;
    if revstart+rev_end>420; rev_range(end) = 420; end
    for ii = 1:rev_end/50+2
        licktmp = lickNum(rev_range(ii):rev_range(ii+1));
        odortmp = odorCue(rev_range(ii):rev_range(ii+1));
        HP_revlickmean(ifile, ii, 1) = mean(licktmp(odortmp==(find(lickEachCue.cue_kind==1)-1)))./1.5;
        HP_revlickmean(ifile, ii, 2) = mean(licktmp(odortmp==(find(lickEachCue.cue_kind==2)-1)))./1.5;
        HP_revlickmean(ifile, ii, 3) = mean(licktmp(odortmp==(find(lickEachCue.cue_kind==3)-1)))./1.5;
        HP_revlickmean(ifile, ii, 4) = mean(licktmp(odortmp==(find(lickEachCue.cue_kind==4)-1)))./1.5;
    end
end
if ismember(lickEachCue.cue_kind,[0 1 3 4])
    behaviortitle = 'rwpnrev_';
    yset = [0 6];
    hset = zeros(1,6);
    for ii = 1:6
        [p,~,stat] = anova1(HP_revlickmean(:,ii,:),[],'off');
        pset(1,ii) = p;
        if p<0.05 & ii<2.5
            s = multcompare(stat,'display','off');
            if s(1,6)<0.1 & s(2,6)<0.1 & s(1,4)>0 & s(2,4)>0
                hset(1,ii) = 1;
            end
        elseif p<0.05 & ii>2.5
            s = multcompare(stat,'display','off');
            if s(2,6)<0.1 & s(3,6)<0.1 & s(2,4)<0 & s(3,4)<0
                hset(1,ii) = 1;
            end
        end
    end
elseif ismember(lickEachCue.cue_kind,[0 1 2 3])
    behaviortitle = 'rwprobrev_';
    yset = [0 5];
    hset = zeros(1,6);
    for ii = 1:6
        [p,~,stat] = anova1(HP_revlickmean(:,ii,:),[],'off');
        pset(1,ii) = p;
        if p<0.05 & ii<2.5
            s = multcompare(stat,'display','off');
            if s(:,6)<0.1 & s(1,4)>0 & s(2,4)>0 & s(3,4)>0
                hset(1,ii) = 1;
            end
        elseif p<0.05 & ii>2.5
            s = multcompare(stat,'display','off');
            if s(:,6)<0.1 & s(1,4)<0 & s(2,4)>0 & s(3,4)>0
                hset(1,ii) = 1;
            end
        end
    end
end
%%
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.5 2]);
for imouse = 1:size(indexHP,2)
hold on
plot(HP_revlickmean(imouse,:,1),'color',cmap1(1,:),'linewidth',0.25)
plot(HP_revlickmean(imouse,:,2),'color',cmap1(2,:),'linewidth',0.25)
plot(HP_revlickmean(imouse,:,3),'color',cmap1(3,:),'linewidth',0.25)
plot(HP_revlickmean(imouse,:,4),'color',cmap1(4,:),'linewidth',0.25)
end
errorbar(mean(HP_revlickmean(:,:,1)),sem(HP_revlickmean(:,:,1)),'color',cmap2(1,:),'linewidth',1)
errorbar(mean(HP_revlickmean(:,:,2)),sem(HP_revlickmean(:,:,2)),'color',cmap2(2,:),'linewidth',1)
errorbar(mean(HP_revlickmean(:,:,3)),sem(HP_revlickmean(:,:,3)),'color',cmap2(3,:),'linewidth',1)
errorbar(mean(HP_revlickmean(:,:,4)),sem(HP_revlickmean(:,:,4)),'color',cmap2(4,:),'linewidth',1)
line([2.95 2.95],[0 10],'linestyle',':','color','k')
xlim([0.5 6.5]); xticks([1:6]); xticklabels({num2str([-100:50:150]')});
% plot(yset(1,2)*ones(1,6)-0.3,'MarkerIndices',find(hset),'Linestyle','none','Marker','*','color','k','markersize',3)
ylim(yset); yticks(yset);
set(gca,'Fontsize',7)
% saveas(f1,['E:\data\vHPC\all\figure\behavior\',behaviortitle,d,'.tif'])
print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\behavior\trace_',behaviortitle,d,'.ai']);


%% cue, reward history effect
clear bset
dirn = dir('*.mat');

for jj = 1:length(dirn)
    a = dirn(jj).name;
    load(a)
    
    lickZ = (lickNum-mean(lickNum))./std(lickNum);
    RwCue = odorCue==(find(lickEachCue.cue_kind==1)-1);
    PnCue = odorCue==(find(lickEachCue.cue_kind==4)-1);
    value = 0.75*(odorCue==(find(lickEachCue.cue_kind==1)-1)) + 0.25*(odorCue==(find(lickEachCue.cue_kind==2)-1));
    Rw = RwCue.*waterReward;
    Pn = PnCue.*waterReward;
    
    y = lickZ(11:420);
    x1=[];x2=[];x3=[];x4=[];
    for ii = 1:11
        if ismember(lickEachCue.cue_kind,[0 1 3 4])
            x1= [x1,RwCue(ii:ii+409)];
            x2= [x2,PnCue(ii:ii+409)];
            x3= [x3,Rw(ii:ii+409)];
            x4= [x4,Pn(ii:ii+409)];
        elseif ismember(lickEachCue.cue_kind,[0 1 2 3])
            x1= [x1,value(ii:ii+409)];
            x2= [x2,Rw(ii:ii+409)];
            
        end
    end
    x=[x1,x2];
    
    [b,dev,stats]=glmfit(x,y,'Normal','link','identity');
    bset(:,jj) = b;
end


[h,p] = ttest2(bset,zeros(size(bset)),'Dim',2);

%% draw history effect-rwpn
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 9]);
subplot(1,3,1)
hold on
%plot(mean(bset(2:12,:),2))
errorbar(mean(bset(2:12,:),2), sem(bset(2:12,:)'),'k')
scatter(1:11, ones(1,11),36*h(2:12)+1e-30,'k','*')
line([1 11],[0 0],'color','k','linestyle',':')
ylim([-0.4 1.4]); xlim([1 11]); xticks([1,6,10,11]); xticklabels([-10 -5 -1 0]);
title('Rw cue effect')
subplot(1,3,2)
hold on
%plot(mean(bset(13:23,:),2))
errorbar(mean(bset(13:23,:),2), sem(bset(13:23,:)'),'k')
scatter(1:11, ones(1,11),50*h(13:23)+1e-30,'k','*')
line([1 11],[0 0],'color','k','linestyle',':')
ylim([-0.4 1.4]); xlim([1 11]); xticks([1,6,10,11]); xticklabels([-10 -5 -1 0]);
title('Pn cue effect')
xlabel('Past trials')
subplot(1,3,3)
hold on
%plot(mean(bset(24:34,:),2))
errorbar(mean(bset(24:34,:),2), sem(bset(24:34,:)'),'k')
scatter(1:11, ones(1,11),50*h(24:34)+1e-30,'k','*')
line([1 11],[0 0],'color','k','linestyle',':')
ylim([-0.4 1.4]); xlim([1 11]); xticks([1,6,10,11]); xticklabels([-10 -5 -1 0]);
title('Rw effect')
% saveas(f1,[cd,'\figure\',d,'history effect.tif'])
%% draw history effect-rwprob
day = datetime('today');
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 6]);
subplot(1,2,1)
hold on
%plot(mean(bset(2:12,:),2))
errorbar(mean(bset(2:12,:),2), sem(bset(2:12,:)'),'k')
scatter(1:11, ones(1,11),36*h(2:12)+1e-30,'k','*')
line([1 11],[0 0],'color','k','linestyle',':')
ylim([-0.4 1.4]); xlim([1 11]); xticks([1,6,10,11]); xticklabels([-10 -5 -1 0]);
title('Cue effect')
xlabel('Past trials')
subplot(1,2,2)
hold on
%plot(mean(bset(13:23,:),2))
errorbar(mean(bset(13:23,:),2), sem(bset(13:23,:)'),'k')
scatter(1:11, ones(1,11),50*h(13:23)+1e-30,'k','*')
line([1 11],[0 0],'color','k','linestyle',':')
ylim([-0.4 1.4]); xlim([1 11]); xticks([1,6,10,11]); xticklabels([-10 -5 -1 0]);
title('Rw effect')
xlabel('Past trials')

saveas(f1,['E:\data\vHPC\all\figure\behavior\','history effect',d,'.tif'])