clear
cd('E:\data\vHPC\Muscimol\rwprob'); exclusionlist ={'dHPc11'};
% cd('E:\data\vHPC\Muscimol\rwpn'); exclusionlist ={'vHPc10'; 'vHPc11'};
% cd('E:\data\vHPC\Muscimol\ctrl_rwprob'); exclusionlist ={'dHPc11'};
% cd('E:\data\vHPC\Muscimol\ctrl_rwpn'); exclusionlist ={'vHPc10'; 'vHPc11'};
dirn = dir('*.mat');
dirname = {dirn.name};
dirnameset = cell(length(dirname),7);
for ii= 1:length(dirname);
    dirnameset(ii,:)= strsplit(dirname{ii},'_');
end
for iexp = 1: length(exclusionlist);
    indexEXCEPT(:,iexp) = strcmp(exclusionlist{iexp},strcat(dirnameset(:,1),dirnameset(:,2)));
end
indexEXCEPT = sum(indexEXCEPT,2)>0;
dirname = dirname(~indexEXCEPT);
dirnameset = dirnameset(~indexEXCEPT,:);

dHPlist = strcmp('dHP',dirnameset(:,1));
vHPlist = strcmp('vHP',dirnameset(:,1));
Muslist  = strcmpi('mus',dirnameset(:,5));
Acsflist  = strcmpi('acsf',dirnameset(:,5));
Ctrllist = strcmpi('ctrl',dirnameset(:,5));
%1 : dHP acsf 2: dHP mus 3: vHP acsf 4: vHP mus
indexHP=1*(dHPlist&Acsflist)+2*(dHPlist&Muslist)+3*(vHPlist&Acsflist)+4*(vHPlist&Muslist)+...
    1*(Ctrllist&dHPlist)+2*(Ctrllist&dHPlist);
% anticipatory lick 
HP_trace = NaN(4,42,size(dirname,2));
for ifile = 1:size(dirname,2)
    load(dirname{1,ifile});
    for ii = 1:4
        HP_trace(ii,:,ifile) = histcounts(timeline{1,1}(timeline{1,1}(:,3)==ii,1),[4000-100:100:8000+100]*1000);
    end
    odorCueset(:,ifile) = odorCue;
    lickNumset(:,ifile) = lickNum;
    allLickNumset(:,ifile) = size(lickTime,1);
    allSpeedset(:,ifile) = size(cylinderTime,1);
    cueKindset(:,ifile) = lickEachCue.cue_kind;
    
    lick_bout_onset = [0; find(diff(lickTime(:,1))>1*10^6)];
    lick_bout_offset = [find(diff(lickTime(:,1))>1*10^6);size(lickTime(:,1),1)];
    lick_bout_offset([diff(lick_bout_onset)<2;false])=[];
    lick_bout_onset([diff(lick_bout_onset)<2;false])=[];    
    
    for ibout = 1:size(lick_bout_onset,1)-1
        lick_bout_Hz(ibout) = (lick_bout_offset(ibout)-lick_bout_onset(ibout))/...
            (lickTime(lick_bout_offset(ibout),1)-lickTime(lick_bout_onset(ibout)+1,1))*10^6;
    end
    lickboutset{:,ifile} = lick_bout_Hz; clear lick_bout_Hz;
    
    if ismember(lickEachCue.cue_kind,[0 1 3 4])
        meanr= [mean(lickEachCue.r80), mean(lickEachCue.p80),mean(lickEachCue.r20)];
        meanz = mean([lickEachCue.r80;lickEachCue.p80;lickEachCue.r20]);
        stdz = std([lickEachCue.r80;lickEachCue.p80;lickEachCue.r20]);
        colormap = [0 28 145; 230 0 1; 128 128 128]./255;
        lickz = ([mean(lickEachCue.r80-meanz),mean(lickEachCue.p80-meanz),mean(lickEachCue.r20-meanz)]./stdz);
        behaviortitle = 'rwpn';
    elseif ismember(lickEachCue.cue_kind,[0 1 2 3])
        meanr= [mean(lickEachCue.r80), mean(lickEachCue.r50),mean(lickEachCue.r20)];
        meanz = mean([lickEachCue.r80;lickEachCue.r50;lickEachCue.r20]);
        stdz = std([lickEachCue.r80;lickEachCue.r50;lickEachCue.r20]);
        colormap = [0 28 145; 0 153 255; 128 128 128]./255;
        lickz = ([mean(lickEachCue.r80-meanz),mean(lickEachCue.r50-meanz),mean(lickEachCue.r20-meanz)]./stdz);
        behaviortitle = 'rwprob';
    end
    HP_lick(ifile,:,1)=meanr; HP_lick(ifile,:,2)=lickz;
end

HP_trace = movmean(HP_trace,5,2)/10;
%% Anticipatory trace lick plot - trace
day = datetime('today');
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 3]);
% cmap1= [30 45 232; 120 177 255; 125 125 125; 255 30 70]./255;
cmap1= [0 0 0; 80 80 80; 150 150 150; 80 80 80;...
%     210 98 32; 244 164 68; 251 212 127; 244 164 68;]./255;
    230 0 1; 240 56 0; 255 140 0; 240 56 0 ]./255;
for iplot = 1:2
    subplot(1,2,iplot)
    hold on
    stdshade(permute(HP_trace(3,:,indexHP==2*iplot),[3 2 1]),0.3, cmap1(7,:),[])
    stdshade(permute(HP_trace(3,:,indexHP==2*iplot-1),[3 2 1]),0.3, cmap1(3,:))
    stdshade(permute(HP_trace(1,:,indexHP==2*iplot),[3 2 1]),0.3, cmap1(5,:),[])
    stdshade(permute(HP_trace(1,:,indexHP==2*iplot-1),[3 2 1]),0.3, cmap1(1,:))    
    if ismember(lickEachCue.cue_kind,[0 1 3 4])
        stdshade(permute(HP_trace(4,:,indexHP==2*iplot),[3 2 1]),0.3, cmap1(8,:),[])
        stdshade(permute(HP_trace(4,:,indexHP==2*iplot-1),[3 2 1]),0.3, cmap1(4,:))        
    elseif ismember(lickEachCue.cue_kind,[0 1 2 3])
        stdshade(permute(HP_trace(2,:,indexHP==2*iplot),[3 2 1]),0.3, cmap1(6,:),[])
        stdshade(permute(HP_trace(2,:,indexHP==2*iplot-1),[3 2 1]),0.3, cmap1(2,:))        
    end
    xlim([1 42]); xticks([2 1+15 1+40]); xticklabels([])
    line([2 2], [0 10], 'linestyle',':','color','k')
    line([16 16], [0 10], 'linestyle',':','color','k')
    line([41 41], [0 10], 'linestyle',':','color','k')
    ylim([0 10]); yticks([0 10]); yticklabels([0 10]);
end

% saveas(f1,['E:\data\vHPC\Muscimol\figure\',behaviortitle,'_totallick.tif'])
print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\',behaviortitle,'_totallick.ai']);

%% Anticipatory mean lick plot mean
% day = datetime('today');
% d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 2.5]);

% hold on
for ii=1:2
    subplot(1,2,ii)
    hold on
    for ibar = 1:3
        bar(ibar-0.25,mean(HP_lick(indexHP==2*ii-1,ibar,1)),0.4,'Facecolor',colormap(ibar,:),'Edgecolor','none')
        errorbar(ibar-0.25,mean(HP_lick(indexHP==2*ii-1,ibar,1)),sem(HP_lick(indexHP==2*ii-1,ibar,1)),'k')
        bar(ibar+0.25,mean(HP_lick(indexHP==2*ii,ibar,1)),0.4,'Facecolor','none','Edgecolor',colormap(ibar,:))
        errorbar(ibar+0.25,mean(HP_lick(indexHP==2*ii,ibar,1)),sem(HP_lick(indexHP==2*ii,ibar,1)),'k')
%         plot([ibar*ones(sum(indexHP==2*ii-1),1)-0.25 ibar*ones(sum(indexHP==2*ii-1),1)+0.25]',...
%             [HP_lick(indexHP==2*ii-1,ibar,1) HP_lick(indexHP==2*ii,ibar,1)]','color',[0.3 0.3 0.3])
    end

    xlim([0.3 3.7]); xticks([]); xticklabels([])
    ylim([0 10]); yticks([0 10]); yticklabels([0 10]);
    
end

% saveas(f1,['E:\data\vHPC\Muscimol\figure\',behaviortitle,'_meanlick.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\',behaviortitle,'_meanlick.ai']);
%% Ctrl graph
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6.5 3]);
indexHP = dHPlist;
cmap1= [0 28 145; 0 153 255; 128 128 128; 230 0 1]./255;
ii=1;
subplot(1,3,[1:2])
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
ylim([0 10]); yticks([0 10]); yticklabels([]);

% axes('Position',[.65 .7 .25 .2])
subplot(1,3,3)
hold on
for ibar = 1:3
    bar(ibar,mean(HP_lick(:,ibar,1)),0.6,'Facecolor',cmap1(ibar,:),'Edgecolor','none')
    errorbar(ibar,mean(HP_lick(:,ibar,1)),sem(HP_lick(:,ibar,1)),'k','Capsize',3)
    scatter(ibar-0.2+0.1*rand(sum(indexHP==1),1),HP_lick(indexHP==1,ibar,1),1,'k','o','MarkerFaceColor',[0 0 0],'linewidth',0.25)
    scatter(ibar+0.1+0.1*rand(sum(indexHP==0),1),HP_lick(indexHP==0,ibar,1),1,'k','o','MarkerFaceColor',[1 1 1],'linewidth',0.25)
%     scatter(1.8+0.1*rand(size(dHPlist,2),1),HP_lick(indexHP==1,2,1),3,'k','o','MarkerFaceColor',[0.3 0.3 0.3])
%     scatter(2.8+0.1*rand(size(dHPlist,2),1),HP_lick(indexHP==1,3,1),3,'k','o','MarkerFaceColor',[0.3 0.3 0.3])
%     scatter(1.1+0.1*rand(size(vHPlist,2),1),HP_lick(indexHP==0,1,1),3,'k','o','MarkerFaceColor',[1 1 1])
%     scatter(2.1+0.1*rand(size(vHPlist,2),1),HP_lick(indexHP==0,2,1),3,'k','o','MarkerFaceColor',[1 1 1])
%     scatter(3.1+0.1*rand(size(vHPlist,2),1),HP_lick(indexHP==0,3,1),3,'k','o','MarkerFaceColor',[1 1 1])

    %         plot([ibar*ones(sum(indexHP==2*ii-1),1)-0.25 ibar*ones(sum(indexHP==2*ii-1),1)+0.25]',...
    %             [HP_lick(indexHP==2*ii-1,ibar,1) HP_lick(indexHP==2*ii,ibar,1)]','color',[0.3 0.3 0.3])
end

xlim([0.3 3.7]); xticks([]); xticklabels([])
ylim([0 10]); yticks([0 10]); yticklabels([]);


% saveas(f1,['E:\data\vHPC\Muscimol\figure\',behaviortitle,'_meanlick.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\',behaviortitle,'_ctrl_meanlick.ai']);


[tbl,rm] = simple_mixed_anova(HP_lick(:,:,1),dHPlist, {'Cue'},{'Mouse'});
[c_cue]=multcompare(rm,'Cue','ComparisonType','tukey-kramer');%cue

%% test
% indexM= [ones(sum(indexHP==1),1);zeros(sum(indexHP==3),1)];
% [tbl,rm] = simple_mixed_anova(cat(3,HP_lick(ismember(indexHP,[1 3]),:,1),HP_lick(ismember(indexHP,[2 4]),:,1)), indexM, {'Cue','Drug'},{'Mouse'});
% [c_cue]=multcompare(rm,'Cue','ComparisonType','tukey-kramer');%cue
% [c_dbyc]=multcompare(rm,'Drug','by','Cue','ComparisonType','tukey-kramer');%cue
% [c_cbyd]=multcompare(rm,'Cue','by','Drug','ComparisonType','tukey-kramer');%cue

[tbl,rm] = simple_mixed_anova(cat(3,HP_lick(ismember(indexHP,[1]),:,1),HP_lick(ismember(indexHP,[2]),:,1)),[], {'Cue','Drug'},{});
[c_cue]=multcompare(rm,'Cue','ComparisonType','tukey-kramer');%cue
[c_dbyc]=multcompare(rm,'Drug','by','Cue','ComparisonType','tukey-kramer');%cue
[c_cbyd]=multcompare(rm,'Cue','by','Drug','ComparisonType','tukey-kramer');%cue
% 
% indexM= [ones(sum(indexHP==1),1);zeros(sum(indexHP==3),1)];
% [tbl,rm] = simple_mixed_anova(HP_lick(ismember(indexHP,[1 3]),:,1),indexM, {'Cue'},{'Mouse'});
% [c_cue]=multcompare(rm,'Cue','ComparisonType','tukey-kramer');%cue
% [c_dbyc]=multcompare(rm,'Drug','by','Cue','ComparisonType','tukey-kramer');%cue
% [c_cbyd]=multcompare(rm,'Cue','by','Drug','ComparisonType','tukey-kramer');%cue
%% All anticipatory lick mean plot
if strcmp(behaviortitle,'rwpn')
    cmap1 = [55 0 102; 179 143 177]./255;
    ylimlist = [0.25 1]; cl = 1/2;
elseif strcmp(behaviortitle,'rwprob')
    cmap1 = [0 56 66; 57 197 187]./255;
    ylimlist = [0.25 0.75]; cl = 1/3;
end
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 3]);
for ii= 1:2
    if ii==1; dHPind=1; vHPind= 3;
    elseif ii==2; dHPind=2; vHPind=4;
    end
    hold on
    dg = mean(lickNumset(:,indexHP==dHPind));
    vg = mean(lickNumset(:,indexHP==vHPind));
    bar(ii-0.23,mean(dg),0.4,'facecolor','none','edgecolor',cmap1(1,:))
    errorbar(ii-0.23,mean(dg),sem(dg),'color','k','capsize',3, 'linewidth',0.5)
    bar(ii+0.23,mean(vg),0.4,'facecolor','none','edgecolor',cmap1(2,:))
    errorbar(ii+0.23,mean(vg),sem(vg),'color','k','capsize',3, 'linewidth',0.5)
    xticks([]); xlim([0.4 2.6]); ylim([0 5]); yticks([0 5]); yticklabels([])
%     [~,p,~,stat] = ttest(mean(lickNumset(:,indexHP==acsfind)),mean(lickNumset(:,indexHP==musind)))

end
print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\all_lick_mean_',behaviortitle,'.ai']);
    meanlick = nanmean(lickNumset)';
    indexM= [ones(sum(indexHP==1),1);zeros(sum(indexHP==3),1)];
    [tbl,rm] = simple_mixed_anova(cat(2,meanlick(ismember(indexHP,[1 3])),meanlick(ismember(indexHP,[2 4]))), indexM, {'Drug'},{'Mouse'});
%% lick rate adjustment
lickNumtmp= lickNumset;
p = 0;
while p < 0.05
    maxv= max(max(lickNumtmp(:,ismember(indexHP,[1 3]))));
    maxnum = sum(sum(lickNumtmp(:,ismember(indexHP,[1 3]))==maxv));
    minv = min(min(lickNumtmp(:,ismember(indexHP,[2 4]))));
    minnum = sum(sum(lickNumtmp(:,ismember(indexHP,[2 4]))==minv));
    selnum= min(maxnum,minnum);
    maxi=find(lickNumtmp==maxv); lickNumtmp(maxi(randperm(selnum,1)))=NaN;
    mini=find(lickNumtmp==minv); lickNumtmp(mini(randperm(selnum,1)))=NaN;
    
    meanlick = nanmean(lickNumtmp)';
    indexM= [ones(sum(indexHP==1),1);zeros(sum(indexHP==3),1)];
    [tbl,rm] = simple_mixed_anova(cat(2,meanlick(ismember(indexHP,[1 3])),meanlick(ismember(indexHP,[2 4]))), indexM, {'Drug'},{'Mouse'});
    p = table2array(tbl(4,5));
end

lickNumset = lickNumtmp; odorCueset(isnan(lickNumtmp))=NaN;
save('adjustlick_rwpn.mat','lickNumset','odorCueset');
%% decoding
mintrial = floor(min([sum(odorCueset==0) sum(odorCueset==1) sum(odorCueset==3)]));
niter = 1;
for itmp = 1:3;
decodingresult_lick = NaN(3,mintrial,niter,size(dirname,2));
for ifile = 1:size(dirname,2)
    odor1 = find(odorCueset(:,ifile)==0);
    odor2 = find(odorCueset(:,ifile)==1);
    odor3 = find(odorCueset(:,ifile)==3);
    odorset = repmat([1 2 3],mintrial-1,1); odorset = reshape(odorset, (mintrial-1)*3,1);
    
    for iter = 1:niter
        randodor1 = randperm(size(odor1,1),mintrial);
        randodor2 = randperm(size(odor2,1),mintrial);
        randodor3 = randperm(size(odor3,1),mintrial);
        randodorset = [odor1(randodor1),odor2(randodor2),odor3(randodor3)];
        for itest = 1:mintrial
            testodor = randodorset(itest,:);
            trainodor = randodorset;
            trainodor(itest,:)=[];
            trainodor = reshape(trainodor,(mintrial-1)*3,1);
            Mdl1 = fitcecoc(lickNumset(trainodor,ifile),odorset); %cue
            for icue = 1:3
                label = predict(Mdl1,lickNumset(testodor(icue),ifile));
                decodingresult_lick(icue,itest,iter,ifile) = label;
            end
        end
    end
    disp(ifile)
end
   
save(['E:\data\vHPC\Muscimol\', behaviortitle,'_adjustlickdecoding',num2str(itmp),'.mat'],'decodingresult_lick')
end

%% decoding -Rw noRw
mintrial = floor(min([sum(odorCueset==0) sum(odorCueset==1) sum(odorCueset==3)]));
niter = 1;
for itmp = 1:3;
decodingresult_lick = NaN(2,mintrial,niter,size(dirname,2));
for ifile = 1:size(dirname,2)
    odor1 = find(odorCueset(:,ifile)==find(cueKindset(:,ifile)==1)-1);
    odor2 = find(ismember(odorCueset(:,ifile),find(cueKindset(:,ifile)>2)-1));
    odorset = repmat([1 2],mintrial-1,1); odorset = reshape(odorset, (mintrial-1)*2,1);    
    for iter = 1:niter
        randodor1 = randperm(size(odor1,1),mintrial);
        randodor2 = randperm(size(odor2,1),mintrial);
        randodorset = [odor1(randodor1),odor2(randodor2)];
        for itest = 1:mintrial
            testodor = randodorset(itest,:);
            trainodor = randodorset;
            trainodor(itest,:)=[];
            trainodor = reshape(trainodor,(mintrial-1)*2,1);
            Mdl1 = fitcecoc(lickNumset(trainodor,ifile),odorset); %cue
            for icue = 1:2
                label = predict(Mdl1,lickNumset(testodor(icue),ifile));
                decodingresult_lick(icue,itest,iter,ifile) = label;
            end
        end
    end
    disp(ifile)
end
   
save(['E:\data\vHPC\Muscimol\', behaviortitle,'_adjustlickdecoding_RwnoRw',num2str(itmp),'.mat'],'decodingresult_lick')
end
%% Decoding plot
% load(['E:\data\vHPC\Muscimol\', behaviortitle,'_lickdecoding1.mat'],'decodingresult_lick')  %decodingresult_lick = decodingresult_lick(:,:,:,~indexEXCEPT);
% decoding_performance = mean([mean(decodingresult_lick(1,:,:,:)==1);mean(decodingresult_lick(2,:,:,:)==2);mean(decodingresult_lick(3,:,:,:)==3);]);
load('E:\data\vHPC\Muscimol\rwpn_lickdecoding_RwnoRw1.mat')
decoding_performance = mean([mean(decodingresult_lick(1,:,:,:)==1);mean(decodingresult_lick(2,:,:,:)==2)]);

decoding_per_dHP = [mean(decoding_performance(:,:,:,indexHP==1),3);mean(decoding_performance(:,:,:,indexHP==2),3)];
decoding_per_vHP = [mean(decoding_performance(:,:,:,indexHP==3),3);mean(decoding_performance(:,:,:,indexHP==4),3)];
decoding_per_dHP = permute(decoding_per_dHP,[4 1 2 3]); decoding_per_vHP = permute(decoding_per_vHP,[4 1 2 3]);
decoding_anova=[decoding_per_dHP;decoding_per_vHP];
mouseG = [zeros(size(decoding_per_dHP,1),1); ones(size(decoding_per_vHP,1),1)];


[tbl,rm] = simple_mixed_anova(decoding_anova, mouseG, {'Drug'},{'Mouse'});
[c]=multcompare(rm,'Drug','by','Mouse');
%%
cmap1=([0 56 66; 57 197 187; 55 0 102; 179 143 177;]./255);
if strcmp(behaviortitle,'rwpn')
    cmap1 = [55 0 102; 179 143 177]./255;
    ylimlist = [0.25 1]; cl = 1/2;
elseif strcmp(behaviortitle,'rwprob')
    cmap1 = [0 56 66; 57 197 187]./255;
    ylimlist = [0.25 0.75]; cl = 1/3;
end
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 3]);
hold on

for ii =1:2
    [~,p1] = ttest(decoding_per_dHP(:,ii),cl,'tail','right');%1/2
    [~,p2] = ttest(decoding_per_vHP(:,ii),cl,'tail','right');%1/2
    pset(:,ii) = [p1,p2];
    bar(ii-0.23,mean(decoding_per_dHP(:,ii)),'barwidth',0.4,'facecolor',cmap1(1,:),'edgecolor','none')
    scatter(ii-0.23-0.15+0.3*rand(size(decoding_per_dHP(:,ii),1),1),decoding_per_dHP(:,ii),5,[0.5 0.5 0.5],'.','MarkerFaceColor',[0.5 0.5 0.5])
    bar(ii+0.23,mean(decoding_per_vHP(:,ii)),'barwidth',0.4,'facecolor',cmap1(2,:),'edgecolor','none')
    scatter(ii+0.23-0.15+0.3*rand(size(decoding_per_vHP(:,ii),1),1),decoding_per_vHP(:,ii),5,[0.5 0.5 0.5],'.','MarkerFaceColor',[0.5 0.5 0.5])
    errorbar(ii-0.23,mean(decoding_per_dHP(:,ii)),sem(decoding_per_dHP(:,ii)),'color','k','capsize',3)
    errorbar(ii+0.23,mean(decoding_per_vHP(:,ii)),sem(decoding_per_vHP(:,ii)),'color','k','capsize',3)
end

% plot([0.77 1.77],[decoding_per_dHP(:,1) decoding_per_dHP(:,2)],'color',cmap1(1,:))
% plot([1.23 2.23],[decoding_per_vHP(:,1) decoding_per_vHP(:,2)],'color',cmap1(2,:))
xlim([0.4 2.6]); xticks([1 2]); xticklabels({})
ylim(ylimlist); yticks([0.25 0.5 0.75 1]); yticklabels([0 0 0 0])
line([0 5],[cl cl], 'color','k','linestyle',':')
% print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\decoding_',behaviortitle,'2.ai']);

%% lick mean
% lickmean = mean(lickNumset)./1.5;
% lickmean = allLickNumset./(4*420);
% lickmean = cellfun(@mean, lickboutset);
lickmean = cellfun(@length,lickboutset);
cmap1=([0 56 66; 57 197 187; 55 0 102; 179 143 177;]./255);
if strcmp(behaviortitle,'rwpn')
    cmap1 = [55 0 102; 179 143 177]./255;
elseif strcmp(behaviortitle,'rwprob')
    cmap1 = [0 56 66; 57 197 187]./255;
end
[~,p1] = ttest(lickmean(indexHP==1),lickmean(indexHP==2));
[~,p2] = ttest(lickmean(indexHP==3),lickmean(indexHP==4));
mouseG= [ones(sum(indexHP==1),1);zeros(sum(indexHP==3),1)];
[tbl,rm] = simple_mixed_anova([[lickmean(indexHP==1),lickmean(indexHP==3)]', [lickmean(indexHP==2),lickmean(indexHP==4)]'],...
    mouseG, {'Drug'},{'Mouse'});
[c]=multcompare(rm,'Mouse','by','Drug');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3]);
hold on
for ii =1:2    
    bar(ii-0.23,mean(lickmean(indexHP==ii)),'barwidth',0.4,'facecolor',cmap1(1,:),'edgecolor','none')
    scatter(ii-0.23-0.15+0.3*rand(size(lickmean(indexHP==ii),2),1),lickmean(indexHP==ii),5,[0.5 0.5 0.5],'.','MarkerFaceColor',[0.5 0.5 0.5])
    bar(ii+0.23,mean(lickmean(indexHP==ii+2)),'barwidth',0.4,'facecolor',cmap1(2,:),'edgecolor','none')
    scatter(ii+0.23-0.15+0.3*rand(size(lickmean(indexHP==ii+2),2),1),lickmean(indexHP==ii+2),5,[0.5 0.5 0.5],'.','MarkerFaceColor',[0.5 0.5 0.5])
    errorbar(ii-0.23,mean(lickmean(indexHP==ii)),sem(lickmean(indexHP==ii)),'color','k','capsize',3)
    errorbar(ii+0.23,mean(lickmean(indexHP==ii+2)),sem(lickmean(indexHP==ii+2)),'color','k','capsize',3)
end
% plot([0.77 1.77],[lickmean(indexHP==1); lickmean(indexHP==2)],'color',cmap1(1,:))
% plot([1.23 2.23],[lickmean(indexHP==3); lickmean(indexHP==4)],'color',cmap1(2,:))
xlim([0.3 2.7]); xticks([1 2]); xticklabels({})
ylim([0 400]); yticks([0 400]);yticklabels([])
% ylim([0 10]); yticks([0 10]); yticklabels([0 10])
print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\lick_bout_num_',behaviortitle,'.ai']);

%% Speed mean
lickmean = allSpeedset*(8*pi/10)./(12*420);
cmap1=([0 56 66; 57 197 187; 55 0 102; 179 143 177;]./255);
if strcmp(behaviortitle,'rwpn')
    cmap1 = [55 0 102; 179 143 177]./255;
elseif strcmp(behaviortitle,'rwprob')
    cmap1 = [0 56 66; 57 197 187]./255;
end
[~,p1] = ttest(lickmean(indexHP==1),lickmean(indexHP==2));
[~,p2] = ttest(lickmean(indexHP==3),lickmean(indexHP==4));
mouseG= [ones(sum(indexHP==1),1);zeros(sum(indexHP==3),1)];
[tbl,rm] = simple_mixed_anova([[lickmean(indexHP==1),lickmean(indexHP==3)]', [lickmean(indexHP==2),lickmean(indexHP==4)]'],...
    mouseG, {'Drug'},{'Mouse'});
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3]);
hold on
for ii =1:2    
    bar(ii-0.23,mean(lickmean(indexHP==ii)),'barwidth',0.4,'facecolor',cmap1(1,:),'edgecolor','none')
    scatter(ii-0.23-0.15+0.3*rand(size(lickmean(indexHP==ii),2),1),lickmean(indexHP==ii),5,[0.5 0.5 0.5],'.','MarkerFaceColor',[0.5 0.5 0.5])
    bar(ii+0.23,mean(lickmean(indexHP==ii+2)),'barwidth',0.4,'facecolor',cmap1(2,:),'edgecolor','none')
    scatter(ii+0.23-0.15+0.3*rand(size(lickmean(indexHP==ii+2),2),1),lickmean(indexHP==ii+2),5,[0.5 0.5 0.5],'.','MarkerFaceColor',[0.5 0.5 0.5])
    errorbar(ii-0.23,mean(lickmean(indexHP==ii)),sem(lickmean(indexHP==ii)),'color','k','capsize',3)
    errorbar(ii+0.23,mean(lickmean(indexHP==ii+2)),sem(lickmean(indexHP==ii+2)),'color','k','capsize',3)
end
% plot([0.77 1.77],[lickmean(indexHP==1); lickmean(indexHP==2)],'color',cmap1(1,:))
% plot([1.23 2.23],[lickmean(indexHP==3); lickmean(indexHP==4)],'color',cmap1(2,:))
xlim([0.3 2.7]); xticks([1 2]); xticklabels({})
ylim([0 16]); yticks([0 15]); yticklabels([0 0 0 0])
% print(f1,'-depsc','-painters',['E:\data\vHPC\Muscimol\figure\speed_',behaviortitle,'.ai']);

