%% Anticipatory lick plot
Lick_dHP = cat(1,lickz_dHP03,lickz_dHP04,lickz_dHP06,lickz_dHP07,lickz_dHP08);
Lick_vHP = cat(1,lickz_vHP06,lickz_vHP07,lickz_vHP08,lickz_vHP10,lickz_vHP11);
Lick_all= cat(1,lickz_dHP03,lickz_dHP04,lickz_dHP06,lickz_dHP07,lickz_dHP08,...
    lickz_vHP06,lickz_vHP07,lickz_vHP08,lickz_vHP10,lickz_vHP11);
day = datetime('today');
d = datestr(day,'yymmdd');
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);
hold on
bar(1,mean(Lick_all(:,1)),0.6,'Facecolor',colormap(1,:))
bar(2,mean(Lick_all(:,2)),0.6,'Facecolor',colormap(2,:))
bar(3,mean(Lick_all(:,3)),0.6,'Facecolor',colormap(3,:))
errorbar(1,mean(Lick_all(:,1)),sem(Lick_all(:,1)),'k')
errorbar(2,mean(Lick_all(:,2)),sem(Lick_all(:,2)),'k')
errorbar(3,mean(Lick_all(:,3)),sem(Lick_all(:,3)),'k')
scatter(0.75+0.5*rand(size(Lick_dHP,1),1),Lick_dHP(:,1),40,'k','o','filled')
scatter(1.75+0.5*rand(size(Lick_dHP,1),1),Lick_dHP(:,2),40,'k','o','filled')
scatter(2.75+0.5*rand(size(Lick_dHP,1),1),Lick_dHP(:,3),40,'k','o','filled')
scatter(0.75+0.5*rand(size(Lick_vHP,1),1),Lick_vHP(:,1),40,'k','x')
scatter(1.75+0.5*rand(size(Lick_vHP,1),1),Lick_vHP(:,2),40,'k','x')
scatter(2.75+0.5*rand(size(Lick_vHP,1),1),Lick_vHP(:,3),40,'k','x')
xticks([1 2 3]); xticklabels({'75%','25%','0%'})
saveas(f1,[cd,'\figure\','behavior_rwprob',d,'.tif'])


%%
triallist= {'trial_80R_ind', 'trial_80P_ind'};
zlist = {'trial_80R_ind', 'trial_80P_ind'};
titlelist = {'Before reversal','After reversal'};
cmap2 = [30 45 232; 255 30 70; 125 125 125;  232 126 58; 120 177 255; 125 125 125]/255; %blue red gray pink skyblue, gray

C_event = reshape(C_raw(oo,set_eventframe(:)),[],12*fr);
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 10]);
for isubplot = 1:1+(sum(rev_b)<nTrial)
    subplot(1,1+(sum(rev_b)<nTrial),isubplot)
    hold on
    for ii = 1:2
        eval(['trialind=',triallist{ii},';',]);
        eval(['zind=', zlist{fix((ii+1)/2)},';']);
        if isubplot==1;
            trialind=logical([trialind(logical(rev_b));zeros(sum(~rev_b),1)]);
        else
            trialind = logical([zeros(sum(rev_b),1);trialind(~logical(rev_b))]);
        end
        Ctemp = mean(C_event(zind,1:0.5*fr),2);
        CM = mean(Ctemp); CSTD=std(Ctemp);
        if rem(ii,2)==1;
            plot(mean(movmean((C_event(trialind,:)-CM)./CSTD,[4,4],2)),'linewidth',1.2,'color', cmap2(1,:))
        else
            plot(mean(movmean((C_event(trialind,:)-CM)./CSTD,[4,4],2)),'linewidth',1.2,'color', cmap2(2,:))
        end
    end
    ylim([-0.5 2])
    lim = axis;
    line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    xlim([0.5*fr 3.5*fr]); xticks(fr*[0.5]); xticklabels({'Cue onset'});
    
    hold off
    
end
%% Draw FON
% best performance dHP03:2, dHP04:2, dHP06:3, vHP06:3, vHP07:3 vHP08:3
load('fon_rwprob_20211118.mat')
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);

mouseG = {'dHP','vHP'};
for ii = 1:2
    eval(['pvalueset = pvalueset_',mouseG{ii},';']);
    day = datetime('today');
    d = datestr(day,'yymmdd');
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
    
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 9]);
    subplot(1,3,1)
    hold on
    p1=plot(FON(1,:), 'color', [34 139 34]./255, 'linewidth', 1);% o
    p4=plot(FON(4,:), 'color', [199 21 133]./255,'linewidth', 1);% v
%     p5=plot(FON(5,:), 'color', [80 80 80]./255,      'linewidth', 1);%Lick
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
    
    
    subplot(1,3,2)
    hold on
    p2=plot(FON(2,:), 'color', [34 139 34]./255, 'linewidth', 1);% o(t-1)
    p5=plot(FON(5,:), 'color', [199 21 133]./255, 'linewidth', 1);% v(t-1)
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
    
    subplot(1,3,3)
    hold on
    p3=plot(FON(3,:), 'color', [34 139 34]./255, 'linewidth', 1);% o(t-1)
    p6=plot(FON(6,:), 'color', [199 21 133]./255, 'linewidth', 1);% v(t-1)
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
    
    %     saveas(f1,['E:\data\vHPC\all\figure\cpd\rwprob\FON_C15_rwprob_',mouseG{ii},d,'.tif'])
end
%% FON shuffle
clear
shuffledir = 'E:\data\vHPC\all\shuffle\seltrial\';
dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
    'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
for imouse = 1:12
    load([shuffledir,dirn{imouse},'_rwprob_fon_shuffle_seltrial_01_1.mat'])
    eval(['srcset_r_',dirn{imouse},' = srcset_r;'])
    eval(['pvalueset_r_',dirn{imouse},' = pvalueset_r;'])    
end
save('fon_rwprob_shuffle01_sel.mat','srcset_r_dHP03','srcset_r_dHP04','srcset_r_dHP06','srcset_r_dHP07','srcset_r_dHP08','srcset_r_dHP10',...
    'srcset_r_vHP06','srcset_r_vHP07','srcset_r_vHP08','srcset_r_vHP11','srcset_r_vHP12','srcset_r_vHP14',...
    'pvalueset_r_dHP03','pvalueset_r_dHP04','pvalueset_r_dHP06','pvalueset_r_dHP07','pvalueset_r_dHP08','pvalueset_r_dHP10',...
    'pvalueset_r_vHP06','pvalueset_r_vHP07','pvalueset_r_vHP08','pvalueset_r_vHP11','pvalueset_r_vHP12','pvalueset_r_vHP14')

% clear
% shuffledir = 'E:\data\vHPC\all\shuffle\fon_shuffle_01_1\';
% dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
%     'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
% for imouse = 1:12
%     load([shuffledir,dirn{imouse},'_rwprob_fon_shuffle_outcome_1.mat'])
%     eval(['srcset_r_',dirn{imouse},' = srcset_r;'])
%     eval(['pvalueset_r_',dirn{imouse},' = pvalueset_r;'])    
% end
% save('fon_rwprob_shuffle_outcome01.mat','srcset_r_dHP03','srcset_r_dHP04','srcset_r_dHP06','srcset_r_dHP07','srcset_r_dHP08','srcset_r_dHP10',...
%     'srcset_r_vHP06','srcset_r_vHP07','srcset_r_vHP08','srcset_r_vHP11','srcset_r_vHP12','srcset_r_vHP14',...
%     'pvalueset_r_dHP03','pvalueset_r_dHP04','pvalueset_r_dHP06','pvalueset_r_dHP07','pvalueset_r_dHP08','pvalueset_r_dHP10',...
%     'pvalueset_r_vHP06','pvalueset_r_vHP07','pvalueset_r_vHP08','pvalueset_r_vHP11','pvalueset_r_vHP12','pvalueset_r_vHP14')
%% Draw value FON all case
for icase = [4]
    clearvars -except icase
    if icase==1; %value(t,t-1,t-2)
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwprob_shuffle01.mat';
        loadshufflename2= 'fon_rwprob_shuffle01.mat';
        %     figureind = 1;
        figureind = [6 7 8 9 10 11];
    elseif icase==2; %Rw
        deltafon=false; lickwithout = true;
        loadshufflename1= 'fon_rwprob_shuffle01_sel.mat';
        loadshufflename2= 'fon_rwprob_shuffle01_sel.mat';
        figureind = 1;
    elseif icase==3; %Rw (t-1,t-2)
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwprob_shuffle01.mat';
        loadshufflename2= 'fon_rwprob_shuffle01.mat';
        figureind = [2 3 4 5];
    elseif icase==4 % lick, speed
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwprob_shuffle01.mat';
        loadshufflename2= 'fon_rwprob_shuffle01.mat';
        figureind = [12 13 14];
    end
    if deltafon; deltaname='_delta'; else; deltaname=''; end
    if lickwithout; lickname='_nolick'; else; lickname='';end
    

    
    cmap = [0 56 66; 57 197 187]./255;
    load('fon_rwprob01.mat')
%     load('fon_rwprob01_sel_20220319.mat')
    pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
    pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
    pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;
    pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;
    
    % binomial test threshold
    % neuronum_dHP = size(pvalueset_dHP,3); % binomial test --> make function
    % neuron_sigbinomial_dHP = binoinv([0.05 0.975], neuronum_dHP, 0.05);
    % bin_percent_dHP = neuron_sigbinomial_dHP(2)./neuronum_dHP;
    % neuronum_vHP = size(pvalueset_vHP,3); % binomial test --> make function
    % neuron_sigbinomial_vHP = binoinv([0.05 0.975], neuronum_vHP, 0.05);
    % bin_percent_vHP = neuron_sigbinomial_vHP(2)./neuronum_vHP;
    % FON_sig_dHP = FON_dHP>bin_percent_dHP;
    % FON_sig_vHP = FON_vHP>bin_percent_vHP;
    % shuffle threshold
    cue = [4 5 6 7 8]; outcome = [1 2 3];
    load(loadshufflename1)
    pvalueset_r_dHP_tmp = cat(3,pvalueset_r_dHP03,pvalueset_r_dHP04,pvalueset_r_dHP06,pvalueset_r_dHP07,pvalueset_r_dHP08,pvalueset_r_dHP10);
    pvalueset_r_vHP_tmp = cat(3,pvalueset_r_vHP06,pvalueset_r_vHP07,pvalueset_r_vHP08,pvalueset_r_vHP11,pvalueset_r_vHP12,pvalueset_r_vHP14);
    pvalueset_r_dHP(cue,:,:) = pvalueset_r_dHP_tmp(cue+1,:,:);
    pvalueset_r_vHP(cue,:,:) = pvalueset_r_vHP_tmp(cue+1,:,:);
    load(loadshufflename2)
    pvalueset_r_dHP_tmp  = cat(3,pvalueset_r_dHP03,pvalueset_r_dHP04,pvalueset_r_dHP06,pvalueset_r_dHP07,pvalueset_r_dHP08,pvalueset_r_dHP10);
    pvalueset_r_vHP_tmp  = cat(3,pvalueset_r_vHP06,pvalueset_r_vHP07,pvalueset_r_vHP08,pvalueset_r_vHP11,pvalueset_r_vHP12,pvalueset_r_vHP14);
    pvalueset_r_dHP(outcome,:,:) = pvalueset_r_dHP_tmp(outcome+1,:,:);
    pvalueset_r_vHP(outcome,:,:) = pvalueset_r_vHP_tmp(outcome+1,:,:);
    pvaluesig_r_dHP = pvalueset_r_dHP(:,:,:)<0.05;
    pvaluesig_r_vHP = pvalueset_r_vHP(:,:,:)<0.05;
    if lickwithout
        load('E:\data\vHPC\all\rwprob\lick_neuron.mat');
        lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
        lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
        pvaluesig_dHP = pvaluesig_dHP(:,:,lick_dHP(:,1)>=0.05);
        pvaluesig_vHP = pvaluesig_vHP(:,:,lick_vHP(:,1)>=0.05);
        pvaluesig_r_dHP = pvaluesig_r_dHP(:,:,lick_dHP(:,1)>=0.05);
        pvaluesig_r_vHP = pvaluesig_r_vHP(:,:,lick_vHP(:,1)>=0.05);
    end
    if deltafon
        FON_dHP = mean(pvaluesig_dHP, 3)-mean(pvaluesig_r_dHP, 3);
        FON_vHP = mean(pvaluesig_vHP, 3)-mean(pvaluesig_r_vHP, 3);
        chis_dHP =sum(pvaluesig_dHP, 3)-sum(pvaluesig_r_dHP, 3); chis_vHP = sum(pvaluesig_vHP, 3)-sum(pvaluesig_r_vHP, 3);
    else
        FON_dHP = mean(pvaluesig_dHP, 3);
        FON_vHP = mean(pvaluesig_vHP, 3);
        bin_percent_dHP = mean(pvaluesig_r_dHP, 3);
        bin_percent_vHP =mean(pvaluesig_r_vHP, 3);
        FON_sig_dHP = FON_dHP>bin_percent_dHP;
        FON_sig_vHP = FON_vHP>bin_percent_vHP;
        chis_dHP =sum(pvaluesig_dHP, 3); chis_vHP = sum(pvaluesig_vHP, 3);
    end

    
    % figurename = {'O(t)','O(t-1)','O(t-2)','V(t)','V(t-1)','V(t-2)','Lick','Speed'};
    psig= [1 2 2 3 3, 4 4 5 5 6 6, 7 8 8]; % for figure
    xlimlist = [51 80; 1 50; 51 80; 1 50; 51 80;...
        1 50; 51 80;1 50; 51 80;1 50; 51 80;...
        35 85; 1 50; 51 80;];
    wlist = (diff(xlimlist')+1)*0.075; %0.075 0.055
    termlist = {'Outcome(t)','Outcome(t-1)1','Outcome(t-1)2','Outcome(t-2)1','Outcome(t-2)2',...
        'Value(t)1','Value(t)2','Value(t-1)1','Value(t-1)2','Value(t-2)1','Value(t-2)2',...
        'lick(t)','speed(t)1','speed(t)2'};
    
    for isig = figureind
        if ismember(psig(isig),1); ymax = 0.4;
        elseif ismember(psig(isig),[2 3]); ymax = 0.3;
        elseif ismember(psig(isig),[4]); ymax = 0.6; 
        elseif ismember(psig(isig),[5 6]); ymax = 0.4;
        else ymax=0.5; end
        ii = psig(isig);
        for tt = 1:120
            chis_set(ii,tt) = chis([chis_dHP(ii,tt),size(pvaluesig_dHP,3)-chis_dHP(ii,tt); chis_vHP(ii,tt),size(pvaluesig_vHP,3)-chis_vHP(ii,tt)]);
        end
%         f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3]);
        f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 wlist(isig) 2.5]);%
        hold on
        plot(FON_dHP(ii,:), 'color', cmap(1,:), 'linewidth', 0.5);% d
        plot(FON_vHP(ii,:), 'color', cmap(2,:), 'linewidth', 0.5);% v
        plot(0.95*ymax*ones(1,120),'color', 'k','linestyle','none','marker','*','markersize',1,'markerindices',find(chis_set(ii,:)<0.05));
        xticks([5 20 35 40 55 80]);xticklabels({});
        xlim([xlimlist(isig,:)]); ylim([-0.04 ymax]); yticks([0 ymax]); yticklabels({});
        lim = axis;
        line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
        line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
        line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
        line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
        line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
        if ~deltafon
            for tt = 1:120
                sigchis_dHP(ii,tt) = chis([sum(pvaluesig_dHP(ii,tt,:), 3) size(pvaluesig_dHP,3)-sum(pvaluesig_dHP(ii,tt,:), 3);...
                    sum(pvaluesig_r_dHP(ii,tt,:), 3) size(pvaluesig_r_dHP,3)-sum(pvaluesig_r_dHP(ii,tt,:), 3)]);
                sigchis_vHP(ii,tt) = chis([sum(pvaluesig_vHP(ii,tt,:), 3) size(pvaluesig_vHP,3)-sum(pvaluesig_vHP(ii,tt,:), 3);...
                    sum(pvaluesig_r_vHP(ii,tt,:), 3) size(pvaluesig_r_vHP,3)-sum(pvaluesig_r_vHP(ii,tt,:), 3)]);
            end
            FON_sig_dHP = sigchis_dHP<0.05; FON_sig_vHP = sigchis_vHP<0.05;
            plot(FON_dHP(ii,:), 'color', cmap(1,:), 'linestyle','none','markersize',4,'marker','.','markerindices',find(FON_sig_dHP(ii,:)));% d
            plot(FON_vHP(ii,:), 'color', cmap(2,:), 'linestyle','none','markersize',4,'marker','.','markerindices',find(FON_sig_vHP(ii,:)));% v
            plot(bin_percent_dHP(ii,:),'color',cmap(1,:), 'linestyle', ':','linewidth', 0.5)
            plot(bin_percent_vHP(ii,:),'color',cmap(2,:), 'linestyle', ':','linewidth', 0.5)
            ylim([0 ymax])
        end
%           saveas(f1,['E:\data\vHPC\all\figure\FON\rwprob\',termlist{isig},deltaname,lickname,'.tif'])
        print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwprob\',termlist{isig},deltaname,lickname,'.ai']);
    end
end

%% FON period - cue
cmap = [0 56 66; 57 197 187]./255;
ymax = 0.5;iterm=4; %value
load('period_reg_rwprob_20220320.mat')
% 7: data 8: shuffle 9: v(t) 0-25 10: v(t) 25-75 11: v(t)scramble
% 12: 0-25 all 13: 25-75 all 14: all scramble
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;
pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;
% scrumble_dHP = mean(pvaluesig_dHP(iterm,[7 9 10],:)); smethodname = 'order_v(t)';
% scrumble_vHP = mean(pvaluesig_vHP(iterm,[7 9 10],:));
% scrumble_dHP = pvaluesig_dHP(iterm,11,:);smethodname = 'scramble_v(t)';
% scrumble_vHP = pvaluesig_vHP(iterm,11,:);
scrumble_dHP = mean(pvaluesig_dHP(iterm,[12 13],:));smethodname = 'order_all';
scrumble_vHP = mean(pvaluesig_vHP(iterm,[12 13],:));
% scrumble_dHP = pvaluesig_dHP(iterm,14,:);
% scrumble_vHP = pvaluesig_vHP(iterm,14,:);smethodname = 'order_scramble';
f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1 1]);%2.5 2.5
hold on
bar(0.25,mean(pvaluesig_dHP(iterm,7,:)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
% bar(0.5,mean(scrumble_dHP),0.2,'facecolor',[160 206 222]./255,'edgecolor','none')
% bar(0.75,mean(pvaluesig_dHP(iterm,8,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')

bar(0.5,mean(pvaluesig_vHP(iterm,7,:)),0.2,'facecolor',cmap(2,:),'edgecolor','none') %1.1
% bar(1.35,mean(scrumble_vHP),0.2,'facecolor',[160 206 222]./255,'edgecolor','none')
% bar(1.6,mean(pvaluesig_vHP(iterm,8,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')

xlim([0 0.75]); xticks([]); xticklabels({});
ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
% saveas(f2,['E:\data\vHPC\all\figure\FON\rwprob\fon_period_',smethodname,'.tif'])
print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwprob\fon_period.ai']);


[p,stat]=chis([sum(pvaluesig_dHP(iterm,7,:)) size(pvaluesig_dHP(iterm,7,:),3)-sum(pvaluesig_dHP(iterm,7,:));...
    sum(pvaluesig_vHP(iterm,7,:)) size(pvaluesig_vHP(iterm,7,:),3)-sum(pvaluesig_vHP(iterm,7,:))])
[p,stat]=chis([sum(pvaluesig_dHP(iterm,7,:)) size(pvaluesig_dHP(iterm,7,:),3)-sum(pvaluesig_dHP(iterm,7,:));...
    sum(scrumble_dHP) size(pvaluesig_dHP(iterm,7,:),3)-sum(scrumble_dHP)])
[p,stat]=chis([sum(pvaluesig_dHP(iterm,7,:)) size(pvaluesig_dHP(iterm,7,:),3)-sum(pvaluesig_dHP(iterm,7,:));...
    sum(pvaluesig_dHP(iterm,8,:)) size(pvaluesig_dHP(iterm,8,:),3)-sum(pvaluesig_dHP(iterm,8,:))])
[p,stat]=chis([sum(scrumble_dHP) size(pvaluesig_dHP(iterm,7,:),3)-sum(scrumble_dHP);...
    sum(pvaluesig_dHP(iterm,8,:)) size(pvaluesig_dHP(iterm,8,:),3)-sum(pvaluesig_dHP(iterm,8,:))])
[p,stat]=chis([sum(pvaluesig_vHP(iterm,7,:)) size(pvaluesig_vHP(iterm,7,:),3)-sum(pvaluesig_vHP(iterm,7,:));...
    sum(scrumble_vHP) size(pvaluesig_vHP(iterm,7,:),3)-sum(scrumble_vHP)])
[p,stat]=chis([sum(pvaluesig_vHP(iterm,7,:)) size(pvaluesig_vHP(iterm,7,:),3)-sum(pvaluesig_vHP(iterm,7,:));...
    sum(pvaluesig_vHP(iterm,8,:)) size(pvaluesig_vHP(iterm,8,:),3)-sum(pvaluesig_vHP(iterm,8,:))])
[p,stat]=chis([sum(scrumble_vHP) size(pvaluesig_vHP(iterm,7,:),3)-sum(scrumble_vHP);...
    sum(pvaluesig_vHP(iterm,8,:)) size(pvaluesig_vHP(iterm,8,:),3)-sum(pvaluesig_vHP(iterm,8,:))])
%% FON period - value in outcome period
cmap = [0 56 66; 57 197 187]./255;
ymax = 0.7;iterm=4; %value
load('period_reg_rwprob_20220320.mat')
% 1: data 15 : scramble 16: shuffle 17 : 25-0 scramble 18; 25-75 srcamble
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;
pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;

smethodname = 'scramble_all';
scrumble_dHP = pvaluesig_dHP(iterm,15,:);
scrumble_vHP = pvaluesig_vHP(iterm,15,:);
% scrumble_dHP = mean(pvaluesig_dHP(iterm,[17 18],:));
% scrumble_dHP = mean(pvaluesig_dHP(iterm,[17 18],:));
f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1 1]);%2.5 2.5
hold on
bar(0.25,mean(pvaluesig_dHP(iterm,1,:)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
% bar(0.5,mean(scrumble_dHP),0.2,'facecolor',[160 206 222]./255,'edgecolor','none')
% bar(0.75,mean(pvaluesig_dHP(iterm,16,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')

bar(0.5,mean(pvaluesig_vHP(iterm,1,:)),0.2,'facecolor',cmap(2,:),'edgecolor','none')%1.1
% bar(1.35,mean(scrumble_vHP),0.2,'facecolor',[160 206 222]./255,'edgecolor','none')
% bar(1.6,mean(pvaluesig_vHP(iterm,16,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')

xlim([0 0.75]); xticks([]); xticklabels({});
ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
% saveas(f2,['E:\data\vHPC\all\figure\FON\rwprob\fon_period_',smethodname,'.tif'])
print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwprob\fon_period_value_outperiod.ai']);

[sum(pvaluesig_vHP(iterm,1,:)) size(pvaluesig_vHP(iterm,1,:),3)-sum(pvaluesig_vHP(iterm,1,:));...
    sum(scrumble_vHP) size(pvaluesig_vHP(iterm,1,:),3)-sum(scrumble_vHP)]

[p,stat]=chis([sum(pvaluesig_dHP(iterm,1,:)) size(pvaluesig_dHP(iterm,1,:),3)-sum(pvaluesig_dHP(iterm,1,:));...
    sum(pvaluesig_vHP(iterm,1,:)) size(pvaluesig_vHP(iterm,1,:),3)-sum(pvaluesig_vHP(iterm,1,:))])
[p,stat]=chis([sum(pvaluesig_dHP(iterm,1,:)) size(pvaluesig_dHP(iterm,1,:),3)-sum(pvaluesig_dHP(iterm,1,:));...
    sum(scrumble_dHP) size(pvaluesig_dHP(iterm,1,:),3)-sum(scrumble_dHP)])
[p,stat]=chis([sum(pvaluesig_dHP(iterm,1,:)) size(pvaluesig_dHP(iterm,1,:),3)-sum(pvaluesig_dHP(iterm,1,:));...
    sum(pvaluesig_dHP(iterm,16,:)) size(pvaluesig_dHP(iterm,16,:),3)-sum(pvaluesig_dHP(iterm,16,:))])
[p,stat]=chis([sum(scrumble_dHP) size(pvaluesig_dHP(iterm,1,:),3)-sum(scrumble_dHP);...
    sum(pvaluesig_dHP(iterm,16,:)) size(pvaluesig_dHP(iterm,16,:),3)-sum(pvaluesig_dHP(iterm,16,:))])
[p,stat]=chis([sum(pvaluesig_vHP(iterm,1,:)) size(pvaluesig_vHP(iterm,1,:),3)-sum(pvaluesig_vHP(iterm,1,:));...
    sum(scrumble_vHP) size(pvaluesig_vHP(iterm,1,:),3)-sum(scrumble_vHP)])
[p,stat]=chis([sum(pvaluesig_vHP(iterm,1,:)) size(pvaluesig_vHP(iterm,1,:),3)-sum(pvaluesig_vHP(iterm,1,:));...
    sum(pvaluesig_vHP(iterm,16,:)) size(pvaluesig_vHP(iterm,16,:),3)-sum(pvaluesig_vHP(iterm,16,:))])
[p,stat]=chis([sum(scrumble_vHP) size(pvaluesig_vHP(iterm,1,:),3)-sum(scrumble_vHP);...
    sum(pvaluesig_vHP(iterm,16,:)) size(pvaluesig_vHP(iterm,16,:),3)-sum(pvaluesig_vHP(iterm,16,:))])
%% FON-perion - Outcome
cmap = [0 56 66; 57 197 187]./255;
load('period_reg_rwprob_20220320.mat')
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;% remove lick p
pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;% remove lick p
load('E:\data\vHPC\all\rwprob\lick_neuron.mat');
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
pvaluesig_dHP = pvaluesig_dHP(:,:,lick_dHP(:,1)>=0.05);
pvaluesig_vHP = pvaluesig_vHP(:,:,lick_vHP(:,1)>=0.05);
iterm = 1; ymax = 0.3;

for iperiod = 1:3
    f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1 2.5]);%
    
    hold on
    bar(0.25,mean(pvaluesig_dHP(iterm,2*iperiod-1,:)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
    bar(0.5,mean(pvaluesig_vHP(iterm,2*iperiod-1,:)),0.2,'facecolor',cmap(2,:),'edgecolor','none')
    
    xlim([0.0 0.75]); xticks([]); xticklabels({});
    ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
    %         saveas(f2,['E:\data\vHPC\all\figure\FON\rwprob\fig4_fon_periodRw',num2str(iperiod),'.tif'])
    print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwprob\fon_periodRw',num2str(iperiod),'.ai']);
    dtmp = sum(pvaluesig_dHP(iterm,2*iperiod-1,:));
    vtmp = sum(pvaluesig_vHP(iterm,2*iperiod-1,:));
    [p,stat] = chis([dtmp length(pvaluesig_dHP(iterm,2*iperiod-1,:))-dtmp;...
         vtmp length(pvaluesig_vHP(iterm,2*iperiod-1,:))-vtmp])
end

%% SRC period
cmap = [0 56 66; 57 197 187]./255;

load('period_reg_rwprob_20220320.mat')
srcset_dHP = cat(3,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
srcset_vHP = cat(3,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);

abs_dHP = permute(abs(srcset_dHP),[3 2 1]);
abs_vHP = permute(abs(srcset_vHP),[3 2 1]);
load('E:\data\vHPC\all\rwprob\lick_neuron.mat');
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
abs_dHP = abs_dHP(lick_dHP(:,1)>=0.05,:,:);
abs_vHP = abs_vHP(lick_vHP(:,1)>=0.05,:,:);
iterm = 1; ymax = 0.4;

for iperiod = 1:3
    f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 1]);%
    
    hold on
    bar(0.25,mean(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
    errorbar(0.25,mean(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm)),sem(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm)),'color','k','Capsize',3)
%     bar(0.5,mean(abs_dHP(:,2*iperiod,iterm)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
%     errorbar(0.5,mean(abs_dHP(:,2*iperiod,iterm)),sem(abs_dHP(:,2*iperiod,iterm)),'color','k','Capsize',3) % shuffle
    
    bar(0.85,mean(abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm)),0.2,'facecolor',cmap(2,:),'edgecolor','none')
    errorbar(0.85,mean(abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm)),sem(abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm)),'color','k','Capsize',3)
%     bar(1.1,mean(abs_vHP(:,2*iperiod,iterm)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
%     errorbar(1.1,mean(abs_vHP(:,2*iperiod,iterm)),sem(abs_vHP(:,2*iperiod,iterm)),'color','k','Capsize',3) % shuffle
    
    xlim([0 1.35]); xticks([0.6 0.9]); xticklabels({});
    ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
%             saveas(f2,['E:\data\vHPC\all\figure\FON\rwprob\fig4_src_periodRw',num2str(iperiod),'.tif'])
    %         print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\src_period',tname,num2str(tmpid),'.ai']);
    
    
%     anova_d = [abs_dHP(:,2*iperiod-1,iterm);abs_dHP(:,2*iperiod,iterm)];
%     anova_v = [abs_vHP(:,2*iperiod-1,iterm);abs_vHP(:,2*iperiod,iterm)];
%     mouseG  = [zeros(size(anova_d));ones(size(anova_v))];
%     termG = [1*ones(size(abs_dHP(:,2*iperiod-1,iterm))); 2*ones(size(abs_dHP(:,2*iperiod,iterm))); ...
%         1*ones(size(abs_vHP(:,2*iperiod-1,iterm))); 2*ones(size(abs_vHP(:,2*iperiod,iterm))); ];
%     [~,tbl,stat] = anovan([anova_d;anova_v],{mouseG,termG},'model','interaction','varnames',{'mouse','term'},'display','off');
%     c = multcompare(stat, 'dimension', [1 2],'display','off');
%     cset{1,iperiod} = c;
[~,p]=ttest2(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm),abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm))
end


%% Draw value src all case
cmap = [0 56 66; 57 197 187]./255;
load('fon_rwprob01.mat')
srcset_dHP = cat(3,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
srcset_vHP = cat(3,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);

abs_dHP = permute(abs(srcset_dHP),[3 2 1]);
abs_vHP = permute(abs(srcset_vHP),[3 2 1]);

load('fon_rwprob_shuffle_outcome01.mat')
srcset_r_dHP = cat(3,srcset_r_dHP03,srcset_r_dHP04,srcset_r_dHP06,srcset_r_dHP07,srcset_r_dHP08,srcset_r_dHP10);
srcset_r_vHP = cat(3,srcset_r_vHP06,srcset_r_vHP07,srcset_r_vHP08,srcset_r_vHP11,srcset_r_vHP12,srcset_r_vHP14);
abs_r_dHP = permute(abs(srcset_r_dHP),[3 2 1]);
abs_r_vHP = permute(abs(srcset_r_vHP),[3 2 1]);
% abs_dHP = abs_dHP-abs_r_dHP;
% abs_vHP = abs_vHP-abs_r_vHP;

load('E:\data\vHPC\all\rwprob\lick_neuron.mat');
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
abs_dHP = abs_dHP(lick_dHP(:,1)>=0.05,:,:);
abs_vHP = abs_vHP(lick_vHP(:,1)>=0.05,:,:);

figurename = {'O(t)','O(t-1)','O(t-2)','V(t)','V(t-1)','V(t-2)','Lick','Speed'};

for ii =1%1:8
    for tt = 1:120
        sig_t(1,tt) = ttest2(abs_dHP(:,tt,ii),abs_vHP(:,tt,ii));
    end
    ymax = 0.2;
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
    hold on
    stdshade(abs_dHP(:,:,ii),0.3,cmap(1,:))
    stdshade(abs_vHP(:,:,ii),0.3,cmap(2,:))
    stdshade(abs_r_dHP(:,:,ii),0.3,cmap(1,:),[],':')
    stdshade(abs_r_vHP(:,:,ii),0.3,cmap(2,:),[],':')
    ylim([0 ymax])
    plot(0.9*ymax*ones(1,120),'color', 'k','linestyle','none','marker','*','markersize',1,'markerindices',find(sig_t));
    lim = axis;
    xlim([50 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    % line([lim(1) lim(2)], [bin_percent_dHP bin_percent_dHP], 'color',[148 69 40]./255, 'linestyle', ':','linewidth', 0.75)
    % line([lim(1) lim(2)], [bin_percent_vHP bin_percent_vHP], 'color',[25 148 123]./255, 'linestyle', ':','linewidth', 0.75)
    saveas(f1,['E:\data\vHPC\all\figure\FON\rwprob\fig4_src-raw',figurename{ii},'.tif'])
    % print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\',figurename{ii},'.ai']);
end
%% fon value in delay1
load('period_reg_rwprob.mat')
pvalueset_dHP = cat(2,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(2,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
p_d = sum(pvalueset_dHP(5,:,1)<0.05,2); n_d = size(pvalueset_dHP,2);
p_v = sum(pvalueset_vHP(5,:,1)<0.05,2); n_v = size(pvalueset_vHP,2);

chis([p_d(1) n_d-p_d(1); p_v(1) n_v-p_v(1);])
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.5 2.5]);
hold on
bar(1,p_d(1)./n_d,'facecolor',cmap(1,:),'barwidth',0.8);
bar(2,p_v(1)./n_v,'facecolor',cmap(2,:),'barwidth',0.8);

xticks([]); xlim([0.4 2.6]); yticks([0 .5 1]); ylim([0 1])
print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwprob\barg.ai']);
%% Draw CPD
% best performance dHP03:2, dHP04:2, dHP06:3, vHP06:3, vHP07:3 vHP08:3
% load('cpd_rwprob'); load('cpd_rwprob_withv')
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
% load('cpd_rwprob_withv')
% sse_set_dHP_withv = cat(3,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
% sse_set_vHP_withv = cat(3,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP13);
mouseG = {'dHP','vHP'};
for ii = 1:2
    eval(['sse_set = sse_set_',mouseG{ii},';']);
    day = datetime('today');
    d = datestr(day,'yymmdd');
    cpd = (sse_set(2:end,:,:)-sse_set(1,:,:))./sse_set(2:end,:,:);
    cpd = permute(cpd,[3 2 1]);
    yl = [0 0.04];
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 17 6]);
    hxa(1) = subplot(1,3,1);
    hold on
    stdshade(cpd(:,:,7), 0.3, [80 80 80]./255)%Lick
%     stdshade(cpd(:,:,6), 0.3, [160 82 45]./255)%speed
    stdshade(cpd(:,:,1), 0.3, [34 139 34]./255)% o
    stdshade(cpd(:,:,4), 0.3, [199 21 133]./255)% v
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [yl(1) yl(2)],'color','k', 'linestyle', '-.'); %2.5 outcome
    title('t','FontSize',14)
    hold off
    
    hxa(2) = subplot(1,3,2);
    hold on
    stdshade(cpd(:,:,2), 0.3, [34 139 34]./255)%o t-1
    stdshade(cpd(:,:,5), 0.3, [199 21 133]./255)%v t-1
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [yl(1) yl(2)],'color','k', 'linestyle', '-.'); %2.5 outcome
    title('t-1','FontSize',14)
    hold off
  
    
    hxa(3) = subplot(1,3,3);
    hold on
    stdshade(cpd(:,:,3), 0.3, [34 139 34]./255)%o t-1
    stdshade(cpd(:,:,6), 0.3, [199 21 133]./255)%v t-1
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [yl(1) yl(2)],'color','k', 'linestyle', '-.'); %2.5 outcome
    title('t-2','FontSize',14)
    hold off

    xticks(hxa,[12 27 47 67 100]); xticklabels(hxa, {'C', 'D1','D2','O','ITI'})
    ylim(hxa, yl); xlim(hxa, [1 120]);
    saveas(f1,['E:\data\vHPC\all\figure\cpd\rwprob\','CPD_rwprob_',mouseG{ii},d,'.tif'])
end
%% CPD anova
%Rwcue, Pncue, Rw, Pn, Lick, (t-1) (t-2)
% load('cpd_rwprob'); load('cpd_rwprob_withv')
% sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07(1:6,:,:),sse_set_dHP08(1:6,:,:),sse_set_dHP10(1:6,:,:));
% sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07(1:6,:,:),sse_set_vHP08(1:6,:,:),sse_set_vHP11(1:6,:,:),sse_set_vHP12(1:6,:,:),sse_set_vHP13(1:6,:,:));
% load('cpd_rwprob_withv')
% sse_set_dHP_withv = cat(3,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
% sse_set_vHP_withv = cat(3,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP13);
% load('cpd_rwprob_20211118.mat')
load('cpd_with rw-101.mat');lickwithout = false; subt='_'; shufName = 'cpd_rwprob_shuffle_ttest_Outcome.mat';
% load('cpd_with rw-101.mat');lickwithout = false; subt='_'; shufName = 'cpd_rwprob_shuffle_ttest.mat';
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
load('lick_neuron.mat'); lickwithout = true;subt='_nolick';
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);
mouseG = {'dHP','vHP'};
for ii = 1:2
    eval(['sse_set = sse_set_',mouseG{ii},';']);
    day = datetime('today');
    d = datestr(day,'yymmdd');
    cpd = (sse_set(2:end,:,:)-sse_set(1,:,:))./sse_set(2:end,:,:);
    eval(['cpd_',mouseG{ii},'= permute(cpd,[3 2 1]);']);    
    a=[];
    for iterm = 1:size(cpd,1)
    a(:,:,iterm) = iterm*ones(size(cpd,3),size(cpd,2));
    end
    eval(['cpd_anova_term_',mouseG{ii},'=a;']);
    clear a
end
cpd_anova = cat(1,cpd_dHP,cpd_vHP);
cpd_anova_term = cat(1,cpd_anova_term_dHP,cpd_anova_term_vHP);
cpd_anova_group = cat(1,zeros(size(cpd_dHP)),ones(size(cpd_vHP)));
a =size(cpd_dHP,3)+1:2*size(cpd_dHP,3)-1; a =flip(a); a = cumsum(a);
pickt = [size(cpd_dHP,3),a+size(cpd_dHP,3)];
panovaset = NaN(size(cpd_dHP,2),size(cpd_dHP,3));
for itime = 1:size(cpd_dHP,2);
    
    cpd_temp = reshape(cpd_anova(:,itime,:),[],1);
    cpd_temp_term = reshape(cpd_anova_term(:,itime,:),[],1);
    cpd_temp_group = reshape(cpd_anova_group(:,itime,:),[],1);
    [p,tlb,stats] = anovan(cpd_temp,{cpd_temp_term,cpd_temp_group},'model','interaction','varnames',{'Term','Group'},...
        'display','off');
    if p(3)<0.05;
        c = multcompare(stats,'Dimension',[1 2],'display','off');
        panovaset(itime,:) = c(pickt,6);
    end
        
end
%% cpd delta
day = datetime('today');
d = datestr(day,'yymmdd');
%Rwcue, Pncue, Rw, Pn, Lick, (t-1) (t-2)
% load('cpd_with rw-101.mat');lickwithout = false; subt='_'; shufName = 'cpd_rwprob_shuffle_ttest_Outcome.mat';
% load('cpd_with rw-101.mat');lickwithout = false; subt='_'; shufName = 'cpd_rwprob_shuffle_ttest.mat';
load('cpd_rwprob01.mat');lickwithout = false; subt='_'; shufName = 'cpd_rwprob_shuffle_ttest_Outcome_01.mat';
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
load('lick_neuron.mat'); lickwithout = true;subt='_nolick';
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);
cpd_dHP = (sse_set_dHP(2:end,:,:)-sse_set_dHP(1,:,:))./sse_set_dHP(2:end,:,:); cpd_dHP = permute(cpd_dHP, [3 2 1]);
cpd_vHP = (sse_set_vHP(2:end,:,:)-sse_set_vHP(1,:,:))./sse_set_vHP(2:end,:,:); cpd_vHP = permute(cpd_vHP, [3 2 1]);
mouseG = {'dHP','vHP'};
load(shufName)
cpd_dHP_r = cpd_dHP_r(lick_dHP(:,1)>=0.05,:,:);
cpd_vHP_r = cpd_vHP_r(lick_vHP(:,1)>=0.05,:,:);

cpd_dHP = cpd_dHP-cpd_dHP_r;
cpd_vHP = cpd_vHP-cpd_vHP_r;
for tt = 1:120
    for iterm = 1:size(cpd_dHP,3)
        [~,p]= ttest2(cpd_dHP(:,tt,iterm),cpd_vHP(:,tt,iterm));
        panovaset(tt,iterm) =p;
    end
end

%% cpd shuffle
shuffledir = 'E:\data\vHPC\all\shuffle\outcome_01\';
dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
    'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
for imouse = 1:12
    load([shuffledir,dirn{imouse},'_rwprob_cpd_shuffle_outcome_1.mat'])
    eval(['sse_set_r_',dirn{imouse},' = sse_set_r;'])
    cpd_r =(sse_set_r(2:end,:,:,:)-sse_set_r(1,:,:,:))./sse_set_r(2:end,:,:,:);
    eval(['cpd_r_',dirn{imouse},'=cpd_r;'])
    
end
save('cpd_rwprob_shuffle_Outcome_01.mat','cpd_r_dHP03','cpd_r_dHP04','cpd_r_dHP06','cpd_r_dHP07','cpd_r_dHP08','cpd_r_dHP10',...
    'cpd_r_vHP06','cpd_r_vHP07','cpd_r_vHP08','cpd_r_vHP11','cpd_r_vHP12','cpd_r_vHP14')

%permutation

% load('cpd_rwprob_shuffle.mat','cpd_r_dHP03','cpd_r_dHP04','cpd_r_dHP06','cpd_r_dHP07','cpd_r_dHP08','cpd_r_dHP10',...
%     'cpd_r_vHP06','cpd_r_vHP07','cpd_r_vHP08','cpd_r_vHP11','cpd_r_vHP12','cpd_r_vHP14');
% cpd_dHP_r = cat(3,cpd_r_dHP03,cpd_r_dHP04,cpd_r_dHP06,cpd_r_dHP07,cpd_r_dHP08,cpd_r_dHP10);
% cpd_vHP_r = cat(3,cpd_r_vHP06,cpd_r_vHP07,cpd_r_vHP08,cpd_r_vHP11,cpd_r_vHP12,cpd_r_vHP14);
% meancpd_d = mean(cpd_dHP,1); meancpdr_d = mean(cpd_dHP_r,3); % neuron mean
% meancpd_v = mean(cpd_vHP,1); meancpdr_v = mean(cpd_vHP_r,3); % neuron mean
% chancelevel = size(sse_set_r,4)*0.975;
% sig_d = sum(permute(meancpd_d,[3 2 1])>meancpdr_d,4)>chancelevel;
% sig_v = sum(permute(meancpd_v,[3 2 1])>meancpdr_v,4)>chancelevel;
% sig_d = sig_d'; sig_v = sig_v';
% save('cpd_rwprob_shuffle_permutation.mat','sig_d','sig_v')

%% ttset
cpd_dHP_r = cat(3,cpd_r_dHP03,cpd_r_dHP04,cpd_r_dHP06,cpd_r_dHP07,cpd_r_dHP08,cpd_r_dHP10);
cpd_vHP_r = cat(3,cpd_r_vHP06,cpd_r_vHP07,cpd_r_vHP08,cpd_r_vHP11,cpd_r_vHP12,cpd_r_vHP14);
cpd_dHP_r = permute(cpd_dHP_r(:,:,:,1),[3 2 1]);
cpd_vHP_r = permute(cpd_vHP_r(:,:,:,1),[3,2,1]);
save('cpd_rwprob_shuffle_ttest_Outcome_01.mat','cpd_dHP_r','cpd_vHP_r')
%% cpd anova plot
load(shufName)
if lickwithout ==true; 
    cpd_dHP_r = cpd_dHP_r(lick_dHP(:,1)>=0.05,:,:); 
    cpd_vHP_r = cpd_vHP_r(lick_vHP(:,1)>=0.05,:,:);
end
for iterm = 1:size(cpd_dHP_r,3)
    for tt = 1:120
        [~,p]= ttest(cpd_dHP(:,tt,iterm),cpd_dHP_r(:,tt,iterm));
        pdHP(tt,iterm) = p;
        [~,p]= ttest(cpd_vHP(:,tt,iterm),cpd_vHP_r(:,tt,iterm));
        pvHP(tt,iterm) = p;
    end
end
%%
cmap = [0 56 66; 57 197 187]./255;
% psig= sum(panovaset<0.05);
% termlist = {'Outcome(t)','Outcome(t-1)','Outcome(t-2)','Value(t)','Value(t-1)','Value(t-2)','Lick','Speed'};
psig= [1 2 2 3 3, 4 4 5 5 6 6]; % for figure
xlimlist = [51 80; 1 50; 51 80; 1 50; 51 80;...
    1 50; 51 80;1 50; 51 80;1 50; 51 80;];
wlist = (diff(xlimlist')+1)*0.06;
termlist = {'Outcome(t)','Outcome(t-1)1','Outcome(t-1)2','Outcome(t-2)1','Outcome(t-2)2',...
    'Value(t)1','Value(t)2','Value(t-1)1','Value(t-1)2','Value(t-2)1','Value(t-2)2'};
for isig =1%1:length(psig)
    if ismember(psig(isig),[1 2 3]); ymax = 0.012;
    elseif ismember(psig(isig),[4 5 6]); ymax = 0.05; end
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
%     f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 wlist(isig) 3]);%
    hold on
    stdshade(cpd_dHP(:,:,psig(isig)), 0.3, cmap(1,:))%dorsal %purple
    stdshade(cpd_vHP(:,:,psig(isig)), 0.3, cmap(2,:))%ventral %green
%     stdshade(cpd_dHP_r(:,:,psig(isig)), 0.3, cmap(1,:),[],':') %dorsal shuffle
%     stdshade(cpd_vHP_r(:,:,psig(isig)), 0.3, cmap(2,:),[],':')%ventral shuffle
    %     title(termlist{isig});
    plot((panovaset(:,psig(isig))<0.05)-1+0.99*ymax,'color',[235 110 30]./255,'Marker','*','MarkerSize',1,'linestyle','none','linewidth',0.25)
%     plot((pdHP(:,psig(isig))<0.05)-1+0.95*ymax,'color',cmap(1,:),'Marker','*','MarkerSize',1,'linestyle','none','linewidth',0.25)
%     plot((pvHP(:,psig(isig))<0.05)-1+0.91*ymax,'color',cmap(2,:),'Marker','*','MarkerSize',1,'linestyle','none','linewidth',0.25)
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
%     line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    xticks([5 20 35 40 55 80]);xticklabels({});  %xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    xlim([xlimlist(isig,:)]); ylim([0 ymax]); yticks([0 0.01]); yticklabels({});  %xlim([1 120]); ylim([0 0.04])
    hold off
%     print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwprob\fig4',termlist{isig},subt,num2str(d),'.ai']);
    saveas(f1,['E:\data\vHPC\all\figure\cpd\rwprob\fig4_delta',termlist{isig},d,'.tif'])
%     close all

end
%% cpd anova plot - period
cmap = [0 56 66; 57 197 187]./255;
%cue
load('cpd_rwprob_period_withrw-101.mat')
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
cpd_dHP = (sse_set_dHP(2:end,:,:)-sse_set_dHP(1,:,:))./sse_set_dHP(2:end,:,:);
cpd_vHP = (sse_set_vHP(2:end,:,:)-sse_set_vHP(1,:,:))./sse_set_vHP(2:end,:,:);

tname = 'value'; iterm = 4; ymax = 0.05;
for ii = 2
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.8 2.5]);%
    hold on
    d_scramble = mean(cpd_dHP(iterm,[ii,ii+4,ii+6],:),2);
    bar(0.25,mean(cpd_dHP(iterm,ii,:),3),0.2,'facecolor',cmap(1,:),'edgecolor','none')
    errorbar(0.25,mean(cpd_dHP(iterm,ii,:),3),sem(cpd_dHP(iterm,ii,:)),'color','k','Capsize',3)
    bar(0.5,mean(d_scramble,3),0.2,'facecolor',[160 206 222]./255,'edgecolor','none')
    errorbar(0.5,mean(d_scramble,3),sem(d_scramble),'color','k','Capsize',3) % scramble
    bar(0.75,mean(cpd_dHP(iterm,ii+2,:),3),0.2,'facecolor',[180 180 180]./255,'edgecolor','none')
    errorbar(0.75,mean(cpd_dHP(iterm,ii+2,:),3),sem(cpd_dHP(iterm,ii+2,:)),'color','k','Capsize',3) % all shuffle
    
    v_scramble = mean(cpd_vHP(iterm,[ii+2,ii+4,ii+6],:),2);
    bar(1.1,mean(cpd_vHP(iterm,ii,:),3),0.2,'facecolor',cmap(2,:),'edgecolor','none')
    errorbar(1.1,mean(cpd_vHP(iterm,ii,:),3),sem(cpd_vHP(iterm,ii,:)),'color','k','Capsize',3)
    bar(1.35,mean(v_scramble,3),0.2,'facecolor',[160 206 222]./255,'edgecolor','none')
    errorbar(1.35,mean(v_scramble,3),sem(v_scramble),'color','k','Capsize',3) % scramble
    bar(1.6,mean(cpd_vHP(iterm,ii+2,:),3),0.2,'facecolor',[180 180 180]./255,'edgecolor','none')
    errorbar(1.6,mean(cpd_vHP(iterm,ii+2,:),3),sem(cpd_vHP(iterm,ii+2,:)),'color','k','Capsize',3) % all shuffle
    xlim([0.1 1.75]); xticks([]); xticklabels({});
    ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
    % saveas(f1,['E:\data\vHPC\all\figure\cpd\rwprob\mean_',tname,'.tif'])
%     print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwprob\mean_',tname,num2str(ii),'.ai']);
end
% 
anova_d = [squeeze(cpd_dHP(iterm,ii,:));squeeze(d_scramble);squeeze(cpd_dHP(iterm,ii+2,:))];
anova_v = [squeeze(cpd_vHP(iterm,ii,:));squeeze(v_scramble);squeeze(cpd_vHP(iterm,ii+2,:))];
mouseG  = [zeros(size(anova_d));ones(size(anova_v))];
termG = [1*ones(size(squeeze(cpd_dHP(iterm,ii,:)))); 2*ones(size(squeeze(cpd_dHP(iterm,ii,:)))); 3*ones(size(squeeze(cpd_dHP(iterm,ii,:))));...
    1*ones(size(squeeze(cpd_vHP(iterm,ii,:)))); 2*ones(size(squeeze(cpd_vHP(iterm,ii,:)))); 3*ones(size(squeeze(cpd_vHP(iterm,ii,:))))];
[~,tbl,stat] = anovan([anova_d;anova_v],{mouseG,termG},'model','interaction','varnames',{'mouse','term'});
c = multcompare(stat, 'dimension', [1 2]);
%% outcome
cmap = [0 56 66; 57 197 187]./255;
load('cpd_rwprob_period_withrw-101.mat')
load('cpd_rwprob_period_withrw01.mat')
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
load('lick_neuron.mat'); 
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);
cpd_dHP = (sse_set_dHP(2:end,:,:)-sse_set_dHP(1,:,:))./sse_set_dHP(2:end,:,:);
cpd_vHP = (sse_set_vHP(2:end,:,:)-sse_set_vHP(1,:,:))./sse_set_vHP(2:end,:,:);

tname = {'outcome_all','outcome_1s','outcome_2s'};
for ii = 1:3
f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 1]);%
 iterm = 1; ymax = 0.02; 

hold on
bar(0.25,mean(cpd_dHP(iterm,9+2*(ii-1),:)-cpd_dHP(iterm,10+2*(ii-1),:),3),0.2,'facecolor',cmap(1,:),'edgecolor','none')
errorbar(0.25,mean(cpd_dHP(iterm,9+2*(ii-1),:)-cpd_dHP(iterm,10+2*(ii-1),:),3),sem(cpd_dHP(iterm,9+2*(ii-1),:)-cpd_dHP(iterm,10+2*(ii-1),:)),'color','k','Capsize',3)
% bar(0.5,mean(cpd_dHP(iterm,10+2*(ii-1),:),3),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
% errorbar(0.5,mean(cpd_dHP(iterm,10+2*(ii-1),:),3),sem(cpd_dHP(iterm,10+2*(ii-1),:)),'color','k','Capsize',3) % all shuffle

bar(0.85,mean(cpd_vHP(iterm,9+2*(ii-1),:)-cpd_vHP(iterm,10+2*(ii-1),:),3),0.2,'facecolor',cmap(2,:),'edgecolor','none')
errorbar(0.85,mean(cpd_vHP(iterm,9+2*(ii-1),:)-cpd_vHP(iterm,10+2*(ii-1),:),3),sem(cpd_vHP(iterm,9+2*(ii-1),:)-cpd_vHP(iterm,10+2*(ii-1),:)),'color','k','Capsize',3)
% bar(1.1,mean(cpd_vHP(iterm,10+2*(ii-1),:),3),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
% errorbar(1.1,mean(cpd_vHP(iterm,10+2*(ii-1),:),3),sem(cpd_vHP(iterm,10+2*(ii-1),:)),'color','k','Capsize',3) % all shuffle

% plot(0.9,mean(cpd_vHP(iterm,9,:),3),'color',cmap(2,:),'marker','.','markersize',5,'markerfacecolor',[8 128 75]./255)

xlim([0 1.35]); xticks([0.6 0.9]); xticklabels({});
ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
saveas(f2,['E:\data\vHPC\all\figure\cpd\rwprob\fig4_mean_',tname{ii},'.tif'])
% print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwprob\mean_',tname{ii},'.ai']);

% anova_d = [squeeze(cpd_dHP(iterm,10+2*(ii-1),:));squeeze(cpd_dHP(iterm,9+2*(ii-1),:))];
% anova_v = [squeeze(cpd_vHP(iterm,10+2*(ii-1),:));squeeze(cpd_vHP(iterm,9+2*(ii-1),:))];
% mouseG  = [zeros(size(anova_d));ones(size(anova_v))];
% termG = [1*ones(size(squeeze(cpd_dHP(iterm,9,:)))); 2*ones(size(squeeze(cpd_dHP(iterm,10,:))));...
%     1*ones(size(squeeze(cpd_vHP(iterm,9,:)))); 2*ones(size(squeeze(cpd_vHP(iterm,10,:))))];
% [~,tbl,stat] = anovan([anova_d;anova_v],{mouseG,termG},'model','interaction','varnames',{'mouse','term'},'display','off');
% tblset{ii} = tbl;
% c = multcompare(stat, 'dimension', [1 2],'display','off');
% cset{ii} = c;
[~,p] = ttest2(cpd_dHP(iterm,9+2*(ii-1),:)-cpd_dHP(iterm,10+2*(ii-1),:),cpd_vHP(iterm,9+2*(ii-1),:)-cpd_vHP(iterm,10+2*(ii-1),:))
end

%% cpd anova plot - speed
%for speed effect - dHP
load('cpd_with rw-101.mat');lickwithout = false;
sse_set_dHP = cat(3,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_dHP  = sse_set_dHP (1:8,:,:);
load('CPD_no_speed');
sse_set_dHP_nospeed = cat(3,sse_set__nospeed_dHP07,sse_set__nospeed_dHP08,sse_set__nospeed_dHP10);
%for speed effect - vHP
load('cpd_with rw-101.mat');lickwithout = false;
sse_set_vHP = cat(3,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
sse_set_vHP  = sse_set_vHP(1:8,:,:);
load('CPD_no_speed');
sse_set_vHP_nospeed = cat(3,sse_set__nospeed_vHP07,sse_set__nospeed_vHP08,sse_set__nospeed_vHP11,sse_set__nospeed_vHP12,sse_set__nospeed_vHP14);
load('lick_neuron.mat'); lickwithout = true;subt='_nolick';
lick_dHP = cat(1, pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
sse_set_dHP_nospeed = sse_set_dHP_nospeed(:,:,lick_dHP(:,1)>=0.05);
sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);
sse_set_vHP_nospeed = sse_set_vHP_nospeed(:,:,lick_vHP(:,1)>=0.05);

mouseG = {'dHP','vHP'};
for ii = 1:2
    eval(['sse_set = sse_set_',mouseG{ii},'_nospeed;']);
    day = datetime('today');
    d = datestr(day,'yymmdd');
    cpd = (sse_set(2:end,:,:)-sse_set(1,:,:))./sse_set(2:end,:,:);
    eval(['cpd_',mouseG{ii},'= permute(cpd,[3 2 1]);']);    
    a=[];
    for iterm = 1:size(cpd,1)
    a(:,:,iterm) = iterm*ones(size(cpd,3),size(cpd,2));
    end
    eval(['cpd_anova_term_',mouseG{ii},'=a;']);
    clear a
end
cpd_anova = cat(1,cpd_dHP,cpd_vHP);
cpd_anova_term = cat(1,cpd_anova_term_dHP,cpd_anova_term_vHP);
cpd_anova_group = cat(1,zeros(size(cpd_dHP)),ones(size(cpd_vHP)));
a =size(cpd_dHP,3)+1:2*size(cpd_dHP,3)-1; a =flip(a); a = cumsum(a);
pickt = [size(cpd_dHP,3),a+size(cpd_dHP,3)];
panovaset = NaN(size(cpd_dHP,2),size(cpd_dHP,3));
for itime = 1:size(cpd_dHP,2);
    
    cpd_temp = reshape(cpd_anova(:,itime,:),[],1);
    cpd_temp_term = reshape(cpd_anova_term(:,itime,:),[],1);
    cpd_temp_group = reshape(cpd_anova_group(:,itime,:),[],1);
    [p,tlb,stats] = anovan(cpd_temp,{cpd_temp_term,cpd_temp_group},'model','interaction','varnames',{'Term','Group'},...
        'display','off');
    if p(3)<0.05;
        c = multcompare(stats,'Dimension',[1 2],'display','off');
        panovaset(itime,:) = c(pickt,6);
    end        
end

psig= [1 2 2 3 3];%, 4 4 5 5 6 6]; % for figure
xlimlist = [51 85; 1 50; 51 85; 1 50; 51 85;...
    1 50; 51 85;1 50; 51 85;1 50; 51 85;];
wlist = (diff(xlimlist')+1)*0.06;
termlist = {'Outcome(t)','Outcome(t-1)1','Outcome(t-1)2','Outcome(t-2)1','Outcome(t-2)2',...
    'Value(t)1','Value(t)2','Value(t-1)1','Value(t-1)2','Value(t-2)1','Value(t-2)2'};
for isig = 1:length(psig)
    if ismember(psig(isig),[1 2 3]); ymax = 0.03;
    elseif ismember(psig(isig),[4 5 6]); ymax = 0.05; end
    f1 = figure;
%     f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 wlist(isig) 2.5]);
    hold on
    stdshade(cpd_dHP(:,:,psig(isig)), 0.3, cmap(1,:))%dorsal %purple
    stdshade(cpd_vHP(:,:,psig(isig)), 0.3, cmap(2,:))%ventral %green
    %     title(termlist{isig});
    plot((panovaset(:,psig(isig))<0.05)-1+0.99*ymax,'color','k','Marker','*','MarkerSize',1,'linestyle','none','linewidth',0.25)

    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    xticks([5 20 35 40 55 80]);xticklabels({});  %xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    xlim([xlimlist(isig,:)]); ylim([0 ymax]); yticks([0 ymax]); yticklabels({});  %xlim([1 120]); ylim([0 0.04])
    hold off
        saveas(f1,['E:\data\vHPC\all\figure\speed\rwprob\nospeed',termlist{isig},d,'.tif'])
%     print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwprob\',termlist{isig},subt,num2str(d),'.ai']);
%     saveas(f1,['E:\data\vHPC\all\figure\cpd\rwprob\',termlist{isig},d,'.tif'])
%     close all

end
%% Reversal - cue period
clear
mouseID = input('mouseID = ');
load('behavior_align.mat','C_z','trial_80R_ind','trial_50R_ind','odorCue','rwProb','rev_b','fr','trial_kind')
C_cue = mean(C_z([1:sum(rev_b), 420-sum(rev_b):end],0.5*fr+1:3.5*fr,:),2);
cueanova = trial_80R_ind([1:sum(rev_b), 420-sum(rev_b):end])+-1*trial_50R_ind([1:sum(rev_b), 420-sum(rev_b):end]); cueanova(cueanova==0)=[];
probanova = (odorCue([1:sum(rev_b), 420-sum(rev_b):end])==find(trial_kind==1)-1).*(rwProb(([1:sum(rev_b), 420-sum(rev_b):end]),find(trial_kind==1)))+...
    (odorCue([1:sum(rev_b), 420-sum(rev_b):end])==find(trial_kind==2)-1).*(rwProb(([1:sum(rev_b), 420-sum(rev_b):end]),find(trial_kind==2))); probanova(probanova==0)=[];


clear panovaset
for ii = 1:size(C_cue,3)
    % [p,tbl,stats] = anovan(C_cue(trial_80R_ind|trial_80P_ind,:,ii),{cueanova, outanova},'varnames',{'Identity','Contingency'},...
    %     'display','off');
    [p,tbl,stats] = anovan(C_cue(trial_80R_ind([1:sum(rev_b), 420-sum(rev_b):end])|trial_50R_ind([1:sum(rev_b), 420-sum(rev_b):end]),:,ii),...
        {cueanova,probanova},'varnames',{'Identity','Contingency'},    'display','off');
    panovaset(ii,:) = p;
end
% endpien = [sum(panovaset(:,3)<0.05),sum(panovaset(:,3)>=0.05&panovaset(:,2)<0.05),...
%     sum(panovaset(:,3)>=0.05&panovaset(:,2)>=0.05&panovaset(:,1)<0.05), ...
%     sum(panovaset(:,3)>=0.05&panovaset(:,2)>=0.05&panovaset(:,1)>=0.05)];
endpien = [sum(panovaset(:,2)>=0.05&panovaset(:,1)<0.05), sum(panovaset(:,2)<0.05&panovaset(:,1)>=0.05),...
    sum(panovaset(:,2)<0.05&panovaset(:,1)<0.05), sum(panovaset(:,2)>=0.05&panovaset(:,1)>=0.05)];
eval(['endpien_',mouseID,' = endpien;']);
save('E:\data\vHPC\all\endpien_rwprob','-regexp','endpien_','-append')
%% figure each mouse
mouselist = {'endpien_dHP03','endpien_dHP04','endpien_dHP06','endpien_dHP07','endpien_dHP08',...
    'endpien_vHP06','endpien_vHP08','endpien_vHP10','endpien_vHP11'};
figurepath = 'E:\data\vHPC\all\figure\';
day = datetime('today');
d = datestr(day,'yymmdd');
for ii = 1:length(mouselist)
f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 15]);
eval(['endpien=',mouselist{1,ii},';'])
pie(endpien)
colormap([65 105 225; 255 99 71;138 43 226;  130 130 130;]./255)
if ~exist('mouseID')
mouseIDtmp = strsplit(mouselist{1,ii},'_');
mouseID = mouseIDtmp{1,2};
end
title([mouseID])
legend('Identity','Contingency','Both','None')
saveas(f_cue,[figurepath,'Contingency_rwprob_',mouseID,'.tif'])
clear mouseID
end
%% figure group
load('E:\data\vHPC\all\endpien_rwprob')
endpien_dHP = sum([endpien_dHP03;endpien_dHP04;endpien_dHP06;endpien_dHP07;endpien_dHP08]);
endpien_vHP = sum([endpien_vHP06;endpien_vHP08;endpien_vHP10;endpien_vHP11]);
mouseG = {'dHP','vHP'};
day = datetime('today');
d = datestr(day,'yymmdd');
for ii = 1:2
    eval(['endpien = endpien_',mouseG{ii},';']);
    figurepath = 'E:\data\vHPC\all\figure\';
    day = datetime('today');
    d = datestr(day,'yymmdd');
    f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 15]);
    pie(endpien)
    colormap([65 105 225; 255 99 71;138 43 226;  130 130 130;]./255)
    title(mouseG{ii})
    legend('Identity','Contingency','Both','None')
    saveas(f_cue,[figurepath,'Contingency_rwprob_',mouseG{ii},'.tif'])    
end
%% glm data
mouselist = {'glmresult_dHP03','glmresult_dHP04','glmresult_dHP06','glmresult_dHP07','glmresult_dHP08',...
    'glmresult_vHP06','glmresult_vHP07','glmresult_vHP08','glmresult_vHP10','glmresult_vHP11'};
figurepath = 'E:\data\vHPC\all\figure\glm_rwprob\';
glmbetaset_dHP = [];
glmbetaset_vHP = [];
glmtvalueset_dHP = [];
glmtvalueset_vHP = [];

for ii = 1:length(mouselist)
    if contains(mouselist{1,ii},'dHP')
        eval(['lengthtmp = length(', mouselist{1,ii},');'])
        eval(['settmp = ', mouselist{1,ii},';'])
        for oo = 1:lengthtmp
            betatmp(:,oo) = settmp{1,oo}.beta;
            ttmp(:,oo) = settmp{1,oo}.t;
        end
        glmbetaset_dHP = [glmbetaset_dHP,betatmp];
        glmtvalueset_dHP = [glmtvalueset_dHP,ttmp];
        
    elseif contains(mouselist{1,ii},'vHP')
        eval(['lengthtmp = length(', mouselist{1,ii},');'])
        eval(['settmp = ', mouselist{1,ii},';'])
        for oo = 1:lengthtmp
            betatmp(:,oo) = settmp{1,oo}.beta;
            ttmp(:,oo) = settmp{1,oo}.t;
        end
        glmbetaset_vHP = [glmbetaset_vHP,betatmp];
        glmtvalueset_vHP = [glmtvalueset_vHP,ttmp];
    end
    clear betatmp ttmp
end

%% glm figure
figurepath = 'E:\data\vHPC\all\figure\glm_rwprob\glm_raw\';
yset = {'glmbetaset_dHP','glmbetaset_vHP','glmtvalueset_dHP','glmtvalueset_vHP'};
ylimset = [0 0.4; 0 0.4; 0 0.2; 0 0.1; 0 0.35;... %beta
           0 3.5;  0 5;   0 3;   0 2.5; 0 5];
% ylimset = [-0.1 0.1; -0.1 0.1;-0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1;... %beta
%     -1 1; -1 1; -1 1; -1 1;   -1 1;   -1 1];   %tvalue

for HP = 1:2
    f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_value = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_rw = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licks = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_lickm = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licke = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    eval(['ytmp=',yset{1,HP},';'])
    nst=1;
    for ii = 1:size(winstep,1)
        if contains(xlist{1,ii},'cue')
            %             cset = [0 0 205; 0 128 255; 80 80 80]./255;
            %             colortmp = cset(contains(xlist{1,ii},'75')+2*contains(xlist{1,ii},'25')+3*contains(xlist{1,ii},'0'),:);
            colortmp = [0 0 205]./255;
            figure(f_cue); hold on; xticks([0 15 30 35]);xticklabels({'0','1.5','3','3.5'});
            ylim(ylimset(floor((HP-1)/2)*6+1,:));
        elseif contains(xlist{1,ii},'Value')
            colortmp = [199 21 133]./255;
            figure(f_value); hold on;
            ylim(ylimset(floor((HP-1)/2)*6+1,:));
        elseif contains(xlist{1,ii},'o_rw')
            cset = [0 0 205; 80 80 80]./255;
            colortmp = cset(contains(xlist{1,ii},'x')+1,:);
            figure(f_rw);hold on; xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*6+2,:));
        elseif contains(xlist{1,ii},'lick')
            colortmp = [0 0 0]./255;
            if contains(xlist{1,ii},'onset'); figure(f_licks);ylim(ylimset(floor((HP-1)/2)*6+3,:));
                xticks([-5 0 10]);xticklabels({'-0.5','0','1'});
            elseif contains(xlist{1,ii},'mid'); figure(f_lickm);ylim(ylimset(floor((HP-1)/2)*6+4,:));
                xticks([0 20 40 60]);xticklabels({'0','2','4','6'});
            elseif contains(xlist{1,ii},'offset'); figure(f_licke);ylim(ylimset(floor((HP-1)/2)*6+5,:));
                xticks([0 20 40 60]);xticklabels({'0','2','4','6'});
            end
            hold on;
        end
        stdshade(abs(ytmp(nst:nst+size(winstep{ii,1},2)-1,:))', 0.3, colortmp,winstep{ii,1})
%         stdshade(ytmp(nst:nst+size(winstep{ii,1},2)-1,:)', 0.3, colortmp,winstep{ii,1})
        xlim([winstep{ii,1}(1,1) winstep{ii,1}(1,end)])
        nst = nst+size(winstep{ii,1},2);
    end
    figure(f_cue); title('Cue'); ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_cue,[figurepath,yset{1,HP},'_cue.tif'])
    figure(f_value); title('Value'); ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_value,[figurepath,yset{1,HP},'_value.tif'])
    figure(f_rw); title('Reward');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([25 25],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_rw,[figurepath,yset{1,HP},'_rw.tif'])
    figure(f_licks); title('Lick onset');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_licks,[figurepath,yset{1,HP},'_licks.tif'])
    figure(f_lickm); title('Lick');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_lickm,[figurepath,yset{1,HP},'_lickm.tif'])
    figure(f_licke); title('Lick offset');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_licke,[figurepath,yset{1,HP},'_licke.tif'])
    close all
end
%% glm figure diff
figurepath = 'E:\data\vHPC\all\figure\glm_rwprob\glm_raw\';
yset = {'glmbetaset_dHP','glmbetaset_vHP','glmtvalueset_dHP','glmtvalueset_vHP'};
% ylimset = [0 0.25; 0 0.5; 0 0.2; 0 0.1; 0 0.35;... %beta
%            0 3.5;  0 5;   0 3;   0 2.5; 0 5];
ylimset = [-0.1 0.1; -0.1 0.1;-0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1;... %beta
    -1 1; -1 1; -1 1; -1 1;   -1 1;   -1 1];   %tvalue


for HP = [1,3]
    f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_value = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_rw = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licks = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_lickm = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licke = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    eval(['dtmp=',yset{1,HP},';'])
    eval(['vtmp=',yset{1,HP+1},';'])
    nst=1;
    glmpvalueset = ones(size(dtmp,1),1);
    for ii = 1:size(winstep,1)
        for jj = 1:size(winstep{ii,1},2)
            [p,h] = ttest2(dtmp(nst+jj-1,:),vtmp(nst+jj-1,:)');
            
            glmpvalueset(nst+jj-1,1) = p
        end
        if contains(xlist{1,ii},'cue')
            colortmp = [0 0 205]./255;
            figure(f_cue); hold on; xticks([0 15 30 35]);xticklabels({'0','1.5','3','3.5'});
            ylim(ylimset(floor((HP-1)/2)*6+1,:));
        elseif contains(xlist{1,ii},'Value')
            colortmp = [199 21 133]./255;
            figure(f_value); hold on;
            ylim(ylimset(floor((HP-1)/2)*6+1,:));
        elseif contains(xlist{1,ii},'o_rw')
            cset = [0 0 205; 80 80 80]./255;
            colortmp = cset(contains(xlist{1,ii},'x')+1,:);
            figure(f_rw);hold on; xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*6+2,:));
        elseif contains(xlist{1,ii},'lick')
            colortmp = [0 0 0]./255;
            if contains(xlist{1,ii},'onset'); figure(f_licks);ylim(ylimset(floor((HP-1)/2)*6+3,:));
                xticks([-5 0 10]);xticklabels({'-0.5','0','1'});
            elseif contains(xlist{1,ii},'mid'); figure(f_lickm);ylim(ylimset(floor((HP-1)/2)*6+4,:));
                xticks([0 20 40 60]);xticklabels({'0','2','4','6'});
            elseif contains(xlist{1,ii},'offset'); figure(f_licke);ylim(ylimset(floor((HP-1)/2)*6+5,:));
                xticks([0 20 40 60]);xticklabels({'0','2','4','6'});
            end
            hold on;
        end
%         stdshade(abs(ytmp(nst:nst+size(winstep{ii,1},2)-1,:))', 0.3, colortmp,winstep{ii,1})
        stdshade(ytmp(nst:nst+size(winstep{ii,1},2)-1,:)', 0.3, colortmp,winstep{ii,1})
        xlim([winstep{ii,1}(1,1) winstep{ii,1}(1,end)])
        nst = nst+size(winstep{ii,1},2);
    end
    figure(f_cue); title('Cue'); ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_cue,[figurepath,yset{1,HP},'_cue.tif'])
    figure(f_value); title('Value'); ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([15 15],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([30 30],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_value,[figurepath,yset{1,HP},'_value.tif'])
    figure(f_rw); title('Reward');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([25 25],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_rw,[figurepath,yset{1,HP},'_rw.tif'])
    figure(f_licks); title('Lick onset');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_licks,[figurepath,yset{1,HP},'_licks.tif'])
    figure(f_lickm); title('Lick');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_lickm,[figurepath,yset{1,HP},'_lickm.tif'])
    figure(f_licke); title('Lick offset');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_licke,[figurepath,yset{1,HP},'_licke.tif'])
    close all
end
%% SRC plot
pvalue_set_dHP = cat(2,pvalue_dHP03,pvalue_dHP04,pvalue_dHP06,pvalue_dHP07,pvalue_dHP08);
src_set_dHP = cat(2,src_dHP03,src_dHP04,src_dHP06,src_dHP07,src_dHP08);
panovaset_dHP = cat(1,panovaset_dHP03,panovaset_dHP04,panovaset_dHP06,panovaset_dHP07,panovaset_dHP08);
pvalue_set_vHP = cat(2,pvalue_vHP06,pvalue_vHP08,pvalue_vHP10,pvalue_vHP11);
src_set_vHP = cat(2,src_vHP06,src_vHP08,src_vHP10,src_vHP11);
panovaset_vHP = cat(1,panovaset_vHP06,panovaset_vHP08,panovaset_vHP10,panovaset_vHP11);
mouseG = {'dHP','vHP'};
colormap = [65 105 225; 255 99 71;138 43 226;]./255;
for ii = 1:2
    eval(['src = src_set_',mouseG{ii},';']);
    eval(['pvalue = pvalue_set_',mouseG{ii},';']);
    eval(['panova = panovaset_',mouseG{ii},';']);
    day = datetime('today');
    d = datestr(day,'yymmdd');
    
    f=figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);
    hold on
    sig_cue = panova(:,1)<0.05;
    sig_contingency = panova(:,2)<0.05;
    sig_non = panova(:,1)>=0.05 & panova(:,2)>=0.05;
    plotset = [src(3,:,1);src(3,:,2)];
    scatter(plotset(1,sig_non),plotset(2,sig_non),10,[0.3 0.3 0.3],'o','filled')
    scatter(plotset(1,sig_cue),plotset(2,sig_cue),20,colormap(1,:),'o','filled')
    scatter(plotset(1,sig_contingency),plotset(2,sig_contingency),30,colormap(2,:),'o','linewidth',1.5)
%     [R ,p]= corrcoef(src(3,:,1),src(3,:,2));
    line([-1 1], [0 0],'Color','k','LineStyle',':')
    line([0 0],[-1 1],'Color','k','LineStyle',':')
    xlabel('Before reversal'); ylabel('After reversal')
    
    saveas(f,[cd,'\figure\SRCrev_rwprob_',mouseG{ii},d,'.tif'])
%     src_pie(:,ii) = [sum(sig_before&sig_after),sum(sig_before&~sig_after),sum(~sig_before&sig_after),sum(~sig_before&~sig_after)];

end
%%
for ii = 1:2
    figurepath = 'E:\data\vHPC\all\figure\';
    day = datetime('today');
    d = datestr(day,'yymmdd');
    f = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
    pie(src_pie(:,ii))
    colormap([138 43 226; 0 0 255 ;255 0 0;130 130 130;]./255)
%     title(mouseG{ii})
    %legend('Both','Before','After','None')
    saveas(f,[figurepath,'SRCpie',mouseG{ii},d,'.tif'])    
end
%% SRC plot, beta delay(t) vs. (t-1)
load('period_reg_rwprob.mat')
beta_set_dHP = cat(2,betaset_dHP03,betaset_dHP04,betaset_dHP06,betaset_dHP07,betaset_dHP08,betaset_dHP10);
pvalue_set_dHP = cat(2,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
src_set_dHP = cat(2,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
beta_set_vHP = cat(2,betaset_vHP06,betaset_vHP07,betaset_vHP08,betaset_vHP11,betaset_vHP12,betaset_vHP14);
pvalue_set_vHP = cat(2,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
src_set_vHP = cat(2,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);
% 1st : coeff 0:a 1:outcome 2:out(t-1) 3:out (t-2) 4:value 5:v(t-1)
% 6:v(t-2) 7:lick 8: speed
% 2nd : neuron
% 3rd : time 1:delay 2:ITI
mouseG = {'dHP','vHP'};
colormap = [65 105 225; 255 99 71;138 43 226;]./255;
plotname = {'V(d,t) V(d,t-1)', 'V(I,t) O(I,t)', 'V(I,t) V(d,t-1)' };
labelname = {'V_{d}(t)','V_{d}(t-1)';'V_{ITI}(t)','O_{ITI}(t)',;'V_{ITI}(t)','V_{d}(t-1)'};
pIdx = [3 1 4 1; 3 2 1 2; 3 2 4 1];
for ip = 1:size(pIdx,1)
    f=figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 12]);
    for ii = 1:2
        eval(['betaset = beta_set_',mouseG{ii},';']);
        eval(['srcset = src_set_',mouseG{ii},';']);
        eval(['pvalueset = pvalue_set_',mouseG{ii},';']);
        day = datetime('today');
        d = datestr(day,'yymmdd');
        subplot(2,2,ii)
        hold on
        scatter(srcset(pIdx(ip,1),:,pIdx(ip,2)),srcset(pIdx(ip,3),:,pIdx(ip,4)),10,[0.3 0.3 0.3],'o','filled')
        line([-1 1], [0 0],'Color','k','LineStyle',':')
        line([0 0],[-1 1],'Color','k','LineStyle',':')
        xlim([-1,1]); ylim([-1 1])
        [r,p]=corrcoef(srcset(pIdx(ip,1),:,pIdx(ip,2)),srcset(pIdx(ip,3),:,pIdx(ip,4)));
        text(-0.8, 0.8, {['r=',num2str(r(1,2),'%.3f')],['p=',num2str(p(1,2),'%.3f')]})
        title([mouseG{ii}])
        xlabel(['Src ',labelname{ip,1}]); ylabel(['Src ',labelname{ip,2}]);
    end
    % resultant vector
    const = [1 2 4];
    numconst = length(const);
    result = cell(4,numconst);
    lgd_vector = {'2\theta','4\theta'};
    colorspec = [1 0 0; 0 0 1];
    idx_con = [2 3]; % plot only 2*theta and 4*theta
    idx_test = [2 1]; % 1: x component, 2: y component
    sig_range = [0.05 0.01 0.001];
    axis_range = [-0.1 0.1];
    axis_tick = [-0.1 0 0.1];
    iset=1;
    numplot = length(idx_con);
    for iHP=1:2
        eval(['srcset = src_set_',mouseG{iHP},';']);
        [theta,rho] = cart2pol(srcset(pIdx(ip,1),:,pIdx(ip,2))',srcset(pIdx(ip,3),:,pIdx(ip,4))');
        for icon = 1:numconst
            tmp_theta = const(icon)*theta;
            tmp_rho = rho;
            % vector of each neurons
            result{iHP,icon} = [real(tmp_rho.*exp(1i*tmp_theta)), imag(tmp_rho.*exp(1i*tmp_theta))];
            % resultant vector of the region
            tmp_rv = sum(tmp_rho.*exp(1i*tmp_theta),1)/length(rho);
            result{iHP,icon} = [result{iHP,icon}; real(tmp_rv) imag(tmp_rv)];
        end
        subplot(2,2,iHP+2)
        for iplot = 1:numplot
            icon = idx_con(iplot);
            iaxis = idx_test(iplot);
            tmp_rv = result{iHP,icon}(end,:); % resultant vector
            hold on
            h(iplot) = plot([0 tmp_rv(1)], [0 tmp_rv(2)],'-','Color',colorspec(iplot,:),'LineWidth',2);
            plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor','w',...
                'MarkerEdgeColor',colorspec(iplot,:),'Color',colorspec(iplot,:),'LineWidth',2);
            
            numcell = size(result{iHP,icon}(1:end-1,:),1);
            [pval,~,stats] = ranksum(result{iHP,icon}(1:end-1,iaxis),zeros(numcell,1));
            idx_sig = find(sig_range > pval);
            if isempty(idx_sig) == 0 % just show p<0.05 or not
                plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor',colorspec(iplot,:),...
                    'MarkerEdgeColor',colorspec(iplot,:));
            end
        end
        plot([0 0],axis_range,'k--');
        plot(axis_range,[0 0],'k--');
        title(mouseG{iHP})
        xlim(axis_range); ylim(axis_range);
        set(gca,'box','off','Tickdir','out','XTick',axis_tick,'XTickLabel',axis_tick,...
            'YTick',axis_tick,'YTickLabel',axis_tick,'LineWidth',1);
        %         xlabel('R\cdotCos'); ylabel('R\cdotSin');
        
        legend(h,lgd_vector,'Box','Off','Position',[0.7 0.2 0.05 0.01],...
            'Orientation','Vertical');
    end
    
    saveas(f,['E:\data\vHPC\all\figure\zscore\src compare ',plotname{1,ip},'.tif'])
end
%% SRC, beta delay vs. ITI
load('period_reg_rwprob.mat')
beta_set_dHP = cat(2,betaset_dHP03,betaset_dHP04,betaset_dHP06,betaset_dHP07,betaset_dHP08,betaset_dHP10);
pvalue_set_dHP = cat(2,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
src_set_dHP = cat(2,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
beta_set_vHP = cat(2,betaset_vHP06,betaset_vHP07,betaset_vHP08,betaset_vHP11,betaset_vHP12,betaset_vHP14);
pvalue_set_vHP = cat(2,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
src_set_vHP = cat(2,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);
% 1st : (intercept) o ot-1 ot-2 v vt-1 vt-2 lick speed
% 2nd : Cell
% 3rd : delay outcome ITI
termlist = [4; 1];

mouseG = {'dHP','vHP'};
colormap = [65 105 225; 255 99 71;138 43 226;]./255;
day = datetime('today');
d = datestr(day,'yymmdd');
f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 10]);
for iiHP = 1:2
    eval(['betaset = beta_set_',mouseG{iiHP},';']);
    eval(['srcset = src_set_',mouseG{iiHP},';']);
    eval(['pvalueset = pvalue_set_',mouseG{iiHP},';']);
    for iterm = 1:2
        subplot(3,2,2*(iiHP-1)+iterm)        
        axis square
        neuronset = [pvalueset(5,:,1)>=0.05&pvalueset(termlist(iterm,1)+1,:,2)>=0.05;...
            pvalueset(5,:,1)<0.05&pvalueset(termlist(iterm,1)+1,:,2)>=0.05;...
            pvalueset(5,:,1)>=0.05&pvalueset(termlist(iterm,1)+1,:,2)<0.05;...
            pvalueset(5,:,1)<0.05&pvalueset(termlist(iterm,1)+1,:,2)<0.05];
        hold on
        scatter(srcset(4,:,1),srcset(termlist(iterm,1),:,2),1,[0.5 0.5 0.5],'o','filled') %non sig
        scatter(srcset(4,neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,2) >0),1),...
            srcset(termlist(iterm,1),neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,2) >0),2),1,[0.58 0 0.83],'o','filled') % sig 1,3
        scatter(srcset(4,neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,2) <0),1),...
            srcset(termlist(iterm,1),neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,2) <0),2),1,[0.85 0.44 0.85],'o','filled') % sig 2,4
        line([-1 1], [0 0],'Color','k','LineStyle',':')
        line([0 0],[-1 1],'Color','k','LineStyle',':')
        xlim([-1,1]); ylim([-1 1])
        [r,p]=corrcoef(srcset(4,:,1),srcset(termlist(iterm,1),:,2));
        pset(iiHP,iterm) = p(1,2);
        if p(1,2)<0.05; f = @(x) r(1,2)*x; ezplot( f, [-1 ,1, -1,1]); end
        text(0.5, 0.5,...
            {['R = ',num2str(round(r(1,2),4))],['p = ',num2str(round(p(1,2),4))]},...
            'Fontsize', 6);
        xlabel([]);
        title([mouseG{iiHP},' Delay - Outcome'],'fontsize',7)
        
        barset(iiHP,2*(iterm-1)+1:2*(iterm-1)+2) = ...
            [mean(neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,2) >0)),mean(neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,2) <0))];
        
%         subplot(3,2,iterm+2)
%         axis equal
%         neuronset = [pvalueset(5,:,1)>=0.05&pvalueset(termlist(iterm,1)+1,:,3)>=0.05;...
%             pvalueset(5,:,1)<0.05&pvalueset(termlist(iterm,1)+1,:,3)>=0.05;...
%             pvalueset(5,:,1)>=0.05&pvalueset(termlist(iterm,1)+1,:,3)<0.05;...
%             pvalueset(5,:,1)<0.05&pvalueset(termlist(iterm,1)+1,:,3)<0.05];
%         hold on
%         scatter(srcset(4,:,1),srcset(termlist(iterm,1),:,3),1,[0.5 0.5 0.5],'o','filled')
%         scatter(srcset(4,neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,3) >0),1),...
%             srcset(termlist(iterm,1),neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,3) >0),3),1,[1 0.27 0],'o','filled') % sig 1,3
%         scatter(srcset(4,neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,3) <0),1),...
%             srcset(termlist(iterm,1),neuronset(4,:) & (srcset(4,:,1) .* srcset(termlist(iterm,1),:,3) <0),3),1,[1 0.85 0],'o','filled') % sig 2,4
%         line([-1 1], [0 0],'Color','k','LineStyle',':')
%         line([0 0],[-1 1],'Color','k','LineStyle',':')
%         xlim([-1,1]); ylim([-1 1]);
%         [r,p]=corrcoef(srcset(4,:,1),srcset(termlist(iterm,1),:,3));
%         if p(1,2)<0.05; f = @(x) r(1,2)*x; ezplot( f, [-1 ,1, -1,1]); end
%         text(0.5, 0.5,...
%             {['R = ',num2str(round(r(1,2),4))],['p = ',num2str(round(p(1,2),4))]},...
%             'Fontsize', 6);
%         xlabel([]);
%         title([mouseG{iiHP},' Delay - ITI'],'fontsize',7)
%         
%         subplot(3,2,iterm+4)
%         axis equal
%         neuronset = [pvalueset(5,:,2)>=0.05&pvalueset(termlist(iterm,1)+1,:,3)>=0.05;...
%             pvalueset(5,:,2)<0.05&pvalueset(termlist(iterm,1)+1,:,3)>=0.05;...
%             pvalueset(5,:,2)>=0.05&pvalueset(termlist(iterm,1)+1,:,3)<0.05;...
%             pvalueset(5,:,2)<0.05&pvalueset(termlist(iterm,1)+1,:,3)<0.05];
%         hold on
%         scatter(srcset(4,:,2),srcset(termlist(iterm,1),:,3),1,[0.5 0.5 0.5],'o','filled')
%         scatter(srcset(4,neuronset(4,:) & (srcset(4,:,2) .* srcset(termlist(iterm,1),:,3) >0),2),...
%             srcset(termlist(iterm,1),neuronset(4,:) & (srcset(4,:,2) .* srcset(termlist(iterm,1),:,3) >0),3),1,[1 0.27 0],'o','filled') % sig 1,3
%         scatter(srcset(4,neuronset(4,:) & (srcset(4,:,2) .* srcset(termlist(iterm,1),:,3) <0),2),...
%             srcset(termlist(iterm,1),neuronset(4,:) & (srcset(4,:,2) .* srcset(termlist(iterm,1),:,3) <0),3),1,[1 0.85 0],'o','filled') % sig 2,4
%         line([-1 1], [0 0],'Color','k','LineStyle',':')
%         line([0 0],[-1 1],'Color','k','LineStyle',':')
%         xlim([-1,1]); ylim([-1 1])
%         [r,p]=corrcoef(srcset(4,:,2),srcset(termlist(iterm,1),:,3));
%         if p(1,2)<0.05; f = @(x) r(1,2)*x; ezplot( f, [-1 ,1, -1,1]); end
%         text(0.5, 0.5,...
%             {['R = ',num2str(round(r(1,2),4))],['p = ',num2str(round(p(1,2),4))]},...
%             'Fontsize', 6);
%         xlabel([]);
%         title([mouseG{iiHP},' Outcome - ITI'],'fontsize',7)      
    end
    
end
subplot(3,2,5)
axis square
hold on
bar(0.5,barset(1,1),'facecolor',[0.58 0 0.83],'Barwidth',0.2)
bar(0.75,barset(1,2),'facecolor',[0.85 0.44 0.85],'Barwidth',0.2)
bar(1.25,barset(2,1),'facecolor',[0.58 0 0.83],'Barwidth',0.2)
bar(1.5,barset(2,2),'facecolor',[0.85 0.44 0.85],'Barwidth',0.2)
xlim([0.2 1.8]); xticklabels({});

subplot(3,2,6)
axis square
hold on
bar(0.5,barset(1,3),'facecolor',[0.58 0 0.83],'Barwidth',0.2)
bar(0.75,barset(1,4),'facecolor',[0.85 0.44 0.85],'Barwidth',0.2)
bar(1.25,barset(2,3),'facecolor',[0.58 0 0.83],'Barwidth',0.2)
bar(1.5,barset(2,4),'facecolor',[0.85 0.44 0.85],'Barwidth',0.2)
xlim([0.2 1.8]); xticklabels({});
%     saveas(f1,['E:\data\vHPC\all\figure\zscore\src compare delay_value ITI_outcome_',mouseG{iiHP},d,'.tif'])
print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\zscore\src compare delay_value ITI_outcome_',d,'.ai']);

% % resultant vector
% const = [1 2 4];
% numconst = length(const);
% result = cell(4,numconst);
% lgd_vector = {'2\theta','4\theta'};
% colorspec = [1 0 0; 0 0 1];
% idx_con = [2 3]; % plot only 2*theta and 4*theta
% idx_test = [2 1]; % 1: x component, 2: y component
% sig_range = [0.05 0.01 0.001];
% axis_range = [-0.1 0.1];
% axis_tick = [-0.1 0 0.1];
% numplot = length(idx_con);
% f=figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 12]);
% for iHP=1:2
%     eval(['srcset = src_set_',mouseG{iHP},';']);
%     for iset=1:2
%         if iset==1;
%             [theta,rho] = cart2pol(srcset(3,:,1)',srcset(3,:,2)'); title_t='V(d)-V(ITI)';% nx1 double
%         elseif iset==2;
%             [theta,rho] = cart2pol(srcset(3,:,1)',srcset(1,:,2)'); title_t='V(d)-O(ITI)';% nx1 double
%         end
%         for icon = 1:numconst
%             tmp_theta = const(icon)*theta;
%             tmp_rho = rho;
%             % vector of each neurons
%             result{iHP+2*(iset-1),icon} = [real(tmp_rho.*exp(1i*tmp_theta)), imag(tmp_rho.*exp(1i*tmp_theta))];
%             % resultant vector of the region
%             tmp_rv = sum(tmp_rho.*exp(1i*tmp_theta),1)/length(rho);
%             result{iHP+2*(iset-1),icon} = [result{iHP+2*(iset-1),icon}; real(tmp_rv) imag(tmp_rv)];
%         end
%         
%         subplot(2,2,iHP+2*(iset-1))
%         for iplot = 1:numplot
%             icon = idx_con(iplot);
%             iaxis = idx_test(iplot);
%             tmp_rv = result{iHP+2*(iset-1),icon}(end,:); % resultant vector
%             hold on
%             h(iplot) = plot([0 tmp_rv(1)], [0 tmp_rv(2)],'-','Color',colorspec(iplot,:),'LineWidth',2);
%             plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor','w',...
%                 'MarkerEdgeColor',colorspec(iplot,:),'Color',colorspec(iplot,:),'LineWidth',2);
%             
%             numcell = size(result{iHP+2*(iset-1),icon}(1:end-1,:),1);
%             [pval,~,stats] = ranksum(result{iHP+2*(iset-1),icon}(1:end-1,iaxis),zeros(numcell,1));
%             idx_sig = find(sig_range > pval);
%             if isempty(idx_sig) == 0 % just show p<0.05 or not
%                 plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor',colorspec(iplot,:),...
%                     'MarkerEdgeColor',colorspec(iplot,:));
%             end
%         end
%         plot([0 0],axis_range,'k--');
%         plot(axis_range,[0 0],'k--');
%         title([mouseG{iHP},' ',title_t])
%         xlim(axis_range); ylim(axis_range);
%         set(gca,'box','off','Tickdir','out','XTick',axis_tick,'XTickLabel',axis_tick,...
%             'YTick',axis_tick,'YTickLabel',axis_tick,'LineWidth',1);
% %         xlabel('R\cdotCos'); ylabel('R\cdotSin');
%         
%         legend(h,lgd_vector,'Box','Off','Position',[0.7 0.2 0.05 0.01],...
%             'Orientation','Vertical');
%     end
% end
% saveas(f,['E:\data\vHPC\all\figure\zscore\src compare delay_valye ITI_outcome vector',d,'.tif'])

%% RPE
%data settting
load('RPE_rwprob.mat')
% load('RPE_rwprob_trialcontrol.mat');
load('lick_neuron.mat')
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
betaset_dHP = cat(2,betaset_RPE_dHP03,betaset_RPE_dHP04,betaset_RPE_dHP06,betaset_RPE_dHP07,betaset_RPE_dHP08,betaset_RPE_dHP10);%
pvalueset_dHP = cat(2,pvalueset_RPE_dHP03,pvalueset_RPE_dHP04,pvalueset_RPE_dHP06,pvalueset_RPE_dHP07,pvalueset_RPE_dHP08,pvalueset_RPE_dHP10);%
srcset_dHP = cat(2,srcset_RPE_dHP03,srcset_RPE_dHP04,srcset_RPE_dHP06,srcset_RPE_dHP07,srcset_RPE_dHP08,srcset_RPE_dHP10);%
betaset_vHP = cat(2,betaset_RPE_vHP06,betaset_RPE_vHP07,betaset_RPE_vHP08,betaset_RPE_vHP11,betaset_RPE_vHP12,betaset_RPE_vHP14);%
pvalueset_vHP = cat(2,pvalueset_RPE_vHP06,pvalueset_RPE_vHP07,pvalueset_RPE_vHP08,pvalueset_RPE_vHP11,pvalueset_RPE_vHP12,pvalueset_RPE_vHP14);%
srcset_vHP = cat(2,srcset_RPE_vHP06,srcset_RPE_vHP07,srcset_RPE_vHP08,srcset_RPE_vHP11,srcset_RPE_vHP12,srcset_RPE_vHP14);%
mouseG = {'dHP','vHP'};
% 1st : (intercept) o ot-1 ot-2 v vt-1 vt-2 lick speed
% 2nd : neuron
% 3rd : ITI, Outcome 1.5s, Outcome all
% 4th : All, Reward, Unreward
day = datetime('today'); 
d = datestr(day,'yymmdd');
f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 14 10]);
siglevel = 0.05;
%figure
for iHP = 1:2;    
    eval(['pvalueset_RPE = pvalueset_',mouseG{iHP},';']);
    eval(['srcset_RPE = srcset_',mouseG{iHP},';']);
    eval(['lick_set = lick_',mouseG{iHP},';']);
    % rw value both sig in normal regression
    % (3rd 1: ITI, 2:Outcome 1.5s,)
%     pvalueset_RPE = pvalueset_RPE(:,lick_set(:,1)>=0.05,:,:);
%     srcset_RPE  = srcset_RPE(:,lick_set(:,1)>=0.05,:,:);
    %rwtrial
%     subplot(2,3,1+3*(iHP-1))
    subplot(2,2,iHP)
    axis square
    RPE_nomi_r=pvalueset_RPE(2,:,2,1)<siglevel & pvalueset_RPE(5,:,2,2)<siglevel;
    quart1=srcset_RPE(1,:,2,1)>0 & srcset_RPE(4,:,2,2)>0;
    quart2=srcset_RPE(1,:,2,1)<0 & srcset_RPE(4,:,2,2)>0;
    quart3=srcset_RPE(1,:,2,1)<0 & srcset_RPE(4,:,2,2)<0;
    quart4=srcset_RPE(1,:,2,1)>0 & srcset_RPE(4,:,2,2)<0;
    quart{1,iHP}=[quart1;quart2;quart3;quart4];
    hold on
    scatter(srcset_RPE(1,:,2,1),srcset_RPE(4,:,2,2),3,[153 153 153]/255,'o','filled')
    scatter(srcset_RPE(1,RPE_nomi_r&(quart1|quart3),2,1),srcset_RPE(4,RPE_nomi_r&(quart1|quart3),2,2),3,[255 8 0]/255,'o','filled') %1,3
    scatter(srcset_RPE(1,RPE_nomi_r&(quart2|quart4),2,1),srcset_RPE(4,RPE_nomi_r&(quart2|quart4),2,2),3,[16 52 166]/255,'o','filled') %2,4
    line([-1 1], [0 0],'Color','k','LineStyle',':')
    line([0 0],[-1 1],'Color','k','LineStyle',':')
    fit_set{iHP,1} = [srcset_RPE(1,:,2,1);srcset_RPE(4,:,2,2)];
    fit_set{iHP,2} = [srcset_RPE(1,RPE_nomi_r&(quart1|quart3),2,1),srcset_RPE(1,RPE_nomi_r&(quart2|quart4),2,1);...
        srcset_RPE(4,RPE_nomi_r&(quart1|quart3),2,2),srcset_RPE(4,RPE_nomi_r&(quart2|quart4),2,2)];
%     title([mouseG{iHP},' Rewarded'],'fontsize',7)
%     xlabel('SRC for Reward'); ylabel('SRC for Value');
    set(gca,'FontSize',7)
    xticks([-1 0 1]);yticks([-1 0 1]);
    
    %unrwtrial    
%     subplot(2,3,2+3*(iHP-1))
    subplot(2,2,iHP+2)
    axis square
    RPE_nomi_ur=pvalueset_RPE(2,:,2,1)<siglevel & pvalueset_RPE(5,:,2,3)<siglevel;
    quart1=srcset_RPE(1,:,2,1)>0 & srcset_RPE(4,:,2,3)>0;
    quart2=srcset_RPE(1,:,2,1)<0 & srcset_RPE(4,:,2,3)>0;
    quart3=srcset_RPE(1,:,2,1)<0 & srcset_RPE(4,:,2,3)<0;
    quart4=srcset_RPE(1,:,2,1)>0 & srcset_RPE(4,:,2,3)<0;
    quart{2,iHP}=[quart1;quart2;quart3;quart4];
    hold on
    scatter(srcset_RPE(1,:,2,1),srcset_RPE(4,:,2,3),3,[153 153 153]/255,'o','filled')
    scatter(srcset_RPE(1,RPE_nomi_ur&(quart1|quart3),2,1),srcset_RPE(4,RPE_nomi_ur&(quart1|quart3),2,3),3,[214 99 154]/255,'o','filled')
    scatter(srcset_RPE(1,RPE_nomi_ur&(quart2|quart4),2,1),srcset_RPE(4,RPE_nomi_ur&(quart2|quart4),2,3),3,[98 127 191]/255,'o','filled')
    line([-1 1], [0 0],'Color','k','LineStyle',':')
    line([0 0],[-1 1],'Color','k','LineStyle',':')
    %     title([mouseG{iHP},' Unrewarded'])
    xticks([-1 0 1]);yticks([-1 0 1]);
    %     xlabel('SRC for Reward'); ylabel('SRC for Value');
    set(gca,'FontSize',7)
    fit_set{iHP+2,1} = [srcset_RPE(1,:,2,1);srcset_RPE(4,:,2,3)];
    fit_set{iHP+2,2} = [srcset_RPE(1,RPE_nomi_ur&(quart1|quart3),2,1),srcset_RPE(1,RPE_nomi_ur&(quart2|quart4),2,1);...
        srcset_RPE(4,RPE_nomi_ur&(quart1|quart3),2,3),srcset_RPE(4,RPE_nomi_ur&(quart2|quart4),2,3)];
    
    RPE_nomi{1,iHP}=[RPE_nomi_r; RPE_nomi_ur];
    upV_neu{1,iHP} = RPE_nomi_ur&(quart1|quart3);
    RPE_neu{1,iHP} = RPE_nomi_ur&(quart2|quart4);
end
%% Regression RPE
% f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 14 12]);
titlename = {'dHP Rw','vHP Rw','dHP noRw','vHP noRw'};
cmap_s = [128 128 128; 200 120 120]./255;
cmap_f = [0 0 0; 255 0 0]./255;
for ifig = 1:2
    for ii = 1:4
        mdl = fitlm(fit_set{ii,ifig}(1,:),fit_set{ii,ifig}(2,:));
        alpha = table2array(mdl.Coefficients(2,1)); beta = table2array(mdl.Coefficients(1,1));
        p =table2array( mdl.Coefficients(2,4));
        subplot(2,2,ii)
        hold on
%         scatter(fit_set{ii,ifig}(1,:),fit_set{ii,ifig}(2,:),3,...
%             'markerfacecolor',cmap_s(ifig,:),'markeredgecolor',cmap_s(ifig,:))
        fplot(@(x) alpha*x+beta, [-0.8 0.8], 'color',cmap_s(ifig,:))
        xlim([-1 1]); ylim([-1 1]);
        text(-0.9,1.2-ifig*0.3, {['b = ', num2str(alpha)]; ['p = ', num2str(p)]},...
            'Color',cmap_f(ifig,:),'fontsize',7)
%         title(titlename{ii})
    end
        
end
print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\RPE\rwprob_',d,'.ai']);

%% fraction of neuron of RPE
totNeu = [size(RPE_nomi{1,1},2),size(RPE_nomi{1,2},2)]; % total neuron
% RPE_nomi: dHP, vHP ; reward, unreward nomi
% quart: {(Rw, uRw),(dHP, vHP)};  1st,2nd,3rd,4th quart

bardataL = [sum(RPE_nomi{1,1}(1,:)&(quart{1,1}(1,:)|quart{1,1}(3,:))),sum(RPE_nomi{1,1}(2,:)&(quart{2,1}(1,:)|quart{2,1}(3,:))),...
    sum(RPE_nomi{1,2}(1,:)&(quart{1,2}(1,:)|quart{1,2}(3,:))),sum(RPE_nomi{1,2}(2,:)&(quart{2,2}(1,:)|quart{2,2}(3,:)));...
    sum(RPE_nomi{1,1}(1,:)&(quart{1,1}(2,:)|quart{1,1}(4,:))),sum(RPE_nomi{1,1}(2,:)&(quart{2,1}(2,:)|quart{2,1}(4,:))),...
    sum(RPE_nomi{1,2}(1,:)&(quart{1,2}(2,:)|quart{1,2}(4,:))),sum(RPE_nomi{1,2}(2,:)&(quart{2,2}(2,:)|quart{2,2}(4,:)))];
subplot(2,3,3) % Value update signal
axis square
hold on
bar(0.6,bardataL(1,1)./totNeu(1),'facecolor',[255 8 0]./255,'Barwidth',0.2)
bar(0.85,bardataL(1,2)./totNeu(1),'facecolor',[214 99 154]./255,'Barwidth',0.2)
bar(1.25,bardataL(1,3)./totNeu(2),'facecolor',[255 8 0]./255,'Barwidth',0.2)
bar(1.5,bardataL(1,4)./totNeu(2),'facecolor',[214 99 154]./255,'Barwidth',0.2)
xlim([0.3 1.8]); xticks([]); ylim([0 0.1]); yticks([0 0.1]);
set(gca,'FontSize',7)
% ylabel('Fraction of neuron')


subplot(2,3,6) % RPE
axis square
hold on
bar(0.6,bardataL(2,1)./totNeu(1),'facecolor',[16 52 166]./255,'Barwidth',0.2)
bar(0.85,bardataL(2,2)./totNeu(1),'facecolor',[98 127 191]./255,'Barwidth',0.2)
bar(1.25,bardataL(2,3)./totNeu(2),'facecolor',[16 52 166]./255,'Barwidth',0.2)
bar(1.5,bardataL(2,4)./totNeu(2),'facecolor',[98 127 191]./255,'Barwidth',0.2)
xlim([0.3 1.8]); xticks([]); ylim([0 0.1]); yticks([0 0.1]);
set(gca,'FontSize',7)
% ylabel('Fraction of neuron')

% text(3,0.1-iHP*0.02,{[mouseG{1,iHP}],['Total Neuron:', num2str(totNeu), '  Sig Neu:',num2str(sum(sum(bardataL(1:2,2:3))))],...
%     ['+RPE:', num2str(sum(bardataL(1,2:3))),'  -RPE:',num2str(sum(bardataL(2,2:3)))]})

% xticks([1 2 4 5 7 8]); xticklabels({'dHP-sig','vHP-sig','dHP-pos','dHP-neg','vHP-pos','vHP-neg'}); xlim([0.3 8.8])



% saveas(f1,['E:\data\vHPC\all\figure\RPE\rwprob_sig0.1_',d,'.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\RPE\rwprob_',d,'.ai']);




%% chance level
TN = [649 1365 1732 282 649 1483 1816 316; 139 154 261 32 139 206 291 54];
TN_fon = TN./[3304;724];

Chancelevel = [TN_fon(:,1).*TN_fon(:,2), TN_fon(:,5).*TN_fon(:,6)];
testset = [TN(:,4),TN(:,8)];
AN= [3304 3304; 724 724];

for ii = 1:2
    p = binocdf(testset(1,ii),AN(1,ii), Chancelevel(1,ii),'upper');
    pset(1,ii) = p;
    p = binocdf(testset(2,ii),AN(2,ii), Chancelevel(2,ii),'upper');
    pset(2,ii) = p;
    
    
end
%% resultant vector - all neuron
% for resultant vector
const = [1 2 4];
numconst = length(const);
result = cell(4,numconst);
lgd_vector = {'2\theta','4\theta'};
colorspec = [1 0 0; 0 0 1];
idx_con = [2 ]; % plot only 2*theta and 4*theta
idx_test = [2 1]; % 1: x component, 2: y component
sig_range = [0.05 0.01 0.001];
axis_range = [-0.1 0.1];
axis_tick = [-0.1 0 0.1];
numplot = length(idx_con);
f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 10]);
for iHP = 1:2
    eval(['srcset_RPE = srcset_',mouseG{iHP},';']);
    srcset_RPE  = srcset_RPE(:,lick_set(:,1)>=0.05,:,:);
    for iset=1:2
        if iset==1;
            [theta,rho] = cart2pol(srcset_RPE(1,:,2,1)',srcset_RPE(4,:,2,2)'); title_t='Rewarded';% nx1 double
        elseif iset==2;
            [theta,rho] = cart2pol(srcset_RPE(1,:,2,1)',srcset_RPE(4,:,2,3)'); title_t='Unrewarded';% nx1 double
        end
        for icon = 1:numconst
            tmp_theta = const(icon)*theta;
            tmp_rho = rho;
            % vector of each neurons
            result{iHP+2*(iset-1),icon} = [real(tmp_rho.*exp(1i*tmp_theta)), imag(tmp_rho.*exp(1i*tmp_theta))];
            % resultant vector of the region
            tmp_rv = sum(tmp_rho.*exp(1i*tmp_theta),1)/length(rho);
            result{iHP+2*(iset-1),icon} = [result{iHP+2*(iset-1),icon}; real(tmp_rv) imag(tmp_rv)];
        end
        
        subplot(2,2,2*(iHP-1)+iset)
        for iplot = 1:numplot
            icon = idx_con(iplot);
            iaxis = idx_test(iplot);
            tmp_rv = result{iHP+2*(iset-1),icon}(end,:); % resultant vector
            hold on
            h(iplot) = plot([0 tmp_rv(1)], [0 tmp_rv(2)],'-','Color',colorspec(iplot,:),'LineWidth',2);
            plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor','w',...
                'MarkerEdgeColor',colorspec(iplot,:),'Color',colorspec(iplot,:),'LineWidth',2);
            
            numcell = size(result{iHP+2*(iset-1),icon}(1:end-1,:),1);
            [pval,~,stats] = ranksum(result{iHP+2*(iset-1),icon}(1:end-1,iaxis),zeros(numcell,1));
            idx_sig = find(sig_range > pval);
            if isempty(idx_sig) == 0 % just show p<0.05 or not
                plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor',colorspec(iplot,:),...
                    'MarkerEdgeColor',colorspec(iplot,:));
            end
        end
        plot([0 0],axis_range,'k--');
        plot(axis_range,[0 0],'k--');
        title([mouseG{iHP},' ',title_t])
        xlim(axis_range); ylim(axis_range);
        set(gca,'box','off','Tickdir','out','XTick',axis_tick,'XTickLabel',axis_tick,...
            'YTick',axis_tick,'YTickLabel',axis_tick,'LineWidth',1);
                xlabel('R\cdotCos'); ylabel('R\cdotSin');
        
        legend(h,lgd_vector,'Box','Off','Position',[0.7 0.2 0.05 0.01],...
            'Orientation','Vertical');
    end
end
% saveas(f1,['E:\data\vHPC\all\figure\RPE\rwprob_rvall_sig0.1_nolick_',d,'.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\RPE\rwprob_lick_control_',d,'.ai']);
%% resultant vector - each RPE,upV
% for resultant vector
const = [1 2];
numconst = length(const);
result = cell(4,numconst);
lgd_vector = {'dCA1','vCA1'};
colorspec = [129 0 129; 1 127 1]./255;
idx_con = [2]; % plot only 2*theta and 4*theta
idx_test = [2 1]; % 1: x component, 2: y component
sig_range = [0.05 0.01 0.001];
axis_range = [-0.5 0.5];
axis_tick = [-0.2 0 0.2];
numplot = length(idx_con);
f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 5]);
for iHP = 1:2
    eval(['srcset_RPE = srcset_',mouseG{iHP},';']);
    eval(['lick_set = lick_',mouseG{iHP},';']);
%     srcset_RPE  = srcset_RPE(:,lick_set(:,1)>=0.05,:,:);
    quart_tmp1 = quart{1,iHP};
    quart_tmp2 = quart{2,iHP};
    RPEnomi_tmp = RPE_nomi{1,iHP};
    for iset=1:4
%         if iset==1; % Positive Upv
%             selN = (RPEnomi_tmp(1,:)&sum(quart_tmp1)==1);            
%             [theta,rho] = cart2pol(srcset_RPE(1,selN,2,1)',srcset_RPE(4,selN,2,2)'); title_t='Rewarded';% nx1 double
%         elseif iset==2; % Negavite Upv
%             selN = (RPEnomi_tmp(2,:)&sum(quart_tmp2)==1);           
%             [theta,rho] = cart2pol(srcset_RPE(1,selN,2,1)',srcset_RPE(4,selN,2,3)'); title_t='Unrewarded';% nx1 double
        if iset==1; % Positive Upv
            selN = RPEnomi_tmp(1,:)&(quart_tmp1(1,:)|quart_tmp1(3,:));
            [theta,rho] = cart2pol(srcset_RPE(1,selN,2,1)',srcset_RPE(4,selN,2,2)'); title_t='Positive Upv';% nx1 double
        elseif iset==2; % Negavite Upv
            selN = RPEnomi_tmp(2,:)&(quart_tmp2(1,:)|quart_tmp2(3,:));
            [theta,rho] = cart2pol(srcset_RPE(1,selN,2,1)',srcset_RPE(4,selN,2,3)'); title_t='Negavite Upv';% nx1 double
        elseif iset==3; %Positive RPE
            selN = RPEnomi_tmp(1,:)&(quart_tmp1(2,:)|quart_tmp1(4,:));
            [theta,rho] = cart2pol(srcset_RPE(1,selN,2,1)',srcset_RPE(4,selN,2,2)'); title_t='Positive RPE';% nx1 double
        elseif iset==4; %Negative RPE
            selN = RPEnomi_tmp(2,:)&(quart_tmp2(2,:)|quart_tmp2(4,:));
            [theta,rho] = cart2pol(srcset_RPE(1,selN,2,1)',srcset_RPE(4,selN,2,3)'); title_t='Negative RPE';% nx1 double
        end
        
        for icon = 1:numconst
            tmp_theta = const(icon)*theta;
            tmp_rho = rho;
            % vector of each neurons
            result{iHP+2*(iset-1),icon} = [real(tmp_rho.*exp(1i*tmp_theta)), imag(tmp_rho.*exp(1i*tmp_theta))];
            % resultant vector of the region
            tmp_rv = sum(tmp_rho.*exp(1i*tmp_theta),1)/length(rho);
            result{iHP+2*(iset-1),icon} = [result{iHP+2*(iset-1),icon}; real(tmp_rv) imag(tmp_rv)];
        end
        
        subplot(2,2,iset)
        for iplot = 1:numplot
            icon = idx_con(iplot);
            iaxis = idx_test(iplot);
            tmp_rv = result{iHP+2*(iset-1),icon}(end,:); % resultant vector
            hold on
            h(iHP) = plot([0 tmp_rv(1)], [0 tmp_rv(2)],'-','Color',colorspec(iHP,:),'LineWidth',2);
            plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor','w',...
                'MarkerEdgeColor',colorspec(iHP,:),'Color',colorspec(iHP,:),'LineWidth',2);
            
            numcell = size(result{iHP+2*(iset-1),icon}(1:end-1,:),1);
            [pval,~,stats] = ranksum(result{iHP+2*(iset-1),icon}(1:end-1,iaxis),zeros(numcell,1));
            idx_sig = find(sig_range > pval);
            if isempty(idx_sig) == 0 % just show p<0.05 or not
                plot(tmp_rv(1),tmp_rv(2),'o','MarkerSize',6,'MarkerFaceColor',colorspec(iHP,:),...
                    'MarkerEdgeColor',colorspec(iHP,:));
            end
        end
        plot([0 0],axis_range,'k--');
        plot(axis_range,[0 0],'k--');
        title(title_t)
        xlim(axis_range); ylim(axis_range);
        set(gca,'box','off','Tickdir','out','XTick',axis_tick,'XTickLabel',axis_tick,...
            'YTick',axis_tick,'YTickLabel',axis_tick,'LineWidth',1);
                xlabel('R\cdotCos'); ylabel('R\cdotSin');
        

    end
end
        legend(h,lgd_vector,'Box','Off','Position',[0.8 0.8 0.05 0.01],...
            'Orientation','Vertical');
% saveas(f1,['E:\data\vHPC\all\figure\RPE\rwprob_rv_allselN_',d,'.tif'])
%% RPE z-firing
% load('RPE_rwprob_forz.mat')
zzset_dHP = cat(2, zzset_dHP03,zzset_dHP04,zzset_dHP06,zzset_dHP07,zzset_dHP08,zzset_dHP10);
zzset_vHP = cat(2, zzset_vHP06,zzset_vHP07,zzset_vHP08,zzset_vHP11,zzset_vHP12,zzset_vHP13);
% 1st : timeline(30Hz)
% 2nd : cell
% 3rd : trial
%    1:reward 2:unreward 3:75R 4:25R 5:0R 6:75U 7:25U 8:0U
cmap = [128 0 128; 30 45 232; 120 177 255; 125 125 125;]./255;
% RPE_nomi 1cell:dorsal 2cell:ventral/ 1st: Reward, Unreward/ 2nd: neuron
% quart 1cell:dorsal 2cell:ventral/ 1st: Rew-positive,Rew-negative,UnR-positive,UnR-negative/ 2nd: neuron
lines = {'-',':'};
mouseG = {'dHP','vHP'}; fr=30;
for iHP = 1:2;
    f=figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 15]);
    eval(['zzset_RPE = zzset_',mouseG{iHP},';']);
    for iRPE = 1:2 % p or n
        for iqur = 1:2 % 2 or 4
            subplot(2,2,iqur+2*(iRPE-1))
            hold on
            stdshade(zzset_RPE(4.5*fr:11*fr,quart{1,iHP}(iqur,:)& RPE_nomi{1,iHP}(iRPE,:),3-iRPE  )',0.3, cmap(1,:), -1:1/30:5.5,lines{3-iRPE})
            stdshade(zzset_RPE(4.5*fr:11*fr,quart{1,iHP}(iqur,:)& RPE_nomi{1,iHP}(iRPE,:),iRPE*3  )',0.3, cmap(2,:), -1:1/30:5.5,lines{iRPE})
            stdshade(zzset_RPE(4.5*fr:11*fr,quart{1,iHP}(iqur,:)& RPE_nomi{1,iHP}(iRPE,:),iRPE*3+1)',0.3, cmap(3,:), -1:1/30:5.5,lines{iRPE})
            if iRPE==2;
            stdshade(zzset_RPE(4.5*fr:11*fr,quart{1,iHP}(iqur,:)& RPE_nomi{1,iHP}(iRPE,:),iRPE*3+2)',0.3, cmap(4,:), -1:1/30:5.5,lines{iRPE})
            end
            line([0 0],[-0.5 1.5],'Color','k','LineStyle',':')
            line([2.5 2.5],[-0.5 1.5],'Color','k','LineStyle','-.')
            xlim([-1 5.5]); xticks([-1 0 2.5 5.5]); xlabel('Time from reward onset(s)')
            ylabel('Normalized Ca trace');
        end
    end
    saveas(f,['E:\data\vHPC\all\figure\RPE\rwprob_trace_',mouseG{iHP},'.tif'])
end


%% history effect
% load('historyeffect_rwprob.mat')
C_z_dHP = cat(3, C_z_all_dHP03,C_z_all_dHP04,C_z_all_dHP06,C_z_all_dHP07,C_z_all_dHP08,C_z_all_dHP10);
C_z_vHP = cat(3, C_z_all_vHP06,C_z_all_vHP07,C_z_all_vHP08,C_z_all_vHP11,C_z_all_vHP12,C_z_all_vHP13);
rwhistory_dHP_raw = cat(3, rwhistory_dHP03,rwhistory_dHP04,rwhistory_dHP06,rwhistory_dHP07,rwhistory_dHP08,rwhistory_dHP10);
rwhistory_vHP_raw = cat(3, rwhistory_vHP06,rwhistory_vHP07,rwhistory_vHP08,rwhistory_vHP11,rwhistory_vHP12,rwhistory_vHP13);
rwhistory_dHP_mean = cellfun(@(x) permute(x,[3 2 1]), cellfun(@mean,rwhistory_dHP_raw,'UniformOutput', 0),'UniformOutput', 0);
rwhistory_vHP_mean = cellfun(@(x) permute(x,[3 2 1]), cellfun(@mean,rwhistory_vHP_raw,'UniformOutput', 0),'UniformOutput', 0);
% r-->r u-->r
% r-->u u-->u
for ii = 1:2
    for jj = 1:2
        rwhistory_dHP{ii,jj} = cat(1, rwhistory_dHP_mean{ii,jj,1},rwhistory_dHP_mean{ii,jj,2},rwhistory_dHP_mean{ii,jj,3},...
            rwhistory_dHP_mean{ii,jj,4},rwhistory_dHP_mean{ii,jj,5},rwhistory_dHP_mean{ii,jj,6});
        rwhistory_vHP{ii,jj} = cat(1, rwhistory_vHP_mean{ii,jj,1},rwhistory_vHP_mean{ii,jj,2},rwhistory_vHP_mean{ii,jj,3},...
            rwhistory_vHP_mean{ii,jj,4},rwhistory_vHP_mean{ii,jj,5},rwhistory_vHP_mean{ii,jj,6});
    end
end
trial_dHP = cat(3, trial_all_dHP03,trial_all_dHP04,trial_all_dHP06,trial_all_dHP07,trial_all_dHP08,trial_all_dHP10);
trial_vHP = cat(3, trial_all_vHP06,trial_all_vHP07,trial_all_vHP08,trial_all_vHP11,trial_all_vHP12,trial_all_vHP13);
waterReward_dHP = cat(2, waterReward_dHP03,waterReward_dHP04,waterReward_dHP06,waterReward_dHP07,waterReward_dHP08,waterReward_dHP10);
waterReward_vHP = cat(2, waterReward_vHP06,waterReward_vHP07,waterReward_vHP08,waterReward_vHP11,waterReward_vHP12,waterReward_vHP13);
cmap1= [30 45 232; 176 196 222; 128 0 128; 221 160 221]./255;
%%
f=figure('PaperUnits','Centimeters','PaperPosition',[2 2 25 8]);
fr=30;
for iHP = 1:2
    if iHP ==1; rwh = rwhistory_dHP; elseif iHP==2; rwh = rwhistory_vHP; end
    subplot(1,2,iHP)
    hold on
    for it0 = 1:2
        for ipt = 1:2
            stdshade(rwh{it0,ipt},0.3, cmap1(2*(it0-1)+ipt,:))
        end
    end
    lim = axis;
    line([0.5*fr 0.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([2*fr 2*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([3.5*fr 3.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([4*fr 4*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([5.5*fr 5.5*fr], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([8*fr 8*fr], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    xlim([0 11*fr]); xticks(fr*[1.2 2.8 4.7 6.7]); xticklabels({'C', 'D1','D2','O'})
    
end

saveas(f,['E:\data\vHPC\all\figure\zscore\reward history.tif'])
%% value sig neuron 
load('period_reg_rwprob_20220228.mat')
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
sign_dHP = pvalueset_dHP(5,7,:)<0.05;sign_vHP = pvalueset_vHP(5,7,:)<0.05;
load('cue_order.mat')
cueorder_dHP = cat(1,cueorder_dHP03,cueorder_dHP04,cueorder_dHP06,cueorder_dHP07,cueorder_dHP08,cueorder_dHP10);
cueorder_vHP = cat(1,cueorder_vHP06,cueorder_vHP07,cueorder_vHP08,cueorder_vHP11,cueorder_vHP12,cueorder_vHP14);
value_dHP = cueorder_dHP(sign_dHP,:); value_vHP = cueorder_vHP(sign_vHP,:); 


mean(value_dHP(:,3))
mean(value_vHP(:,3))


