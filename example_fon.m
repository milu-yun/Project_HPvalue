
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
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP10,pvalueset_vHP11);
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
    
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 25 6]);
    subplot(1,3,1);
    hold on
    p1=plot(FON(1,:), 'color', [129 212 250]./255, 'linewidth', 1);%cue1
    p2=plot(FON(2,:), 'color', [239 154 154]./255, 'linewidth', 1);%cue2
    p3=plot(FON(3,:), 'color', [0 0 205]./255,'linewidth', 1);%Rw
    p4=plot(FON(4,:), 'color', cmap1(4,:), 'linewidth', 1);%Pn
    p5=plot(FON(5,:), 'color', [80 80 80]./255,      'linewidth', 1);%Lick
    ylim([0 0.7])
    xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
    title('t','FontSize',14)
    %legend([p1,p2,p3,p4,p5],'Rw Cue', 'Pn Cue','Rw','Pn','Lick');
    
    subplot(1,3,2)
    hold on
    p6=plot(FON(6,:), 'color', [129 212 250]./255,'linewidth', 1);%cue1(Rwcue t-1)
    p7=plot(FON(7,:), 'color', [239 154 154]./255,'linewidth', 1);%cue2(Pncue t-1)
    p8=plot(FON(8,:), 'color', [0 0 205]./255,'linewidth', 1);%Rw(t-1)
    p9=plot(FON(9,:), 'color', cmap1(4,:),'linewidth', 1);%Pn(t-1)
    ylim([0 0.7])
    xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    lim = axis; ax = gca; ax.TickLength = [0 0];
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
    title('t-1','FontSize',14)
    
    
    subplot(1,3,3)
    hold on
    p10=plot(FON(10,:), 'color', [129 212 250]./255,'linewidth', 1);%cue1(Rwcue t-1)
    p11=plot(FON(11,:), 'color', [239 154 154]./255,'linewidth', 1);%cue2(Pncue t-1)
    p12=plot(FON(12,:), 'color', [0 0 205]./255,'linewidth', 1);%Rw(t-1)
    p13=plot(FON(13,:), 'color', cmap1(4,:),'linewidth', 1);%Pn(t-1)
    ylim([0 0.7])
    xlim([1 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', '-.'); %2.5 outcome
    line([lim(1) lim(2)], [bin_percent bin_percent], 'color',[145 84 71]./255, 'linestyle', ':','linewidth', 0.75)
    title('t-2','FontSize',14)
    
    saveas(f1,[cd,'\figure\','FON_C15',mouseG{ii},d,'.tif'])
end
%% FON shuffle
clear
shuffledir = 'E:\data\vHPC\all\shuffle\seltrial\';
dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
    'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
for imouse = 1:12
    load([shuffledir,dirn{imouse},'_rwpn_fon_shuffle_pnseltrial_01_1.mat'])
    eval(['srcset_r_',dirn{imouse},' = srcset_r;'])
    eval(['pvalueset_r_',dirn{imouse},' = pvalueset_r;'])    
end
save('fon_rwpn_shuffle_01_selpn.mat','srcset_r_dHP03','srcset_r_dHP04','srcset_r_dHP06','srcset_r_dHP07','srcset_r_dHP08','srcset_r_dHP10',...
    'srcset_r_vHP06','srcset_r_vHP07','srcset_r_vHP08','srcset_r_vHP11','srcset_r_vHP12','srcset_r_vHP14',...
    'pvalueset_r_dHP03','pvalueset_r_dHP04','pvalueset_r_dHP06','pvalueset_r_dHP07','pvalueset_r_dHP08','pvalueset_r_dHP10',...
    'pvalueset_r_vHP06','pvalueset_r_vHP07','pvalueset_r_vHP08','pvalueset_r_vHP11','pvalueset_r_vHP12','pvalueset_r_vHP14')
% clear
% shuffledir = 'E:\data\vHPC\all\shuffle\outcome_01\';
% dirn = {'dHP03','dHP04','dHP06','dHP07','dHP08','dHP10',...
%     'vHP06','vHP07','vHP08','vHP11','vHP12','vHP14'};
% for imouse = 1:12
%     load([shuffledir,dirn{imouse},'_rwpn_fon_shuffle_outcome_1.mat'])
%     eval(['srcset_r_',dirn{imouse},' = srcset_r;'])
%     eval(['pvalueset_r_',dirn{imouse},' = pvalueset_r;'])    
% end
% save('fon_rwpn_shuffle_outcome01.mat','srcset_r_dHP03','srcset_r_dHP04','srcset_r_dHP06','srcset_r_dHP07','srcset_r_dHP08','srcset_r_dHP10',...
%     'srcset_r_vHP06','srcset_r_vHP07','srcset_r_vHP08','srcset_r_vHP11','srcset_r_vHP12','srcset_r_vHP14',...
%     'pvalueset_r_dHP03','pvalueset_r_dHP04','pvalueset_r_dHP06','pvalueset_r_dHP07','pvalueset_r_dHP08','pvalueset_r_dHP10',...
%     'pvalueset_r_vHP06','pvalueset_r_vHP07','pvalueset_r_vHP08','pvalueset_r_vHP11','pvalueset_r_vHP12','pvalueset_r_vHP14')

%% Draw value FON all case
for icase = [6]
    clearvars -except icase
    if icase==1; %cue(t,t-1,t-2)
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwpn_shuffle01.mat';
%         loadshufflename2= 'fon_rwpn_shuffle_outcome01.mat';
        loadshufflename2= 'fon_rwpn_shuffle01.mat';
        loaddataname = 'fon_rwpn_01_20220302.mat';
        %     figureind = 1;
        figureind = [1 2 3 4 7 8 9 10 15 16 17 18];        
    elseif icase==2; %Rw
        deltafon=false; lickwithout = true;
        loadshufflename1= 'fon_rwpn_shuffle_01_selrw.mat';
        loadshufflename2= 'fon_rwpn_shuffle_01_selrw.mat';
        loaddataname = 'fon_rwpn_01_selrw_20220319.mat';
        figureind = 5;
    elseif icase==3; %Pn
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwpn_shuffle_01_selpn.mat';
        loadshufflename2= 'fon_rwpn_shuffle_01_selpn.mat';
        loaddataname = 'fon_rwpn_01_selpn_20220319.mat';
        figureind = 6;
    elseif icase==4; %Rw (t-1,t-2)
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwpn_shuffle01.mat';
%         loadshufflename2= 'fon_rwpn_shuffle_outcome01.mat';
        loadshufflename2= 'fon_rwpn_shuffle01.mat';
        loaddataname = 'fon_rwpn_01_20220302.mat';
        figureind = [11 12 19 20];
    elseif icase==5; %Pn (t-1,t-2)
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwpn_shuffle01.mat';
%         loadshufflename2= 'fon_rwpn_shuffle_outcome01.mat';
        loadshufflename2= 'fon_rwpn_shuffle01.mat';
        loaddataname = 'fon_rwpn_01_20220302.mat';
        figureind = [13 14 21 22];
    elseif icase==6; %lick speed
        deltafon=false; lickwithout = false;
        loadshufflename1= 'fon_rwpn_shuffle01.mat';
        loadshufflename2= 'fon_rwpn_shuffle01.mat';
        loaddataname = 'fon_rwpn_01_20220302.mat';
        figureind = [23 24 25];
    end
    if deltafon; deltaname='_delta'; else; deltaname=''; end
    if lickwithout; lickname='_nolick'; else; lickname='';end
    
    
    cmap = [55 0 102; 179 143 177]./255;
    load(loaddataname)
    pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
    pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
    pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;% remove lick p
    pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;% remove lick p
    if size(pvaluesig_dHP,1)<14
        pvaluesig_dHP = cat(1,zeros(2,size(pvaluesig_dHP,2),size(pvaluesig_dHP,3)),pvaluesig_dHP);
        pvaluesig_vHP = cat(1,zeros(2,size(pvaluesig_vHP,2),size(pvaluesig_vHP,3)),pvaluesig_vHP);
    end
    
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
    cue = [1 2 5 6 7 8 11 12]; outcome = [3 4 9 10 13 14];
    % cue = 1:14;
    load(loadshufflename1)    
    pvalueset_r_dHP_tmp = cat(3,pvalueset_r_dHP03,pvalueset_r_dHP04,pvalueset_r_dHP06,pvalueset_r_dHP07,pvalueset_r_dHP08,pvalueset_r_dHP10);
    pvalueset_r_vHP_tmp = cat(3,pvalueset_r_vHP06,pvalueset_r_vHP07,pvalueset_r_vHP08,pvalueset_r_vHP11,pvalueset_r_vHP12,pvalueset_r_vHP14);
    if size(pvalueset_r_dHP_tmp,1)<15
        pvalueset_r_dHP_tmp = cat(1,ones(2,size(pvalueset_r_dHP_tmp,2),size(pvalueset_r_dHP_tmp,3)),pvalueset_r_dHP_tmp);
        pvalueset_r_vHP_tmp = cat(1,ones(2,size(pvalueset_r_vHP_tmp,2),size(pvalueset_r_vHP_tmp,3)),pvalueset_r_vHP_tmp);
    end
    pvalueset_r_dHP(cue,:,:) = pvalueset_r_dHP_tmp(cue+1,:,:);
    pvalueset_r_vHP(cue,:,:) = pvalueset_r_vHP_tmp(cue+1,:,:);
    load(loadshufflename2)
    pvalueset_r_dHP_tmp  = cat(3,pvalueset_r_dHP03,pvalueset_r_dHP04,pvalueset_r_dHP06,pvalueset_r_dHP07,pvalueset_r_dHP08,pvalueset_r_dHP10);
    pvalueset_r_vHP_tmp  = cat(3,pvalueset_r_vHP06,pvalueset_r_vHP07,pvalueset_r_vHP08,pvalueset_r_vHP11,pvalueset_r_vHP12,pvalueset_r_vHP14);
    if size(pvalueset_r_dHP_tmp,1)<15
        pvalueset_r_dHP_tmp = cat(1,ones(2,size(pvalueset_r_dHP_tmp,2),size(pvalueset_r_dHP_tmp,3)),pvalueset_r_dHP_tmp);
        pvalueset_r_vHP_tmp = cat(1,ones(2,size(pvalueset_r_vHP_tmp,2),size(pvalueset_r_vHP_tmp,3)),pvalueset_r_vHP_tmp);
    end
    pvalueset_r_dHP(outcome,:,:) = pvalueset_r_dHP_tmp(outcome+1,:,:);
    pvalueset_r_vHP(outcome,:,:) = pvalueset_r_vHP_tmp(outcome+1,:,:);
    pvaluesig_r_dHP = pvalueset_r_dHP(:,:,:)<0.05;
    pvaluesig_r_vHP = pvalueset_r_vHP(:,:,:)<0.05;
    if lickwithout
        load('E:\data\vHPC\all\rwpn\lick_neuron.mat');
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
        chis_dHP =sum(pvaluesig_dHP, 3); chis_vHP = sum(pvaluesig_vHP, 3);
    end

    % figurename = {'RC(t)','PC(t)','R(t)','P(t)','Lick','Speed','RC(t-1)','PC(t-1)','R(t-1)','P(t-1)','RC(t-2)','PC(t-2)','R(t-2)','P(t-2)'};
    psig= [1 1 2 2 3 4 7 7 8 8 9 9 10 10, 11 11 12 12 13 13 14 14 5 6 6]; % for figure
    xlimlist = [1 50; 51 80; 1 50; 51 80; 51 80; 51 80;...% (t)
        1 50; 51 80; 1 50; 51 80; 1 50; 51 80; 1 50; 51 80;...% (t-1)
        1 50; 51 80; 1 50; 51 80; 1 50; 51 80; 1 50; 51 80;...
        35 85'; 1 50; 51 80;]; % (t-2)
    wlist = (diff(xlimlist')+1)*0.075; %0.075  0.055
    
    termlist = {'Rw cue1','Rw cue2','Pn cue1','Pn cue2','Rw','Pn',...
        'Rw cue(t-1)1','Rw cue(t-1)2','Pn cue(t-1)1','Pn cue(t-1)2','Rw(t-1)1','Rw(t-1)2','Pn(t-1)1','Pn(t-1)2',...
        'Rw cue(t-2)1','Rw cue(t-2)2','Pn cue(t-2)1','Pn cue(t-2)2','Rw(t-2)1','Rw(t-2)2','Pn(t-2)1','Pn(t-2)2',...
        'lick(t)','speed(t)1','speed(t)2'};    
    
    for isig = figureind
        if ismember(psig(isig),[1,7,11]); ymax = 0.6;
        elseif ismember(psig(isig),[2,8,12]); ymax = 0.2;
        elseif ismember(psig(isig),[3,9,13]); ymax = 0.3;
        elseif ismember(psig(isig),[10,14]); ymax = 0.2;
        elseif psig(isig)==4; ymax=0.3;
        else ymax = 0.5; end
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
        %xlim([50 120]); xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
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
%         saveas(f1,['E:\data\vHPC\all\figure\FON\rwpn\',termlist{isig},deltaname,lickname,'.tif'])
        print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\',termlist{isig},deltaname,lickname,'.ai']);
    end
end
%% FON period - cue
cmap = [55 0 102; 179 143 177]./255;
load('fon_period_220320.mat')
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;% remove lick p
pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;% remove lick p

for ii=1:2
    f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1 1]);%2.5 2.5
    if ii==1; ymax = 0.2; iterm = 2; tname = 'Pncue';
    elseif ii==2; ymax = 0.4; iterm = 1; tname = 'Rwcue';
    end
    hold on
    bar(0.25,mean(pvaluesig_dHP(iterm,7,:)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
%     bar(0.5,mean(pvaluesig_dHP(iterm,8,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
    
    bar(0.5,mean(pvaluesig_vHP(iterm,7,:)),0.2,'facecolor',cmap(2,:),'edgecolor','none')
%     bar(1.1,mean(pvaluesig_vHP(iterm,8,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
    
    xlim([0 0.75]); xticks([]); xticklabels({});
    ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
%     saveas(f2,['E:\data\vHPC\all\figure\FON\rwpn\fig4_fon_period',tname,'.tif'])
    print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\fon_period_',tname,'.ai']);
    [p,stat]=chis([sum(pvaluesig_dHP(iterm,7,:)) size(pvaluesig_dHP(iterm,7,:),3)-sum(pvaluesig_dHP(iterm,7,:));...
        sum(pvaluesig_vHP(iterm,7,:)) size(pvaluesig_vHP(iterm,7,:),3)-sum(pvaluesig_vHP(iterm,7,:))])
    [p,stat]=chis([sum(pvaluesig_dHP(iterm,7,:)) size(pvaluesig_dHP(iterm,7,:),3)-sum(pvaluesig_dHP(iterm,7,:));...
        sum(pvaluesig_dHP(iterm,8,:)) size(pvaluesig_dHP(iterm,8,:),3)-sum(pvaluesig_dHP(iterm,8,:))])
    [p,stat]= chis([sum(pvaluesig_vHP(iterm,7,:)) size(pvaluesig_vHP(iterm,7,:),3)-sum(pvaluesig_vHP(iterm,7,:));...
        sum(pvaluesig_vHP(iterm,8,:)) size(pvaluesig_vHP(iterm,8,:),3)-sum(pvaluesig_vHP(iterm,8,:))])
end
%% FON period - value outcome period
cmap = [55 0 102; 179 143 177]./255;
load('fon_period_220320.mat')
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;% remove lick p
pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;% remove lick p

for ii=1:2
    f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1 1]);%2.5 2.5
    if ii==1; ymax = 0.2; iterm = 2; tname = 'Pncue';
    elseif ii==2; ymax = 0.6; iterm = 1; tname = 'Rwcue';
    end
    hold on
    bar(0.25,mean(pvaluesig_dHP(iterm,1,:)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
%     bar(0.5,mean(pvaluesig_dHP(iterm,9,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
    
    bar(0.5,mean(pvaluesig_vHP(iterm,1,:)),0.2,'facecolor',cmap(2,:),'edgecolor','none')
%     bar(1.1,mean(pvaluesig_vHP(iterm,9,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
    
    xlim([0 0.75]); xticks([]); xticklabels({});
    ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
%     saveas(f2,['E:\data\vHPC\all\figure\FON\rwpn\fig4_fon_period',tname,'.tif'])
    print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\fon_outcome_period_',tname,'.ai']);
    [p,stat]=chis([sum(pvaluesig_dHP(iterm,1,:)) size(pvaluesig_dHP(iterm,1,:),3)-sum(pvaluesig_dHP(iterm,1,:));...
        sum(pvaluesig_vHP(iterm,1,:)) size(pvaluesig_vHP(iterm,1,:),3)-sum(pvaluesig_vHP(iterm,1,:))])
    [p,stat]=chis([sum(pvaluesig_dHP(iterm,1,:)) size(pvaluesig_dHP(iterm,1,:),3)-sum(pvaluesig_dHP(iterm,1,:));...
        sum(pvaluesig_dHP(iterm,9,:)) size(pvaluesig_dHP(iterm,9,:),3)-sum(pvaluesig_dHP(iterm,9,:))])
    [p,stat]= chis([sum(pvaluesig_vHP(iterm,1,:)) size(pvaluesig_vHP(iterm,1,:),3)-sum(pvaluesig_vHP(iterm,1,:));...
        sum(pvaluesig_vHP(iterm,9,:)) size(pvaluesig_vHP(iterm,9,:),3)-sum(pvaluesig_vHP(iterm,9,:))])
end


%% FON period - outcome
cmap = [55 0 102; 179 143 177]./255;
% load('fon_period.mat')
load('fon_period_220320.mat')
pvalueset_dHP = cat(3,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(3,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
pvaluesig_dHP = pvalueset_dHP(2:end,:,:)<0.05;% remove lick p
pvaluesig_vHP = pvalueset_vHP(2:end,:,:)<0.05;% remove lick p
%10 11 12 13 14 15 for rw
%16 17 18 19 20 21 for pn
nlcheck=true;
for ii=1:2
    for iperiod = 1
        f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1 2.5]);%
        if ii==1; ymax = 0.2; periodlist = [16 17 18 19 20 21];iterm = 4; tname = 'Pn';
        elseif ii==2; ymax = 0.2; periodlist = [10 11 12 13 14 15];iterm = 3; tname = 'Rw';
            if nlcheck==true;
                load('E:\data\vHPC\all\rwpn\lick_neuron.mat');
                lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
                lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
                pvaluesig_dHP = pvaluesig_dHP(:,:,lick_dHP(:,1)>=0.05);
                pvaluesig_vHP = pvaluesig_vHP(:,:,lick_vHP(:,1)>=0.05);
                nlcheck = false;
            end
        end
        hold on
        bar(0.25,mean(pvaluesig_dHP(iterm,periodlist(2*iperiod-1),:)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
        %         bar(0.5,mean(pvaluesig_dHP(iterm,2*iperiod,:))- mean(pvaluesig_dHP(iterm,2*iperiod,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
        
        bar(0.5,mean(pvaluesig_vHP(iterm,periodlist(2*iperiod-1),:)),0.2,'facecolor',cmap(2,:),'edgecolor','none')
        %         bar(1.1,mean(pvaluesig_vHP(iterm,2*iperiod,:))-mean(pvaluesig_vHP(iterm,2*iperiod,:)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
        
        xlim([0 0.75]); xticks([]); xticklabels({});
        ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
        %         saveas(f2,['E:\data\vHPC\all\figure\FON\rwpn\fig4_fon_period',tname,num2str(iperiod),'.tif'])
        print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\fon_period',tname,num2str(iperiod),'.ai']);
        dtmp =sum(pvaluesig_dHP(iterm,periodlist(2*iperiod-1),:));
        vtmp =sum(pvaluesig_vHP(iterm,periodlist(2*iperiod-1),:));
        [p,stat] = chis([dtmp length(pvaluesig_dHP(iterm,periodlist(2*iperiod-1),:))-dtmp;...
             vtmp length(pvaluesig_vHP(iterm,periodlist(2*iperiod-1),:))-vtmp])
    end
end
%% SRC pereiod
clear
cmap = [55 0 102; 179 143 177]./255;
% load('fon_period.mat')
load('fon_period01.mat')
srcset_dHP = cat(3,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
srcset_vHP = cat(3,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);

abs_dHP = permute(abs(srcset_dHP),[3 2 1]);
abs_vHP = permute(abs(srcset_vHP),[3 2 1]);
nlcheck=true;
for ii=1:2
    for iperiod = 1:3
        f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 1]);%
        if ii==1; ymax = 0.1; iterm = 4; tname = 'Pn';
        elseif ii==2; ymax = 0.2; iterm = 3; tname = 'Rw';
            if nlcheck==true;
                load('E:\data\vHPC\all\rwpn\lick_neuron.mat');
                lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
                lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
                abs_dHP = abs_dHP(lick_dHP(:,1)>=0.05,:,:);
                abs_vHP = abs_vHP(lick_vHP(:,1)>=0.05,:,:);
                nlcheck = false;
            end
        end
        hold on
        bar(0.25,mean(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm)),0.2,'facecolor',cmap(1,:),'edgecolor','none')
        errorbar(0.25,mean(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm)),sem(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm)),'color','k','Capsize',3)
        %         bar(0.5,mean(abs_dHP(:,2*iperiod,iterm)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
        %         errorbar(0.5,mean(abs_dHP(:,2*iperiod,iterm)),sem(abs_dHP(:,2*iperiod,iterm)),'color','k','Capsize',3) % shuffle
        
        bar(0.85,mean(abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm)),0.2,'facecolor',cmap(2,:),'edgecolor','none')
        errorbar(0.85,mean(abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm)),sem(abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm)),'color','k','Capsize',3)
        %         bar(1.1,mean(abs_vHP(:,2*iperiod,iterm)),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
        %         errorbar(1.1,mean(abs_vHP(:,2*iperiod,iterm)),sem(abs_vHP(:,2*iperiod,iterm)),'color','k','Capsize',3) % shuffle
        
        xlim([0 1.35]); xticks([0.6 0.9]); xticklabels({});
        ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
                saveas(f2,['E:\data\vHPC\all\figure\FON\rwpn\fig4_src_period',tname,num2str(iperiod),'.tif'])
        %         print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\src_period',tname,num2str(tmpid),'.ai']);
        
        
        %         anova_d = [abs_dHP(:,2*iperiod-1,iterm);abs_dHP(:,2*iperiod,iterm)];
        %         anova_v = [abs_vHP(:,2*iperiod-1,iterm);abs_vHP(:,2*iperiod,iterm)];
        %         mouseG  = [zeros(size(anova_d));ones(size(anova_v))];
        %         termG = [1*ones(size(abs_dHP(:,2*iperiod-1,iterm))); 2*ones(size(abs_dHP(:,2*iperiod,iterm))); ...
        %             1*ones(size(abs_vHP(:,2*iperiod-1,iterm))); 2*ones(size(abs_vHP(:,2*iperiod,iterm))); ];
        %         [~,tbl,stat] = anovan([anova_d;anova_v],{mouseG,termG},'model','interaction','varnames',{'mouse','term'},'display','off');
        %         c = multcompare(stat, 'dimension', [1 2],'display','off');
        %         cset{ii,iperiod} = c;
        [~,p]=ttest2(abs_dHP(:,2*iperiod-1,iterm)-abs_dHP(:,2*iperiod,iterm),abs_vHP(:,2*iperiod-1,iterm)-abs_vHP(:,2*iperiod,iterm));
        pset(ii,iperiod) = p;
    end
end

%% Draw value src all case
clear
cmap = [55 0 102; 179 143 177]./255;
load('fon_rwpn_01_20220302.mat')
% load('fon_rwpn_-101_20220302.mat')
srcset_dHP = cat(3,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
srcset_vHP = cat(3,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);

abs_dHP = permute(abs(srcset_dHP),[3 2 1]);
abs_vHP = permute(abs(srcset_vHP),[3 2 1]);

load('fon_rwpn_shuffle_outcome01.mat')
srcset_r_dHP = cat(3,srcset_r_dHP03,srcset_r_dHP04,srcset_r_dHP06,srcset_r_dHP07,srcset_r_dHP08,srcset_r_dHP10);
srcset_r_vHP = cat(3,srcset_r_vHP06,srcset_r_vHP07,srcset_r_vHP08,srcset_r_vHP11,srcset_r_vHP12,srcset_r_vHP14);
abs_r_dHP = permute(abs(srcset_r_dHP),[3 2 1]);
abs_r_vHP = permute(abs(srcset_r_vHP),[3 2 1]);
abs_dHP = abs_dHP-abs_r_dHP;
abs_vHP = abs_vHP-abs_r_vHP;

load('E:\data\vHPC\all\rwpn\lick_neuron.mat');
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
abs_dHP = abs_dHP(lick_dHP(:,1)>=0.05,:,:);
abs_vHP = abs_vHP(lick_vHP(:,1)>=0.05,:,:);

figurename = {'RC(t)','PC(t)','R(t)','P(t)','Lick','Speed','RC(t-1)','PC(t-1)','R(t-1)','P(t-1)','RC(t-2)','PC(t-2)','R(t-2)','P(t-2)'};


for ii =3%:14
    for tt = 1:120  
        sig_t(1,tt) = ttest2(abs_dHP(:,tt,ii),abs_vHP(:,tt,ii));
    end
    if ismember(ii,[1,7,11]); ymax = 0.3;
    elseif ismember(ii,[2,8,12]); ymax = 0.15;
    elseif ismember(ii,[4,10,14]); ymax = 0.1;
    elseif ii==6; ymax = 0.15;
    else ymax = 0.2; end
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
    hold on
    stdshade(abs_dHP(:,:,ii),0.3,cmap(1,:))
    stdshade(abs_vHP(:,:,ii),0.3,cmap(2,:))
%     stdshade(abs_r_dHP(:,:,ii),0.3,cmap(1,:),[],':')
%     stdshade(abs_r_vHP(:,:,ii),0.3,cmap(2,:),[],':')
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
%     saveas(f1,['E:\data\vHPC\all\figure\FON\rwpn\fig4_src_raw',figurename{ii},'.tif'])
    % print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\',figurename{ii},'.ai']);
end



%% fon value in delay1
load('period_reg_rwpn.mat')
pvalueset_dHP = cat(2,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
pvalueset_vHP = cat(2,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
p_d = sum(pvalueset_dHP([2,3],:,1)<0.05,2); n_d = size(pvalueset_dHP,2);
p_v = sum(pvalueset_vHP([2,3],:,1)<0.05,2); n_v = size(pvalueset_vHP,2);

chis([p_d(1) n_d-p_d(1); p_v(1) n_v-p_v(1);])
chis([p_d(2) n_d-p_d(2); p_v(2) n_v-p_v(2);])
f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 2.5]);
hold on
bar(1.1,p_d(1)./n_d,'facecolor',[148 69 40]./255,'barwidth',0.8);
bar(2.1,p_v(1)./n_v,'facecolor',[25 148 123]./255,'barwidth',0.8);
bar(3.4,p_d(2)./n_d,'facecolor',[148 69 40]./255,'barwidth',0.8);
bar(4.4,p_v(2)./n_v,'facecolor',[25 148 123]./255,'barwidth',0.8);
xticks([]); xlim([0.4 5.1]); yticks([0 .5 1]); ylim([0 1])
% print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\FON\rwpn\barg.ai']);
%% Draw CPD
% load('cpd_rwpn'); load('cpd_rwpn_withv');
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);

mouseG = {'dHP','vHP'};
for ii = 1:2
    eval(['sse_set = sse_set_',mouseG{ii},';']);
    day = datetime('today');
    d = datestr(day,'yymmdd');
    cpd = (sse_set(2:end,:,:)-sse_set(1,:,:))./sse_set(2:end,:,:);
    cpd = permute(cpd,[3 2 1]);
    yl = [0 0.08];
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 25 6]);
    hxa(1) = subplot(1,3,1);
    hold on
%     stdshade(cpd(:,:,14), 0.3, [160 82 45]./255)%speed
    stdshade(cpd(:,:,5), 0.3, [80 80 80]./255)%Lick
    stdshade(cpd(:,:,4), 0.3, [255 30 70]./255)%pn
    stdshade(cpd(:,:,3), 0.3, [0 0 205]./255)%rw
    stdshade(cpd(:,:,1), 0.3, [129 212 250]./255)%rwcue
    stdshade(cpd(:,:,2), 0.3, [239 154 154]./255)%pncue
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
    stdshade(cpd(:,:,6), 0.3, [129 212 250]./255)%rwcue
    stdshade(cpd(:,:,7), 0.3, [239 154 154]./255)%pncue
    stdshade(cpd(:,:,8), 0.3, [0 0 205]./255)%rw
    stdshade(cpd(:,:,9), 0.3, [255 30 70]./255)%pn
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
    stdshade(cpd(:,:,10), 0.3, [129 212 250]./255)%rwcue
    stdshade(cpd(:,:,11), 0.3, [239 154 154]./255)%pncue
    stdshade(cpd(:,:,12), 0.3, [0 0 205]./255)%rw
    stdshade(cpd(:,:,13), 0.3, [255 30 70]./255)%pn
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %-0.5 cue
    line([35 35], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [yl(1) yl(2)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [yl(1) yl(2)],'color','k', 'linestyle', '-.'); %2.5 outcome
    title('t-2','FontSize',14)
    hold off
    xticks(hxa,[12 27 47 67 100]); xticklabels(hxa, {'C', 'D1','D2','O','ITI'})
    ylim(hxa, yl); xlim(hxa, [1 120]);
%     saveas(f1,['E:\data\vHPC\all\figure\cpd\rwpn\','CPD_rwpn_',mouseG{ii},d,'.tif'])
end
%% CPD anova
%Rwcue, Pncue, Rw, Pn, Lick, (t-1) (t-2)
load('cpd_rwpn_20211118');lickwithout = false;subt='_';shufName = 'cpd_rwpn_shuffle_ttest_Outcome.mat';
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
load('lick_neuron.mat'); lickwithout = true;subt='_nolick_';
lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);

%for speed effect - dHP
% load('cpd_rwpn_20211118');lickwithout = false;subt='_dHP_comparespeed_';
% sse_set_dHP = cat(3,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
% sse_set_dHP  = sse_set_dHP (1:14,:,:);
% load('CPD_no_speed');
% sse_set_vHP = cat(3,sse_set__nospeed_dHP07,sse_set__nospeed_dHP08,sse_set__nospeed_dHP10);
% sse_set_vHP  = sse_set_vHP (1:14,:,:);
%for speed effect - vHP
% load('cpd_rwpn_20211118');lickwithout = false;subt='_vHP_comparespeed_';
% sse_set_dHP = cat(3,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
% load('CPD_no_speed');
% sse_set_vHP = cat(3,sse_set__nospeed_vHP07,sse_set__nospeed_vHP08,sse_set__nospeed_vHP11,sse_set__nospeed_vHP12,sse_set__nospeed_vHP14);


mouseG = {'dHP','vHP'};
day = datetime('today');
d = datestr(day,'yymmdd');
for ii = 1:2
    eval(['sse_set = sse_set_',mouseG{ii},';']);
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

% for iterm =1:13
% for itime= 1:size(cpd_dHP,2);
% [~,p] = ttest2(cpd_dHP(:,itime,iterm),cpd_vHP(:,itime,iterm));
% panovaset(itime,iterm) = p;
% end
% end
%% cpd delta
day = datetime('today');
d = datestr(day,'yymmdd');
%Rwcue, Pncue, Rw, Pn, Lick, (t-1) (t-2)
load('cpd_rwpn01_20220303');lickwithout = false; subt='_'; shufName = 'cpd_rwpn_shuffle_ttest_Outcome01.mat';
% load('cpd_with rw-101.mat');lickwithout = false; subt='_'; shufName = 'cpd_rwprob_shuffle_ttest.mat';
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
% load('lick_neuron.mat'); lickwithout = true;subt='_nolick';
% lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
% lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
% sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
% sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);
cpd_dHP = (sse_set_dHP(2:end,:,:)-sse_set_dHP(1,:,:))./sse_set_dHP(2:end,:,:); cpd_dHP = permute(cpd_dHP, [3 2 1]);
cpd_vHP = (sse_set_vHP(2:end,:,:)-sse_set_vHP(1,:,:))./sse_set_vHP(2:end,:,:); cpd_vHP = permute(cpd_vHP, [3 2 1]);
mouseG = {'dHP','vHP'};
load(shufName)
% cpd_dHP_r = cpd_dHP_r(lick_dHP(:,1)>=0.05,:,:);
% cpd_vHP_r = cpd_vHP_r(lick_vHP(:,1)>=0.05,:,:);

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
    load([shuffledir,dirn{imouse},'_rwpn_cpd_shuffle_011.mat'])
    eval(['sse_set_r_',dirn{imouse},' = sse_set_r;'])
    cpd_r =(sse_set_r(2:end,:,:,:)-sse_set_r(1,:,:,:))./sse_set_r(2:end,:,:,:);
    eval(['cpd_r_',dirn{imouse},'=cpd_r;'])
    
end
save('cpd_rwpn_shuffle01.mat','cpd_r_dHP03','cpd_r_dHP04','cpd_r_dHP06','cpd_r_dHP07','cpd_r_dHP08','cpd_r_dHP10',...
    'cpd_r_vHP06','cpd_r_vHP07','cpd_r_vHP08','cpd_r_vHP11','cpd_r_vHP12','cpd_r_vHP14')

% %permutation
% load('cpd_rwpn_shuffle.mat','cpd_r_dHP03','cpd_r_dHP04','cpd_r_dHP06','cpd_r_dHP07','cpd_r_dHP08','cpd_r_dHP10',...
%     'cpd_r_vHP06','cpd_r_vHP07','cpd_r_vHP08','cpd_r_vHP11','cpd_r_vHP12','cpd_r_vHP14');
% cpd_dHP_r = cat(3,cpd_r_dHP03,cpd_r_dHP04,cpd_r_dHP06,cpd_r_dHP07,cpd_r_dHP08,cpd_r_dHP10);
% cpd_vHP_r = cat(3,cpd_r_vHP06,cpd_r_vHP07,cpd_r_vHP08,cpd_r_vHP11,cpd_r_vHP12,cpd_r_vHP14);
% meancpd_d = mean(cpd_dHP,1); meancpdr_d = mean(cpd_dHP_r,3); % neuron mean
% meancpd_v = mean(cpd_vHP,1); meancpdr_v = mean(cpd_vHP_r,3); % neuron mean
% chancelevel = size(cpd_dHP_r,4)*0.975;
% sig_d = sum(permute(meancpd_d,[3 2 1])>meancpdr_d,4)>chancelevel;
% sig_v = sum(permute(meancpd_v,[3 2 1])>meancpdr_v,4)>chancelevel;
% sig_d = sig_d'; sig_v = sig_v';
% sig_d  = sig_d(:,[1:5,7:14,6]);sig_v = sig_v(:,[1:5,7:14,6]);
% save('cpd_rwpn_shuffle_permutation_sort.mat','sig_d','sig_v')

%% ttset
cpd_dHP_r = cat(3,cpd_r_dHP03,cpd_r_dHP04,cpd_r_dHP06,cpd_r_dHP07,cpd_r_dHP08,cpd_r_dHP10);
cpd_vHP_r = cat(3,cpd_r_vHP06,cpd_r_vHP07,cpd_r_vHP08,cpd_r_vHP11,cpd_r_vHP12,cpd_r_vHP14);
cpd_dHP_r = permute(cpd_dHP_r(:,:,:),[3 2 1]);
cpd_vHP_r = permute(cpd_vHP_r(:,:,:),[3,2,1]);
cpd_dHP_r = cpd_dHP_r(:,:,[1:5,7:14,6]);
cpd_vHP_r = cpd_vHP_r(:,:,[1:5,7:14,6]);
save('cpd_rwpn_shuffle_ttest_01.mat','cpd_dHP_r','cpd_vHP_r')
%% cpd anova plot
load(shufName)
%remove speed in shuffle  - sorting

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
cmap = [55 0 102; 179 143 177]./255;
%psig= 1:13; xlimlist = repmat([1,120],13,1); wlist = repmat([8,6],13,1);
% termlist = {'Rw cue','Pn cue','Rw','Pn','Lick','Rw cue(t-1)','Pn cue(t-1)','Rw(t-1)','Pn(t-1)',...
%     'Rw cue(t-2)','Pn cue(t-2)','Rw(t-2)','Pn(t-2)','Speed'};
psig= [1 1 2 2 3 4, 6 6 7 7 8 8 9 9, 10 10 11 11 12 12 13 13]; % for figure
xlimlist = [1 50; 51 80; 1 50; 51 80; 51 80; 51 80;...% (t)
    1 50; 51 80; 1 50; 51 80; 1 50; 51 80; 1 50; 51 80;...% (t-1)
    1 50; 51 80; 1 50; 51 80; 1 50; 51 80; 1 50; 51 80;]; % (t-2) 
wlist = (diff(xlimlist')+1)*0.04; 

termlist = {'Rw cue1','Rw cue2','Pn cue1','Pn cue2','Rw','Pn',...
    'Rw cue(t-1)1','Rw cue(t-1)2','Pn cue(t-1)1','Pn cue(t-1)2','Rw(t-1)1','Rw(t-1)2','Pn(t-1)1','Pn(t-1)2',...
    'Rw cue(t-2)1','Rw cue(t-2)2','Pn cue(t-2)1','Pn cue(t-2)2','Rw(t-2)1','Rw(t-2)2','Pn(t-2)1','Pn(t-2)2'};
for isig = 6%1:length(psig)
    if ismember(psig(isig),[1,6,10]); ymax = 0.06;
    elseif ismember(psig(isig),[2,7,11]); ymax = 0.02;
    elseif ismember(psig(isig),[4,9,13]); ymax = 0.01;
    elseif psig(isig)==14; ymax = 0.1;
    else ymax = 0.01; end
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);%
%     f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 wlist(isig) 2]);%
    hold on
    stdshade(cpd_dHP(:,:,psig(isig)), 0.3, cmap(1,:))%dorsal %pink
    stdshade(cpd_vHP(:,:,psig(isig)), 0.3, cmap(2,:))%ventral %green
%     stdshade(cpd_dHP_r(:,:,psig(isig)), 0.3, cmap(1,:),[],':') %dorsal shuffle
%     stdshade(cpd_vHP_r(:,:,psig(isig)), 0.3, cmap(2,:),[],':')%ventral shuffle
    %     title(termlist{psig(isig)});
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
    xticks([5 20 35 40 55 80]);xticklabels({}); % xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    xlim([xlimlist(isig,:)]); ylim([-0.001 ymax]); yticks([0 ymax]); yticklabels({});
    
    hold off
    saveas(f1,['E:\data\vHPC\all\figure\cpd\rwpn\fig4',termlist{isig},d,'.tif'])
%     print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwpn\',termlist{isig},subt,num2str(d),'.ai']);
%     close
    
end
%% cpd anova plot - period
load('cpd_rwpn_period.mat')
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
cpd_dHP = (sse_set_dHP(2:end,:,:)-sse_set_dHP(1,:,:))./sse_set_dHP(2:end,:,:);
cpd_vHP = (sse_set_vHP(2:end,:,:)-sse_set_vHP(1,:,:))./sse_set_vHP(2:end,:,:);
cmap = [55 0 102; 179 143 177]./255;
% cue
for iterm = 1:2 %rw cue; pn cue
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.8 2.5]);%
    if iterm==1; ymax = 0.05; tname = 'Rwcue';
    elseif iterm==2; ymax = 0.02; tname = 'Pncue'; end
    ii = 2;
    hold on
    bar(0.25,mean(cpd_dHP(iterm,ii,:),3),0.2,'facecolor',cmap(1,:),'edgecolor','none')
    errorbar(0.25,mean(cpd_dHP(iterm,ii,:),3),sem(cpd_dHP(iterm,ii,:)),'color','k','Capsize',3)
    bar(0.5,mean(cpd_dHP(iterm,ii+2,:),3),0.2,'facecolor',[180 180 180]./255,'edgecolor','none')
    errorbar(0.5,mean(cpd_dHP(iterm,ii+2,:),3),sem(cpd_dHP(iterm,ii+2,:)),'color','k','Capsize',3)
    bar(0.85,mean(cpd_vHP(iterm,ii,:),3),0.2,'facecolor',cmap(2,:),'edgecolor','none')
    errorbar(0.85,mean(cpd_vHP(iterm,ii,:),3),sem(cpd_vHP(iterm,ii,:)),'color','k','Capsize',3)
    bar(1.1,mean(cpd_vHP(iterm,ii+2,:),3),0.2,'facecolor',[180 180 180]./255,'edgecolor','none')
    errorbar(1.1,mean(cpd_vHP(iterm,ii+2,:),3),sem(cpd_vHP(iterm,ii+2,:)),'color','k','Capsize',3)
    xlim([0.1 1.25]); xticks([]); xticklabels({});
    ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
    %     saveas(f1,['E:\data\vHPC\all\figure\cpd\rwpn\mean_',tname,'.tif'])
    print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwpn\mean_delay_',tname,'.ai']);
end
%%
anova_d = [squeeze(cpd_dHP(iterm,ii,:));squeeze(cpd_dHP(iterm,ii+2,:))];
anova_v = [squeeze(cpd_vHP(iterm,ii,:));squeeze(cpd_vHP(iterm,ii+2,:))];
mouseG  = [zeros(size(anova_d));ones(size(anova_v))];
termG = [1*ones(size(squeeze(cpd_dHP(iterm,ii,:)))); 2*ones(size(squeeze(cpd_dHP(iterm,ii,:)))); ...
    1*ones(size(squeeze(cpd_vHP(iterm,ii,:)))); 2*ones(size(squeeze(cpd_vHP(iterm,ii,:))))];
[~,tbl,stat] = anovan([anova_d;anova_v],{mouseG,termG},'model','interaction','varnames',{'mouse','term'});
c = multcompare(stat, 'dimension', [1 2]);
%% outcome
load('E:\data\vHPC\all\rwpn\cpd_rwpn_period.mat')
sse_set_dHP = cat(3,sse_set_dHP03,sse_set_dHP04,sse_set_dHP06,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_vHP = cat(3,sse_set_vHP06,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
cpd_dHP = (sse_set_dHP(2:end,:,:)-sse_set_dHP(1,:,:))./sse_set_dHP(2:end,:,:);
cpd_vHP = (sse_set_vHP(2:end,:,:)-sse_set_vHP(1,:,:))./sse_set_vHP(2:end,:,:);
cmap = [55 0 102; 179 143 177]./255;
nlcheck = true;
for tmpid = 1:3 % all 1s 2s
    for ii = 1:2 %pn rw
        f2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 1]);%
        if ii==1; ymax = 0.01; termid = 4; tname = 'Pn';
        elseif ii==2; ymax = 0.02; termid = 3; tname = 'Rw';
            if nlcheck==true;
                load('E:\data\vHPC\all\rwpn\lick_neuron.mat');
                lick_dHP = cat(1, pset_lick_dHP03,pset_lick_dHP04,pset_lick_dHP06,pset_lick_dHP07,pset_lick_dHP08,pset_lick_dHP10);
                lick_vHP = cat(1, pset_lick_vHP06,pset_lick_vHP07,pset_lick_vHP08,pset_lick_vHP11,pset_lick_vHP12,pset_lick_vHP14);
                sse_set_dHP = sse_set_dHP(:,:,lick_dHP(:,1)>=0.05);
                sse_set_vHP = sse_set_vHP(:,:,lick_vHP(:,1)>=0.05);
                nlcheck = false;
            end
        end
        hold on
        bar(0.25,mean(cpd_dHP(termid,5+2*(tmpid-1),:)-cpd_dHP(termid,5+2*(tmpid-1)+1,:),3),0.2,'facecolor',cmap(1,:),'edgecolor','none')
        errorbar(0.25,mean(cpd_dHP(termid,5+2*(tmpid-1),:)-cpd_dHP(termid,5+2*(tmpid-1)+1,:),3),sem(cpd_dHP(termid,5+2*(tmpid-1),:)-cpd_dHP(termid,5+2*(tmpid-1)+1,:)),'color','k','Capsize',3)
%         bar(0.5,mean(cpd_dHP(termid,5+2*(tmpid-1)+1,:),3),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
%         errorbar(0.5,mean(cpd_dHP(termid,5+2*(tmpid-1)+1,:),3),sem(cpd_dHP(termid,5+2*(tmpid-1)+1,:)),'color','k','Capsize',3) % shuffle
        
        bar(0.85,mean(cpd_vHP(termid,5+2*(tmpid-1),:)-cpd_vHP(termid,5+2*(tmpid-1)+1,:),3),0.2,'facecolor',cmap(2,:),'edgecolor','none')
        errorbar(0.85,mean(cpd_vHP(termid,5+2*(tmpid-1),:)-cpd_vHP(termid,5+2*(tmpid-1)+1,:),3),sem(cpd_vHP(termid,5+2*(tmpid-1),:)-cpd_vHP(termid,5+2*(tmpid-1)+1,:)),'color','k','Capsize',3)
%         bar(1.1,mean(cpd_vHP(termid,5+2*(tmpid-1)+1,:),3),0.2,'facecolor',[179 179 179]./255,'edgecolor','none')
%         errorbar(1.1,mean(cpd_vHP(termid,5+2*(tmpid-1)+1,:),3),sem(cpd_vHP(termid,5+2*(tmpid-1)+1,:)),'color','k','Capsize',3) % shuffle
        
        xlim([0 1.35]); xticks([0.6 0.9]); xticklabels({});
        ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
            saveas(f2,['E:\data\vHPC\all\figure\cpd\rwpn\fig4_mean_',tname,num2str(tmpid),'.tif'])
%         print(f2,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwpn\mean_',tname,num2str(tmpid),'.ai']);
        
%         anova_d = [squeeze(cpd_dHP(termid,5+2*(tmpid-1),:));squeeze(cpd_dHP(termid,6+2*(tmpid-1),:))];
%         anova_v = [squeeze(cpd_vHP(termid,5+2*(tmpid-1),:));squeeze(cpd_vHP(termid,6+2*(tmpid-1),:))];
%         mouseG  = [zeros(size(anova_d));ones(size(anova_v))];
%         termG = [1*ones(size(squeeze(cpd_dHP(termid,5+2*(tmpid-1),:)))); 2*ones(size(squeeze(cpd_dHP(termid,6+2*(tmpid-1),:))));...
%             1*ones(size(squeeze(cpd_vHP(termid,5+2*(tmpid-1),:)))); 2*ones(size(squeeze(cpd_vHP(termid,6+2*(tmpid-1),:))))];
%         [~,tbl,stat] = anovan([anova_d;anova_v],{mouseG,termG},'model','interaction','varnames',{'mouse','term'},'display','off');
%         tblset{ii,tmpid} = tbl;
%         c = multcompare(stat, 'dimension', [1 2]);%,'display','off');
%         cset{ii,tmpid} = c;
[~,p] = ttest2(cpd_dHP(termid,5+2*(tmpid-1),:)-cpd_dHP(termid,5+2*(tmpid-1)+1,:),cpd_vHP(termid,5+2*(tmpid-1),:)-cpd_vHP(termid,5+2*(tmpid-1)+1,:))
pset(ii,tmpid) = p;
        
    end
end
% termid = 1; termid = 2; 


%% cpd anova plot - speed
%for speed effect - dHP
load('cpd_rwpn_20211118');lickwithout = false;
sse_set_dHP = cat(3,sse_set_dHP07,sse_set_dHP08,sse_set_dHP10);
sse_set_dHP  = sse_set_dHP (1:14,:,:);
load('CPD_no_speed');
sse_set_dHP_nospeed = cat(3,sse_set__nospeed_dHP07,sse_set__nospeed_dHP08,sse_set__nospeed_dHP10);
%for speed effect - vHP
load('cpd_rwpn_20211118');lickwithout = false;
sse_set_vHP = cat(3,sse_set_vHP07,sse_set_vHP08,sse_set_vHP11,sse_set_vHP12,sse_set_vHP14);
sse_set_vHP  = sse_set_vHP(1:14,:,:);
load('CPD_no_speed');
sse_set_vHP_nospeed = cat(3,sse_set__nospeed_vHP07,sse_set__nospeed_vHP08,sse_set__nospeed_vHP11,sse_set__nospeed_vHP12,sse_set__nospeed_vHP14);


mouseG = {'dHP','vHP'};
day = datetime('today');
d = datestr(day,'yymmdd');
for ii = 1:2
%     eval(['sse_set = sse_set_',mouseG{ii},';']);
    eval(['sse_set = sse_set_',mouseG{ii},';']);
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

%psig= 1:13; xlimlist = repmat([1,120],13,1); wlist = repmat([8,6],13,1);
% termlist = {'Rw cue','Pn cue','Rw','Pn','Lick','Rw cue(t-1)','Pn cue(t-1)','Rw(t-1)','Pn(t-1)',...
%     'Rw cue(t-2)','Pn cue(t-2)','Rw(t-2)','Pn(t-2)','Speed'};
psig= [1 1 2 2 3 4, 6 6 7 7 8 8 9 9, 10 10 11 11 12 12 13 13]; % for figure
xlimlist = [1 50; 51 85; 1 50; 51 85; 51 85; 51 85;...% (t)
    1 50; 51 85; 1 50; 51 85; 1 50; 51 85; 1 50; 51 85;...% (t-1)
    1 50; 51 85; 1 50; 51 85; 1 50; 51 85; 1 50; 51 85;]; % (t-2) 
wlist = (diff(xlimlist')+1)*0.06; 

termlist = {'Rw cue1','Rw cue2','Pn cue1','Pn cue2','Rw','Pn',...
    'Rw cue(t-1)1','Rw cue(t-1)2','Pn cue(t-1)1','Pn cue(t-1)2','Rw(t-1)1','Rw(t-1)2','Pn(t-1)1','Pn(t-1)2',...
    'Rw cue(t-2)1','Rw cue(t-2)2','Pn cue(t-2)1','Pn cue(t-2)2','Rw(t-2)1','Rw(t-2)2','Pn(t-2)1','Pn(t-2)2'};
for isig = 1:length(psig)
    if ismember(psig(isig),[1,6,10]); ymax = 0.07;
    elseif ismember(psig(isig),[2,7,11]); ymax = 0.02;
    elseif ismember(psig(isig),[4,9,13]); ymax = 0.01;
    elseif psig(isig)==14; ymax = 0.1;
    else ymax = 0.03; end

%     f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 wlist(isig) 2.5]);%
    f1 = figure;
    hold on
    stdshade(cpd_dHP(:,:,psig(isig)), 0.3, [128 31 14]./255)%dorsal %brown
    stdshade(cpd_vHP(:,:,psig(isig)), 0.3, [8 128 75]./255)%ventral %green
    %     title(termlist{psig(isig)});
    plot((panovaset(:,psig(isig))<0.05)-1+0.99*ymax,'color','k','Marker','*','MarkerSize',1,'linestyle','none','linewidth',0.25)
    lim = axis;ax = gca; ax.TickLength = [0 0];
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([80 80], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    xticks([5 20 35 40 55 80]);xticklabels({}); % xticks([12 27 47 67 100]); xticklabels({'C', 'D1','D2','O','ITI'})
    xlim([xlimlist(isig,:)]); ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
    
    hold off
    saveas(f1,['E:\data\vHPC\all\figure\speed\speed',termlist{isig},d,'.tif'])
%     print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\cpd\rwpn\',termlist{isig},subt,num2str(d),'.ai']);
%     close
    
end

%% Reversal - cue period
clear
mouseID = input('mouseID = ');
load('behavior_align.mat','C_z','trial_80R_ind','trial_80P_ind','odorCue','outcomeContingency','rev_b','fr','trial_kind')
C_cue = mean(C_z([1:sum(rev_b), 420-sum(rev_b):end],0.5*fr+1:3.5*fr,:),2);
cueanova = trial_80R_ind([1:sum(rev_b), 420-sum(rev_b):end])+-1*trial_80P_ind([1:sum(rev_b), 420-sum(rev_b):end]); cueanova(cueanova==0)=[];
outanova = (odorCue([1:sum(rev_b), 420-sum(rev_b):end])==find(trial_kind==1)-1).*(outcomeContingency(([1:sum(rev_b), 420-sum(rev_b):end]),find(trial_kind==1)))+...
    (odorCue([1:sum(rev_b), 420-sum(rev_b):end])==find(trial_kind==4)-1).*(outcomeContingency(([1:sum(rev_b), 420-sum(rev_b):end]),find(trial_kind==4)));
outanova(outanova==0)=[];
clear panovaset
for ii = 1:size(C_cue,3)
    % [p,tbl,stats] = anovan(C_cue(trial_80R_ind|trial_80P_ind,:,ii),{cueanova, outanova},'varnames',{'Identity','Contingency'},...
    %     'display','off');
    [p,tbl,stats] = anovan(C_cue(trial_80R_ind([1:sum(rev_b), 420-sum(rev_b):end])|trial_80P_ind([1:sum(rev_b), 420-sum(rev_b):end]),:,ii),...
        {cueanova,outanova},'varnames',{'Identity','Contingency'},    'display','off');
    panovaset(ii,:) = p;
end
% endpien = [sum(panovaset(:,3)<0.05),sum(panovaset(:,3)>=0.05&panovaset(:,2)<0.05),...
%     sum(panovaset(:,3)>=0.05&panovaset(:,2)>=0.05&panovaset(:,1)<0.05), ...
%     sum(panovaset(:,3)>=0.05&panovaset(:,2)>=0.05&panovaset(:,1)>=0.05)];
endpien = [sum(panovaset(:,2)>=0.05&panovaset(:,1)<0.05), sum(panovaset(:,2)<0.05&panovaset(:,1)>=0.05),...
    sum(panovaset(:,2)<0.05&panovaset(:,1)<0.05), sum(panovaset(:,2)>=0.05&panovaset(:,1)>=0.05)];
eval(['endpien_',mouseID,' = endpien;']);
save('E:\data\vHPC\all\endpien_rwpn','-regexp','endpien_','-append')
%% figure each mouse
mouselist = {'endpien_dHP03','endpien_dHP04','endpien_dHP06','endpien_dHP08',...
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
saveas(f_cue,[figurepath,'Contingency_rwpn_',mouseID,'.tif'])
clear mouseID
end
%% figure group
load('E:\data\vHPC\all\endpien_rwpn')
endpien_dHP = sum([endpien_dHP03;endpien_dHP04;endpien_dHP06;endpien_dHP08]);
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
    saveas(f_cue,[figurepath,'Contingency_rwpn_',mouseG{ii},'.tif'])    
end

%% glm data
mouselist = {'glmresult_dHP03','glmresult_dHP04','glmresult_dHP06','glmresult_dHP07','glmresult_dHP08',...
    'glmresult_vHP06','glmresult_vHP07','glmresult_vHP08','glmresult_vHP10','glmresult_vHP11'};

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
figurepath = 'E:\data\vHPC\all\figure\glm_rwpn\glm_raw\';
yset = {'glmbetaset_dHP','glmbetaset_vHP','glmtvalueset_dHP','glmtvalueset_vHP'};
ylimset = [0 0.3; 0 1; 0 0.5; 0 0.1; 0 0.5;... %beta
    0 5; 0 5; 0 5;   0 2;   0 5];   %tvalue
% ylimset = [-0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1;... %beta
%     -1 1; -1 1; -1 1;   -1 1;   -1 1];   %tvalue
for HP = 1:4
    f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_outcome = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licks = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_lickm = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licke = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    eval(['ytmp=',yset{1,HP},';'])
    nst=1;
    for ii = 1:size(winstep,1)
        if contains(xlist{1,ii},'cue')
            cset = [0 0 205; 255 30 70; 80 80 80]./255;
            colortmp = cset(contains(xlist{1,ii},'rw')+2*contains(xlist{1,ii},'pn')+3*contains(xlist{1,ii},'0'),:);
            figure(f_cue); hold on; xticks([0 15 30 35]); xticks([0 15 30 35]);xticklabels({'0','1.5','3','3.5'});
            ylim(ylimset(floor((HP-1)/2)*5+1,:));
        elseif contains(xlist{1,ii},'o_rw')
            cset = [0 0 205; 129 212 250]./255;
            colortmp = cset(contains(xlist{1,ii},'x')+1,:);
            figure(f_outcome);hold on; xticks([0 25]); xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*5+2,:));
        elseif contains(xlist{1,ii},'o_pn')
            cset = [255 30 70; 239 154 154]./255;
            colortmp = cset(contains(xlist{1,ii},'x')+1,:);
            figure(f_outcome);hold on; xticks([0 25]); xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*5+2,:));
        elseif contains(xlist{1,ii},'o_no')
            colortmp = [80 80 80]./255;
            figure(f_outcome);hold on; xticks([0 25]); xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*5+2,:));
        elseif contains(xlist{1,ii},'lick')
            colortmp = [0 0 0]./255;
            if contains(xlist{1,ii},'onset'); figure(f_licks); ylim(ylimset(floor((HP-1)/2)*5+3,:));
                xticks([-5 0 10]);xticklabels({'-0.5','0','1'});
            elseif contains(xlist{1,ii},'mid'); figure(f_lickm); ylim(ylimset(floor((HP-1)/2)*5+4,:));
                xticks([0 20 40 60]);xticklabels({'0','2','4','6'});
            elseif contains(xlist{1,ii},'offset'); figure(f_licke); ylim(ylimset(floor((HP-1)/2)*5+5,:));
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
    
    figure(f_outcome); title('Outcome');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([25 25],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_outcome,[figurepath,yset{1,HP},'_outcome.tif'])
    
    
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
%% glm figure - diff
figurepath = 'E:\data\vHPC\all\figure\glm_rwpn\glm_raw\';
yset = {'glmbetaset_dHP','glmbetaset_vHP','glmtvalueset_dHP','glmtvalueset_vHP'};
% ylimset = [0 0.3; 0 0.3; 0 0.3; 0 0.1; 0 0.3;... %beta
%     0 3.5; 0 3.5; 0 3;   0 2;   0 3.5];   %tvalue
ylimset = [-0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1; -0.1 0.1;... %beta
    -1 1; -1 1; -1 1;   -1 1;   -1 1];   %tvalue
for HP = [1 3]
    f_cue = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_outcome = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licks = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_lickm = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    f_licke = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4]);
    eval(['dtmp=',yset{1,HP},';'])
    eval(['vtmp=',yset{1,HP+1},';'])
    nst=1;
    glmpvalueset = ones(size(dtmp,1),1);
    for ee = 1:size(dtmp,1)
        [h,p] = ttest2(dtmp(ee,:),vtmp(ee,:)');
        glmpvalueset(ee,1) = p;
    end
    for ii = 1:size(winstep,1)
        if contains(xlist{1,ii},'cue')
            figure(f_cue); hold on; xticks([0 15 30 35]); xticks([0 15 30 35]);xticklabels({'0','1.5','3','3.5'});
            ylim(ylimset(floor((HP-1)/2)*5+1,:));
        elseif contains(xlist{1,ii},'o_rw')
            figure(f_outcome);hold on; xticks([0 25]); xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*5+2,:));
        elseif contains(xlist{1,ii},'o_pn')
            cset = [255 30 70; 239 154 154]./255;
            colortmp = cset(contains(xlist{1,ii},'x')+1,:);
            figure(f_outcome);hold on; xticks([0 25]); xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*5+2,:));
        elseif contains(xlist{1,ii},'o_no')
            colortmp = [80 80 80]./255;
            figure(f_outcome);hold on; xticks([0 25]); xticks([0 25]);xticklabels({'0','2.5'});
            ylim(ylimset(floor((HP-1)/2)*5+2,:));
        elseif contains(xlist{1,ii},'lick')
            colortmp = [0 0 0]./255;
            if contains(xlist{1,ii},'onset'); figure(f_licks); ylim(ylimset(floor((HP-1)/2)*5+3,:));
                xticks([-5 0 10]);xticklabels({'-0.5','0','1'});
            elseif contains(xlist{1,ii},'mid'); figure(f_lickm); ylim(ylimset(floor((HP-1)/2)*5+4,:));
                xticks([0 20 40 60]);xticklabels({'0','2','4','6'});
            elseif contains(xlist{1,ii},'offset'); figure(f_licke); ylim(ylimset(floor((HP-1)/2)*5+5,:));
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
    
    figure(f_outcome); title('Outcome');ax = axis;
    line([0 0],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([25 25],[ax(1,3) ax(1,4)],'Color','k','LineStyle',':');
    line([ax(1,1) ax(1,2)],[0 0],'Color','k','LineStyle','--');
    saveas(f_outcome,[figurepath,yset{1,HP},'_outcome.tif'])
    
    
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
pvalue_set_dHP = cat(2,pvalue_dHP03,pvalue_dHP04,pvalue_dHP06,pvalue_dHP08);
src_set_dHP = cat(2,src_dHP03,src_dHP04,src_dHP06,src_dHP08);
panovaset_dHP = cat(1,panovaset_dHP03,panovaset_dHP04,panovaset_dHP06,panovaset_dHP08);
pvalue_set_vHP = cat(2,pvalue_vHP06,pvalue_vHP08,pvalue_vHP10,pvalue_vHP11);
src_set_vHP = cat(2,src_vHP06,src_vHP08,src_vHP10,src_vHP11);
panovaset_vHP = cat(1,panovaset_vHP06,panovaset_vHP08,panovaset_vHP10,panovaset_vHP11);
mouseG = {'dHP','vHP'};

colormap = [65 105 225; 255 99 71;138 43 226;]./255;
plotname = {'Rw°ÊPn cue','Pn°ÊRw cue','Rw','Pn'};
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
    plotset = [src(1,:,1)-src(2,:,1);src(1,:,2)-src(2,:,2)];
    scatter(plotset(1,sig_non),plotset(2,sig_non),10,[0.3 0.3 0.3],'o','filled')
    scatter(plotset(1,sig_cue),plotset(2,sig_cue),20,colormap(1,:),'o','filled')
    scatter(plotset(1,sig_contingency),plotset(2,sig_contingency),30,colormap(2,:),'o','linewidth',1.5)
    %         sig_before = pvalue(plotset(iset,1)+1,:,1)<0.05;
    %         sig_after = pvalue(plotset(iset,2)+1,:,2)<0.05;
    %         sig_non = pvalue(plotset(iset,1)+1,:,1)>0.05&pvalue(plotset(iset,2)+1,:,2)>0.05;
    %         [R ,p]= corrcoef(src(plotset(iset,1),:,1),src(plotset(iset,2),:,2))
    line([-1 1], [0 0],'Color','k','LineStyle',':')
    line([0 0],[-1 1],'Color','k','LineStyle',':')
    xlabel('Before reversal'); ylabel('After reversal')
    saveas(f,[cd,'\figure\SRCrev_',mouseG{ii},'_cuediff_',d,'.tif'])
    
    %         src_pie(:,ii,iset) = [sum(sig_before&sig_after),sum(sig_before&~sig_after),sum(~sig_before&sig_after),sum(~sig_before&~sig_after)];

end
%%
figurepath = 'E:\data\vHPC\all\figure\';
day = datetime('today');
d = datestr(day,'yymmdd');
plotset = [1,1;2,2;1,2;2,1];
plotname = {'Rw°ÊPn cue','Pn°ÊRw cue','Rw','Pn'};
for ii = 1:2
    for iset=1:4
    f = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
    pie(src_pie(:,ii,iset))
    colormap([138 43 226; 0 0 255 ;255 0 0;130 130 130;]./255)
    saveas(f,[figurepath,'SRCpie_',mouseG{ii},'_',plotname{iset},'_',d,'.tif'])
    end
end
%% SRC, beta delay vs. ITI
load('period_reg_rwpn.mat')
beta_set_dHP = cat(2,betaset_dHP03,betaset_dHP04,betaset_dHP06,betaset_dHP07,betaset_dHP08,betaset_dHP10);
pvalue_set_dHP = cat(2,pvalueset_dHP03,pvalueset_dHP04,pvalueset_dHP06,pvalueset_dHP07,pvalueset_dHP08,pvalueset_dHP10);
src_set_dHP = cat(2,srcset_dHP03,srcset_dHP04,srcset_dHP06,srcset_dHP07,srcset_dHP08,srcset_dHP10);
beta_set_vHP = cat(2,betaset_vHP06,betaset_vHP07,betaset_vHP08,betaset_vHP11,betaset_vHP12,betaset_vHP14);
pvalue_set_vHP = cat(2,pvalueset_vHP06,pvalueset_vHP07,pvalueset_vHP08,pvalueset_vHP11,pvalueset_vHP12,pvalueset_vHP14);
src_set_vHP = cat(2,srcset_vHP06,srcset_vHP07,srcset_vHP08,srcset_vHP11,srcset_vHP12,srcset_vHP14);
% 1st : (intercept) Rwcue Pncue Rw Pn lick speed Rwcue(t-1) Pncue(t-1)
%                   Rw(t-1) Pn(t-1) Rwcue(t-2) Pncue(t-2) Rw(t-2) Pn(t-2)
% 2nd : Cell
% 3rd : delay outcome outcome(1s)
termlist = [0; 2]; % cue, outcome

mouseG = {'dHP','vHP'};
cueG = {'Rwcue','Pncue'};
colormap = [65 105 225; 255 99 71;138 43 226;]./255;
day = datetime('today');
d = datestr(day,'yymmdd');

for iicue = 1:2
    f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 10]);
    for iiHP = 1:2
        eval(['betaset = beta_set_',mouseG{iiHP},';']);
        eval(['srcset = src_set_',mouseG{iiHP},';']);
        eval(['pvalueset = pvalue_set_',mouseG{iiHP},';']);
        for iterm = 1:2
            if iicue==2 & iterm == 2; outcome1s = 1; else; outcome1s = 0; end;
            if iicue ==1 & iterm ==1; xyrange = [-2 2]; else; xyrange = [-1 1]; end
            subplot(3,2,2*(iiHP-1)+iterm)
            axis square
            neuronset = [pvalueset(iicue+1,:,1)>=0.05&pvalueset(iicue+termlist(iterm,1)+1,:,2+outcome1s)>=0.05;...
                pvalueset(iicue+1,:,1)<0.05&pvalueset(iicue+termlist(iterm,1)+1,:,2+outcome1s)>=0.05;...
                pvalueset(iicue+1,:,1)>=0.05&pvalueset(iicue+termlist(iterm,1)+1,:,2+outcome1s)<0.05;...
                pvalueset(iicue+1,:,1)<0.05&pvalueset(iicue+termlist(iterm,1)+1,:,2+outcome1s)<0.05];
            hold on
            scatter(srcset(iicue,:,1),srcset(iicue+termlist(iterm,1),:,2+outcome1s),1,[0.5 0.5 0.5],'o','filled') %non sig
            scatter(srcset(iicue,neuronset(4,:) & (srcset(iicue,:,1) .* srcset(iicue+termlist(iterm,1),:,2+outcome1s) >0),1),...
                srcset(iicue+termlist(iterm,1),neuronset(4,:) & (srcset(iicue,:,1) .* srcset(iicue+termlist(iterm,1),:,2+outcome1s) >0),2),1,[1 0.27 0],'o','filled') % sig 1,3
            scatter(srcset(iicue,neuronset(4,:) & (srcset(iicue,:,1) .* srcset(iicue+termlist(iterm,1),:,2+outcome1s) <0),1),...
                srcset(iicue+termlist(iterm,1),neuronset(4,:) & (srcset(iicue,:,1) .* srcset(iicue+termlist(iterm,1),:,2+outcome1s) <0),2),1,[1 0.85 0],'o','filled') % sig 2,4
            line([-10 10], [0 0],'Color','k','LineStyle',':')
            line([0 0],[-10 10],'Color','k','LineStyle',':')
            xlim(xyrange); ylim(xyrange)
            [r,p]=corrcoef(srcset(iicue,:,1),srcset(iicue+termlist(iterm,1),:,2+outcome1s));
            pset(iiHP,iterm,iicue) = p(1,2);
            if p(1,2)<0.05; 
                mdl = fitlm(srcset(iicue,:,1),srcset(iicue+termlist(iterm,1),:,2+outcome1s));
                Coe = table2array(mdl.Coefficients);
                f = @(x) Coe(2,1)*x + Coe(1,1); ezplot( f, [xyrange, xyrange]); end
            text(0.5, 0.5,...
                {['R = ',num2str(round(r(1,2),4))],['p = ',num2str(round(p(1,2),4))]},...
                'Fontsize', 6);
            xlabel([]);
            title([mouseG{iiHP},' Delay - Outcome'],'fontsize',7)
            
            barsetsum(iiHP,2*(iterm-1)+1:2*(iterm-1)+2) = ...
                [sum(neuronset(4,:) & (srcset(iicue,:,1) .* srcset(iicue+termlist(iterm,1),:,2) >0)),sum(neuronset(4,:) & (srcset(iicue,:,1) .* srcset(iicue+termlist(iterm,1),:,2) <0))];
            
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
    barset(1,:) = barsetsum(1,:)./size(beta_set_dHP,2);
    barset(2,:) = barsetsum(2,:)./size(beta_set_vHP,2);
    subplot(3,2,5)
        axis square
        hold on
        bar(0.5,barset(1,1),'facecolor',[1 0.27 0],'Barwidth',0.2)
        bar(0.75,barset(1,2),'facecolor',[1 0.85 0],'Barwidth',0.2)
        bar(1.25,barset(2,1),'facecolor',[1 0.27 0],'Barwidth',0.2)
        bar(1.5,barset(2,2),'facecolor',[1 0.85 0],'Barwidth',0.2)
        xlim([0.2 1.8]); xticklabels({});
        
        subplot(3,2,6)
        axis square
        hold on
        bar(0.5,barset(1,3),'facecolor',[1 0.27 0],'Barwidth',0.2)
        bar(0.75,barset(1,4),'facecolor',[1 0.85 0],'Barwidth',0.2)
        bar(1.25,barset(2,3),'facecolor',[1 0.27 0],'Barwidth',0.2)
        bar(1.5,barset(2,4),'facecolor',[1 0.85 0],'Barwidth',0.2)
        xlim([0.2 1.8]); xticklabels({});
%         print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\zscore\src compare delay_',cueG{iicue},' ITI_outcome_',d,'.ai']);
%         saveas(f1,['E:\data\vHPC\all\figure\zscore\src compare delay_Rwcue ITI_outcome_',mouseG{iiHP},d,'.tif'])
end
%% RPE
%data settting
load('RPE_rwpn.mat')
betaset_dHP = cat(2,betaset_RPE_dHP03,betaset_RPE_dHP04,betaset_RPE_dHP06,betaset_RPE_dHP07,betaset_RPE_dHP08,betaset_RPE_dHP10);%
pvalueset_dHP = cat(2,pvalueset_RPE_dHP03,pvalueset_RPE_dHP04,pvalueset_RPE_dHP06,pvalueset_RPE_dHP07,pvalueset_RPE_dHP08,pvalueset_RPE_dHP10);%
srcset_dHP = cat(2,srcset_RPE_dHP03,srcset_RPE_dHP04,srcset_RPE_dHP06,srcset_RPE_dHP07,srcset_RPE_dHP08,srcset_RPE_dHP10);%
betaset_vHP = cat(2,betaset_RPE_vHP06,betaset_RPE_vHP07,betaset_RPE_vHP08,betaset_RPE_vHP11,betaset_RPE_vHP12,betaset_RPE_vHP14);%
pvalueset_vHP = cat(2,pvalueset_RPE_vHP06,pvalueset_RPE_vHP07,pvalueset_RPE_vHP08,pvalueset_RPE_vHP11,pvalueset_RPE_vHP12,pvalueset_RPE_vHP14);%
srcset_vHP = cat(2,srcset_RPE_vHP06,srcset_RPE_vHP07,srcset_RPE_vHP08,srcset_RPE_vHP11,srcset_RPE_vHP12,srcset_RPE_vHP14);%
mouseG = {'dHP','vHP'};
% 1st : (intercept) rwcue pncue rw pn ....
% 2nd : neuron
% 3rd : Outcome 1.5s, Outcome 1s

day = datetime('today');
d = datestr(day,'yymmdd');
f1=figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 10]);
%figure
for iHP = 1:2;    
    eval(['pvalueset_RPE = pvalueset_',mouseG{iHP},';']);
    eval(['srcset_RPE = srcset_',mouseG{iHP},';']);
    % rw value both sig in normal regression
    % (3rd 1: ITI, 2:Outcome 1.5s,)
    
    %rwtrial
    subplot(2,3,iHP)    
    axis square
    RPE_nomi_r=pvalueset_RPE(4,:,1)<0.05 & pvalueset_RPE(2,:,1)<0.05;
    quart1=srcset_RPE(3,:,1)>0 & srcset_RPE(1,:,1)>0;
    quart2=srcset_RPE(3,:,1)<0 & srcset_RPE(1,:,1)>0;
    quart3=srcset_RPE(3,:,1)<0 & srcset_RPE(1,:,1)<0;
    quart4=srcset_RPE(3,:,1)>0 & srcset_RPE(1,:,1)<0;
    quart{1,iHP}=[quart1;quart2;quart3;quart4];
    hold on
    scatter(srcset_RPE(3,:,1),srcset_RPE(1,:,1),1,[0.3 0.3 0.3],'o','filled')
    scatter(srcset_RPE(3,RPE_nomi_r&(quart1|quart3),1),srcset_RPE(1,RPE_nomi_r&(quart1|quart3),1),1,[0 0 205]/255,'o','filled') %1,3
    scatter(srcset_RPE(3,RPE_nomi_r&(quart2|quart4),1),srcset_RPE(1,RPE_nomi_r&(quart2|quart4),1),1,[100 185 255]/255,'o','filled') %2,4
    line([-2 2], [0 0],'Color','k','LineStyle',':')
    line([0 0],[-2 2],'Color','k','LineStyle',':')
    title([mouseG{iHP},' Reward'],'fontsize',7)
    xlabel('SRC for Reward'); ylabel('SRC for Cue');
    set(gca,'FontSize',7)
    
    %pntrial    
    subplot(2,3,3+iHP)
    axis square
    RPE_nomi_ur=pvalueset_RPE(5,:,2)<0.05 & pvalueset_RPE(3,:,2)<0.05;
    quart1=srcset_RPE(4,:,2)>0 & srcset_RPE(2,:,2)>0;
    quart2=srcset_RPE(4,:,2)<0 & srcset_RPE(2,:,2)>0;
    quart3=srcset_RPE(4,:,2)<0 & srcset_RPE(2,:,2)<0;
    quart4=srcset_RPE(4,:,2)>0 & srcset_RPE(2,:,2)<0;
    quart{2,iHP}=[quart1;quart2;quart3;quart4];
    hold on
    scatter(srcset_RPE(4,:,2),srcset_RPE(2,:,2),1,[0.3 0.3 0.3],'o','filled')
    scatter(srcset_RPE(4,RPE_nomi_ur&(quart1|quart3),2),srcset_RPE(2,RPE_nomi_ur&(quart1|quart3),2),1,[255 30 70]/255,'o','filled')
    scatter(srcset_RPE(4,RPE_nomi_ur&(quart2|quart4),2),srcset_RPE(2,RPE_nomi_ur&(quart2|quart4),2),1,[239 154 154]/255,'o','filled')
    line([-1 1], [0 0],'Color','k','LineStyle',':')
    line([0 0],[-1 1],'Color','k','LineStyle',':')
    title([mouseG{iHP},' Punishment'])
    xlabel('SRC for Punishment'); ylabel('SRC for Cue');
    set(gca,'FontSize',7)
    
    RPE_nomi{1,iHP}=[RPE_nomi_r; RPE_nomi_ur];    
end
% fraction of neuron of RPE
totNeu = [size(RPE_nomi{1,1},2),size(RPE_nomi{1,2},2)]; % total neuron
% RPE_nomi: dHP, vHP ; reward, punishment nomi
% quart: {(Rw, Pn),(dHP, vHP)};  1st,2nd,3rd,4th quart

bardataL = [sum(RPE_nomi{1,1}(1,:)&(quart{1,1}(1,:)|quart{1,1}(3,:))),sum(RPE_nomi{1,1}(1,:)&(quart{1,1}(2,:)|quart{2,1}(4,:))),...
    sum(RPE_nomi{1,2}(1,:)&(quart{1,2}(1,:)|quart{1,2}(3,:))),sum(RPE_nomi{1,2}(1,:)&(quart{1,2}(2,:)|quart{1,2}(4,:)));...
    sum(RPE_nomi{1,1}(2,:)&(quart{2,1}(1,:)|quart{2,1}(3,:))),sum(RPE_nomi{1,1}(2,:)&(quart{2,1}(2,:)|quart{2,1}(4,:))),...
    sum(RPE_nomi{1,2}(2,:)&(quart{2,2}(1,:)|quart{2,2}(3,:))),sum(RPE_nomi{1,2}(2,:)&(quart{2,2}(2,:)|quart{2,2}(4,:)))];
subplot(2,3,3) % Value update signal
axis square
hold on
bar(0.5,bardataL(1,1)./totNeu(1),'facecolor',[0 0 205]./255,'Barwidth',0.2)
bar(0.75,bardataL(1,2)./totNeu(1),'facecolor',[100 185 255]./255,'Barwidth',0.2)
bar(1.25,bardataL(1,3)./totNeu(2),'facecolor',[0 0 205]./255,'Barwidth',0.2)
bar(1.5,bardataL(1,4)./totNeu(2),'facecolor',[100 185 255]./255,'Barwidth',0.2)
xlim([0.2 1.8]); xticks([]); ylim([0 0.1]); yticks([0 0.1]);
set(gca,'FontSize',7)
% ylabel('Fraction of neuron')

subplot(2,3,6) % RPE
axis square
hold on
bar(0.5,bardataL(2,1)./totNeu(1),'facecolor',[255 30 70]./255,'Barwidth',0.2)
bar(0.75,bardataL(2,2)./totNeu(1),'facecolor',[239 154 154]./255,'Barwidth',0.2)
bar(1.25,bardataL(2,3)./totNeu(2),'facecolor',[255 30 70]./255,'Barwidth',0.2)
bar(1.5,bardataL(2,4)./totNeu(2),'facecolor',[239 154 154]./255,'Barwidth',0.2)
xlim([0.2 1.8]); xticks([]); ylim([0 0.1]); yticks([0 0.1]);
set(gca,'FontSize',7)
% ylabel('Fraction of neuron')

% text(3,0.1-iHP*0.02,{[mouseG{1,iHP}],['Total Neuron:', num2str(totNeu), '  Sig Neu:',num2str(sum(sum(bardataL(1:2,2:3))))],...
%     ['+RPE:', num2str(sum(bardataL(1,2:3))),'  -RPE:',num2str(sum(bardataL(2,2:3)))]})

% xticks([1 2 4 5 7 8]); xticklabels({'dHP-sig','vHP-sig','dHP-pos','dHP-neg','vHP-pos','vHP-neg'}); xlim([0.3 8.8])



% saveas(f1,['E:\data\vHPC\all\figure\RPE\rwprob_no_rev_',mouseG{iHP},'.tif'])
% print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\RPE\rwpn_',d,'.ai']);

