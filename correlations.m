clear
dirn = {'dHP03\rwprob','dHP04\rwprob','dHP06\rwprob','dHP07\rwprob','dHP08\rwprob','dHP10\rwprob',...
    'vHP06\rwprob','vHP07\rwprob','vHP08\rwprob','vHP11\rwprob','vHP12\rwprob','vHP14\rwprob',...
    'dHP03\rwpn02','dHP04\rwpn02','dHP06\rwpn03','dHP07\rwpn03','dHP08\rwpn03','dHP10\rwpn03',...
    'vHP06\rwpn03','vHP07\rwpn03','vHP08\rwpn03','vHP11\rwpn03','vHP12\rwpn03','vHP14\rwpn02'};

for imouse = 1:24
    cd(['E:\data\vHPC\',dirn{1,imouse}])
    mousename = strsplit(dirn{1,imouse},'\'); mousename = mousename{1,1};
    clear cylinderspeed
    % clearvars -except imouse dirn mousename
    load('behavior_align.mat','C_z_all','set_eventframe','timeline','fr','trial_kind','waterReward','odorCue','outcomeContingency','rwProb','trial_20R_ind')
    load('behavior_align.mat','cylinderspeed')
    if ~exist('cylinderspeed','var')
        cylinderspeed = zeros(420,130);
    end
    lickinfo = timeline{1,1};
    licksort = sortrows(lickinfo(:,1:3), 2);
    for iii = 1:licksort(end,2);
        licktemp = licksort((licksort(:,2)==iii),1);
        lick_reg_temp(iii,:) = histcounts(licktemp, 0:100000:12000000);
    end
    lick_sum = sum(lick_reg_temp,2);
    lick_reg = 2*movsum(lick_reg_temp,5,2); %movmean for 500ms
    lick_reg = movmean(lick_reg, 5,2);
    lick= reshape(lick_reg(:,41:80)',420*40,[]);
    cs = reshape(cylinderspeed(:,1:120)',420*120,[]);
    cs = movmean(cs,5);
  
    clear rew_regTmp pun_regTmp value_regTmp
    set_rectrial= 1:420;
    % cue, reward, punish
    outcome = waterReward;
    outcome(waterReward==0 & ismember(odorCue, find(trial_kind==1|trial_kind==2|trial_kind==4)-1)) =-1;
    for ii = 1:size(odorCue,1);
        rew_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==1)& (rwProb(ii,odorCue(ii)+1)==75);
        pun_regTmp(ii) = (outcomeContingency(ii,odorCue(ii)+1)==2)& (rwProb(ii,odorCue(ii)+1)==75);
        value_regTmp(ii,1) = rwProb(ii,odorCue(ii)+1)./100;
    end
    set_crp = [odorCue, outcome.*rew_regTmp', outcome.*pun_regTmp'];
    set_crp = set_crp(set_rectrial,:);
    if ismember(trial_kind,[0 1 3 4]);
        %     rp reg
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
        corset = [cue1_t1_reg(:,1),cue1_t2_reg(:,1),cue2_t1_reg(:,1),cue2_t2_reg(:,1),rew_1_reg(:,1),rew_2_reg(:,1),pun_1_reg(:,1),pun_2_reg(:,1)];
        for icor = 1:size(corset,2)
            [r,p] = corrcoef(lick_sum,corset(:,icor)); % lick vs v(t-1)
            temp_vif = diag(inv(r))';
            r_lick_2(imouse-12, icor) = r(1,2); p_lick_2(imouse-12, icor) = p(1,2);
        end
    elseif ismember(trial_kind,[0 1 2 3]);
        %     prob reg
        outcome = waterReward;
        outcome(waterReward==0)=-1;
        outcome(trial_20R_ind)=0;
        value_regTmp= value_regTmp(randperm(420)); %for scramble
        set_crp = [outcome, value_regTmp];
        set_crp = set_crp(set_rectrial,:);
        outcome_reg = repmat(set_crp(:,1), 1, 120);
        outcome_1_reg = repmat([0; set_crp(1:end-1,1)],1,120);
        outcome_2_reg = repmat([0;0; set_crp(1:end-2,1)],1,120);
        value_reg = repmat(set_crp(:,2), 1, 120);
        value_1_reg = repmat([0; set_crp(1:end-1,2)],1,120);
        value_2_reg = repmat([0;0; set_crp(1:end-2,2)],1,120);
        corset = [value_1_reg(:,1),value_2_reg(:,1),outcome_1_reg(:,1),outcome_2_reg(:,1)];
        for icor = 1:size(corset,2)
            [r,p] = corrcoef(lick_sum,corset(:,icor)); % lick vs v(t-1)
            r_lick_1(imouse, icor) = r(1,2); p_lick_1(imouse, icor) = p(1,2);
        end
    end
    
    
%         ctmp=movmean(C_raw,15,2);
%         ctmp = ctmp(:,set_eventframe(:));
%         ctmp = reshape(ctmp', [],12*fr,size(ctmp,1));
%     C_raw = movmean(C_raw, 15,2);
%     C_z_all = (C_raw-mean(C_raw(:,set_eventframe(1,1):set_eventframe(420,360)),2))./std(C_raw(:,set_eventframe(1,1):set_eventframe(420,360)),0,2);
%     C_z_all = C_z_all(:,set_eventframe(:));
%     C_z_all = reshape(C_z_all', [],12*fr,size(C_z_all,1));
    ctmp = movmean(C_z_all,15,2);
    craw_lick = permute(ctmp(:,4*fr+1:3:8*fr,:),[3 2 1]);
    craw_lick = reshape(craw_lick,size(craw_lick,1),[420*40]);
    craw_lick = craw_lick';
    
    craw_speed = permute(ctmp(:,2:3:end,:),[3 2 1]);
    craw_speed = reshape(craw_speed,size(craw_speed,1),[420*120]);
    craw_speed = craw_speed';
    
    for oo = 1:size(craw_speed,2)
        [r,p] = corrcoef(craw_speed(:,oo),cs); %speed vs c_raw correlation
        rset(oo,1) = r(1,2);
        pset(oo,1) = p(1,2);
        [r,p] = corrcoef(craw_lick(:,oo),lick); % lick vs c_raw correlation
        rset(oo,2) = r(1,2);
        pset(oo,2) = p(1,2);
    end
    rset_all{imouse} = rset;
    pset_all{imouse} = pset;
    clear rset pset

end
save('E:\data\vHPC\all\speed\All_correlation.mat','rset_all','pset_all','r_lick_1','p_lick_1','r_lick_2','p_lick_2','dirn')
    
%% figure
% load('E:\data\vHPC\all\speed\All_correlation.mat')

rset_mat = cell2mat(rset_all');
figure
subplot(2,1,1)
histogram(rset_mat(:,2),-1:2/51:1,'edgecolor','k','facecolor','w')
xticks([-1 -0.5 0 0.5 1])

subplot(2,1,2)
histogram(rset_mat(:,1),-1:2/51:1,'edgecolor','k','facecolor','w')
xticks([-1 -0.5 0 0.5 1])


%%
load('E:\data\vHPC\all\speed\vif.mat')
task = {'rwpn','rwprob'};
mlist = {'dHP07','dHP08','dHP10','vHP07','vHP08','vHP11','vHP12','vHP14'};
for itask = 1:2
    if itask==1; ind = [5 6];
    elseif itask==2; ind=[3 4];
    end
    for ilist = 1:8
        eval(['temp= max(vif_',task{itask},'_',mlist{ilist},',[],2);'])
        maxlist(ilist+8*(itask-1),:)=temp(ind);
    end
end

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 2.5]);
hold on
for ii =1:2    
    bar(ii,mean(maxlist(:,ii)),'barwidth',0.8,'facecolor','w','edgecolor','k')
    scatter(ii-0.3+0.6*rand(size(maxlist(:,ii),1),1),maxlist(:,ii),3,[0.5 0.5 0.5],'o','MarkerFaceColor',[0.5 0.5 0.5])
    errorbar(ii,mean(maxlist(:,ii)),sem(maxlist(:,ii)),'color','k','capsize',5)
end
% plot([0.77 1.77],[lickmean(indexHP==1); lickmean(indexHP==2)],'color',cmap1(1,:))
% plot([1.23 2.23],[lickmean(indexHP==3); lickmean(indexHP==4)],'color',cmap1(2,:))
xlim([0.4 2.6]); xticks([1 2]); xticklabels({})
ylim([0 12]); yticks([0 5 10 15]); yticklabels([0 0 0 0])
print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\speed\speedvif.ai']);

%%
load('rwpn_fon.mat')
load('rwprob_fon.mat')
pvaluesig = cat(3,pvaluesig_t1_dHP,pvaluesig_t2_dHP,pvaluesig_t1_vHP,pvaluesig_t2_vHP);
pvaluesig_r= cat(3,pvaluesig_r_t1_dHP,pvaluesig_r_t2_dHP,pvaluesig_r_t1_vHP,pvaluesig_r_t2_vHP);

FON= mean(pvaluesig, 3);
bin_percent = mean(pvaluesig_r, 3);
chis_HP=sum(pvaluesig, 3);

psig= [1 2 2]; % for figure
xlimlist = [40 80; 1 50; 51 80;];
wlist = (diff(xlimlist')+1)*0.055; %0.075 0.055
termlist = {'Lick(t)','Speed(t)1','Speed(t)2'};
cmap = [0 0 0; 0 0 0];
ymax = 0.4;
for isig = 1:3
    ii = psig(isig);
    f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 wlist(isig) 2.5]);%
    hold on
    plot(FON(ii,:), 'color', cmap(1,:), 'linewidth', 0.5);
    xticks([5 20 35 40 55 80]);xticklabels({});
    xlim([xlimlist(isig,:)]); ylim([0 ymax]); yticks([0 ymax]); yticklabels({});
    lim = axis;
    line([5 5], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %0.5 baseline
    line([20 20], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([35 35], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    line([40 40], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %1.5 cue
    line([55 55], [lim(3) lim(4)],'color','k', 'linestyle', ':'); %2.5 outcome
    
    for tt = 1:120
        sigchis_HP(ii,tt) = chis([sum(pvaluesig(ii,tt,:), 3) size(pvaluesig,3)-sum(pvaluesig(ii,tt,:), 3);...
            sum(pvaluesig_r(ii,tt,:), 3) size(pvaluesig_r,3)-sum(pvaluesig_r(ii,tt,:), 3)]);
    end
    FON_sig_HP = sigchis_HP<0.05;
    plot(FON(ii,:), 'color', cmap(1,:), 'linestyle','none','markersize',4,'marker','.','markerindices',find(FON_sig_HP(ii,:)));
    plot(bin_percent(ii,:),'color',cmap(1,:), 'linestyle', ':','linewidth', 0.5)
    ylim([0 ymax])
    
    %           saveas(f1,['E:\data\vHPC\all\figure\FON\rwprob\',termlist{isig},deltaname,lickname,'.tif'])
    print(f1,'-depsc','-painters',['E:\data\vHPC\all\figure\speed\',termlist{isig},'.ai']);
    
end









