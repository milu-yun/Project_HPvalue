mintrial = floor(120*0.7);
mincell = floor(42*0.7);
niter = 100;
dirn = {'dHP03_rwpn','dHP04_rwpn','dHP06_rwpn','dHP07_rwpn','dHP08_rwpn','dHP10_rwpn',...
    'vHP06_rwpn','vHP07_rwpn','vHP08_rwpn','vHP11_rwpn','vHP12_rwpn','vHP14_rwpn',...
    'dHP03_rwprob','dHP04_rwprob','dHP06_rwprob','dHP07_rwprob','dHP08_rwprob','dHP10_rwprob',...
    'vHP06_rwprob','vHP07_rwprob','vHP08_rwprob','vHP11_rwprob','vHP12_rwprob','vHP14_rwprob'};
for imouse = 1:24
    disp(['imouse = ',num2str(imouse)])
    clearvars -except imouse mintrial mincell niter dirn
    mousename = strsplit(dirn{1,imouse},'='); mousename = mousename{1,1};
    load(dirn{1,imouse},'C_raw','set_eventframe','fr','trial_80R_ind','trial_50R_ind','trial_20R_ind','trial_80P_ind')
    
    C_raw_mean = movmean(C_raw,15,2);
    C_event = C_raw_mean(:,set_eventframe(:));
    C_event = reshape(C_event', [],12*fr,size(C_event,1));
    C_event = permute(C_event,[1 3 2]);
    C_cue = mean(C_event(:,:,0.5*fr+1:2*fr),3);
    C_delay = mean(C_event(:,:,2*fr+1:3.5*fr),3);
    
    odor1 = find(trial_80R_ind);
    odor2 = [find(trial_50R_ind),find(trial_80P_ind)];
    odor3 = find(trial_20R_ind);
    odorset = repmat([1 2 3],83,1); odorset = reshape(odorset, 83*3,1);
    decodingresult_period = NaN(3,2,mintrial,niter);
    
    parfor iter = 1:niter
        randcell = randperm(size(C_raw,1),mincell);
        randodor1 = randperm(size(odor1,1),mintrial);
        randodor2 = randperm(size(odor2,1),mintrial);
        randodor3 = randperm(size(odor3,1),mintrial);
        randodorset = [odor1(randodor1),odor2(randodor2),odor3(randodor3)];
        for itest = 1:mintrial
            testodor = randodorset(itest,:);
            trainodor = randodorset;
            trainodor(itest,:)=[];
            trainodor = reshape(trainodor,(mintrial-1)*3,1);
            for tt = 1:2
                if tt==1; C_train= C_cue(trainodor,randcell); C_test = C_cue(testodor,randcell);
                elseif tt==2; C_train= C_delay(trainodor,randcell); C_test = C_delay(testodor,randcell);
                end
                Mdl1 = fitcecoc(C_train,odorset); %cue
                for icue = 1:3
                    label = predict(Mdl1,C_test(icue,:));
                    decodingresult_period(icue,tt,itest,iter) = label;
                end
            end
        end
    end
    
    save ([dirn{1,imouse},'_decoding_period.mat'],'decodingresult_period')
end