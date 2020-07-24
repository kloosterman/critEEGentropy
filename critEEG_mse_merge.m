function mseavg = critEEG_mse_merge(PREIN)

if ismac
  %   basepath = '/Users/kloosterman/gridmaster2012/kloosterman';
  basepath = '/Users/kloosterman/beegfs';
else
  basepath = '/home/beegfs/kloosterman'; % on the cluster
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJ  = {
  'Sub1'     ...
  'Sub2'    'Sub3'    ...
  'Sub4'   ...
  'Sub5'    'Sub6'    ...
  'Sub7'       'Sub8'  'Sub9' ...
  'Sub10'   ...
  'Sub11' ...
  'Sub12'    'Sub13'      'Sub14'   ...
  'Sub15'   ...
  'Sub16'
  }; %
nsub = length(SUBJ);

% triggers = {'stim' 'resp'};
triggers = {'stim' };
timelim = [-0.5 1; ... % stim  % timelim = [-0.5 1.5]; % for time courses
  -0.75 0.5]; % resp
% timelim = [-0.25 0.6; ... % stim  % timelim = [-0.5 1.5]; % for time courses
%   -0.25 0.25]; % resp

if nargin == 0
  % PREIN = '/Users/kloosterman/beegfs/projectdata/critEEG/mse_ERP-percond_robdetrord2/filtskip_r1_delERPyes';
  % PREIN = '/Users/kloosterman/beegfs/projectdata/critEEG/mse_robdetrord2_NoErpremoval/filtskip_r1_delERPno';
  % PREIN = '/Users/kloosterman/beegfs/projectdata/critEEG/mse_regressERP_robdetrord2/filtskip_r1_delERPyes';
  
  % PREIN = '/Users/kloosterman/beegfs/projectdata/critEEG/mse_robdetrord0_eqtrldur_eqtrlcnt/filtskip_r1_delERPno';
  % PREIN = '/Users/kloosterman/beegfs/projectdata/critEEG/mse_robdetrord1_eqtrldur/filtskip_r1_delERPno';
  
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/hpf_eqtrldur/filtskip_r_perscale_toi_sp_delERPyes';
  
  % % % % USED IN eLife 1st submission:
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes';
  
%   % % % % % Not removed ERP, per reviewer req, in elife paper2
  PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_';

  % % % % % % Theta (2-6 Hz) bandstop filtered
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_bs2-6hz';
%   % % % % % % Theta (7-11 Hz) bandstop filtered
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_bs7-11hz';

%   % % % % % % Theta (12-16 Hz) bandstop filtered
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_bs12-16hz';

%   % % % % % % 1 Hz highpass filtered
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf1_Hz';
%   % % % % % % 1 Hz phase 
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_1hzphase';
%   % % % % % % 2 Hz highpass filtered
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf2_hp2hz';
%   % % % % % % 3 Hz highpass filtered
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf3_Hz';
  
  % % % % % % Theta filtered, per reviewer request
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_filt2-6hz';
  
  % % % % % % 7-11 Hz filtered
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_filt7-11hz';
  
  % % % % % % Theta phase, out of interest
  PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_4hzphase';
  
  % % %  splithalf control analysis odd even, noERP:
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_alltrials_vardur_ntrlnotmatched_hpf0_oddeven/even';
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_alltrials_vardur_ntrlnotmatched_hpf0_oddeven/odd';
  
%   % % %  splithalf control analysis odd even, yesERP:
%   PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_oddeven/even';
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPno_alltrials_vardur_ntrlnotmatched_hpf0_oddeven/odd';

  % % only hits, see if it changes
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_onlyhits';
  
  % % only hits, fixed trl length, equal trial N
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_onlyhits_fixeddur_matchedntrl';
  % % point averaging. fixed r: control analysis
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/pointavg_r_per_toi_delERPyes';
  
  % % rpara = 0.1, vartrl length, vartrialN
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_alltrials_vardur_ntrlnotmatched_r0.1';
  
  % % 60 Hz highpassfilter to look at gamma range entropy r=0.5
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_alltrials_vardur_ntrlnotmatched_hpf60';
  
  % 60 Hz highpassfilter to look at gamma range entropy, r=0.1
  % PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/mse/hpf_vartrldur/filtskip_r_perscale_toi_sp_delERPyes_alltrials_vardur_ntrlnotmatched_hpf600.1';
end
cd(PREIN);

scalelim = [1 165]; % scale 1-42
mseavg.dat = [];
behavperrun.criterion = nan(nsub, 4, 11); %sub block cond
behavperrun.dprime = nan(nsub, 4, 11); %sub block cond

for itrig = 1:length(triggers)
  % load example mse to get nchan etc
  runname = sprintf('*_ses1_mse_icond%d_run%d_%s.mat', 1, 1, triggers{itrig});
  run = dir(runname);
  load(run(1).name)
  
  tind = mse.time >= timelim(itrig,1) & mse.time <= timelim(itrig,2);
  ntim = size(find(tind),2);
  
  nchan = size(mse.sampen,1); % should be the same for stim and resp
  scind = mse.timescales >= scalelim(1) & mse.timescales <= scalelim(2);
  nscales = size(scind,2);
  
  mseavg.dat{itrig} = nan(nsub,7, nchan,nscales,ntim, 4,11); % subj meas    toi channel scale    cond run
  mseavg.time{itrig} = mse.time(tind);
  mseavg.timescales = mse.timescales(scind);
  
  if itrig == 1
    mseavg.baseline = nan(nsub,nchan,nscales,4,11); % subj
  end
  for icond = 1:2 % 1:2 % crit
    for isub = 1:length(SUBJ)
      
      runname = sprintf('%s_ses*_mse_icond%d_run*_%s.mat', SUBJ{isub}, icond,  triggers{itrig});
      
      run = dir(runname);
      
      if isempty(run)
        fprintf('%s not found\n', runname)
        continue
      end
      
      for irun = 1:length(run)
        
        disp(run(irun).name)
        load(run(irun).name)
        
        trialinfo = mse.trialinfo; % data.trialinfo(mse.cfg.trials,:);
        chronolist = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10]; % 1 and 2 always together etc
        runno = chronolist(trialinfo(1,8)); % runno per condition
        if itrig == 1
          
          % compute behavior per run
          hitrate = sum(trialinfo(:,2) == 1 & trialinfo(:,3) == 1) / ...  % Fig report
            sum(trialinfo(:,2) == 1); % N fig trials
          farate = sum(trialinfo(:,2) == 2 & trialinfo(:,3) == 1) / ...  % Hom report
            sum(trialinfo(:,2) == 2); % N hom trials
          % correct 0 and 1 FA/H rates
          if hitrate == 1;  hitrate = 1 - 1/(2*sum(trialinfo(:,2) == 1));   end
          if farate == 0;   farate  = 1 / (2*sum(trialinfo(:,2) == 2));  end
          if hitrate == 0;  hitrate = 1 / (2*sum(trialinfo(:,2) == 1));   end
          if farate == 1;   farate  = 1 - 1 / (2*sum(trialinfo(:,2) == 2));  end
          
          behavperrun.criterion(isub,icond,runno) = -0.5 * (norminv(hitrate) + norminv(farate));
          behavperrun.dprime(isub,icond,runno) = norminv(hitrate) - norminv(farate);
          
          %                     basetind = mse.time >= -0.45 & mse.time <= -0.25;
          %                     basetind = mse.time >= -0.4 & mse.time <= -0.2;
          basetind = mse.time >= -0.2 & mse.time <= 0;
          baseline = squeeze(mean(mse.sampen(:,:, basetind),3));
          mseavg.baseline(isub,:,:,icond,runno) = baseline; % sub ses chan scales cond stim resp
          baseline_r = squeeze(mean(mse.r(:,:, basetind),3));
          
        else
          %           basetind = mse.time >= -0.75 & mse.time <= -0.5;
          %           %           basetind = mse.time >= -0.25 & mse.time <= 0;
          %           baseline = squeeze(mean(mse.sampen(:,:, basetind),3));
          % take corresponding stimlocked
          baseline = squeeze(mseavg.baseline(isub,:,:,icond,runno)); % sub ses chan scales cond stim resp
        end
        try
          mseavg.dat{itrig}(isub,1, :,:,:,  icond, runno) = mse.sampen(:,scind,tind) ; % chan_scale_time now
          mseavg.dat{itrig}(isub,2, :,:,:,  icond, runno) = ( mse.sampen(:,scind,tind) - baseline) ./ baseline * 100; % chan_scale_time now
          mseavg.dat{itrig}(isub,3, :,:,:,  icond, runno) = mse.r(:,scind,tind) ; % chan_scale_time now
          mseavg.dat{itrig}(isub,4, :,:,:,  icond, runno) = ( mse.r(:,scind,tind) - baseline_r) ./ baseline_r * 100; % chan_scale_time now
        catch
          disp('Not enough timepoints')
        end
        
      end
    end
    
    mseavg.dat{itrig}(:,:,:, :,:,icond, 11) = nanmean(mseavg.dat{itrig}(:,:,:, :,:,icond, 1:10),7); % subj raworpsc    toi channel scale    cond run
    
    %     disp('control MSE for r by regression across subjects')
    %     msedat = reshape(squeeze( mseavg.dat{itrig}(:,1, :,:,:,  icond, 11) ), nsub, []);
    %     rdat =   reshape(squeeze( mseavg.dat{itrig}(:,3, :,:,:,  icond, 11) ), nsub, []);
    %     temp = nan(size(rdat));
    %     for i = 1:size(msedat,2)
    % %       [~,~,residual] = regress(msedat(:,i), rdat(:,i) );
    %       [~,~,residual] = regress(msedat(:,i), [rdat(:,i) ones(size(rdat(:,i)))] ); % add intercept
    %       temp(:,i) = residual;
    %     end
    %     temp = squeeze(reshape(temp, size(mseavg.dat{itrig}(:,5,:, :,:,icond, 11 ))));
    %     mseavg.dat{itrig}(:,5,:, :,:,icond, 11 ) = temp;
    %
    %     % do the same for psc mse vs psc r
    %     msedat = reshape(squeeze( mseavg.dat{itrig}(:,2, :,:,:,  icond, 11) ), nsub, []);
    %     rdat =   reshape(squeeze( mseavg.dat{itrig}(:,4, :,:,:,  icond, 11) ), nsub, []);
    %     temp = nan(size(rdat));
    %     for i = 1:size(msedat,2)
    % %       [~,~,residual] = regress(msedat(:,i), rdat(:,i) );
    %       [~,~,residual] = regress(msedat(:,i), [rdat(:,i) ones(size(rdat(:,i)))] ); % add intercept
    %       temp(:,i) = residual;
    %     end
    %     temp = squeeze(reshape(temp, size(mseavg.dat{itrig}(:,6,:, :,:,icond, 11 ))));
    %     mseavg.dat{itrig}(:,6,:, :,:,icond, 11 ) = temp;
  end
  
  mseavg.dat{itrig}(:,:,:, :,:,3, : ) = nanmean(mseavg.dat{itrig},6); % avg over cond
  mseavg.dat{itrig}(:,:,:, :,:,4, : ) = mseavg.dat{itrig}(:,:,:, :,:,2,:) - mseavg.dat{itrig}(:,:,:, :,:, 1,:); %lib-cons
  
  disp('control lib-cons MSE for r by regression across subjects')
  for irun = 11
    msedat = reshape(squeeze( mseavg.dat{itrig}(:,1, :,:,:,  4, irun) ), nsub, []);
    rdat =   reshape(squeeze( mseavg.dat{itrig}(:,3, :,:,:,  4, irun) ), nsub, []);
    temp = nan(size(rdat));
    temp2 = nan(size(rdat));
    for i = 1:size(msedat,2)
      [~,~,residual] = regress(msedat(:,i), [rdat(:,i) ones(size(rdat(:,i)))] ); % add intercept
      [~,~,residual2] = regress(rdat(:,i), [msedat(:,i) ones(size(rdat(:,i)))] ); % also control r for mse
      temp(:,i) = residual;
      temp2(:,i) = residual2;
    end
    temp = squeeze(reshape(temp, size(mseavg.dat{itrig}(:,5,:, :,:,4, irun ))));
    temp2 = squeeze(reshape(temp2, size(mseavg.dat{itrig}(:,5,:, :,:,4, irun ))));
    mseavg.dat{itrig}(:,5,:, :,:,4, irun ) = temp;
    mseavg.dat{itrig}(:,7,:, :,:,4, irun ) = temp2; %mse-controlled r
    
    % do the same for psc mse vs psc r
    msedat = reshape(squeeze( mseavg.dat{itrig}(:,2, :,:,:,  4, irun) ), nsub, []);
    rdat =   reshape(squeeze( mseavg.dat{itrig}(:,4, :,:,:,  4, irun) ), nsub, []);
    temp = nan(size(rdat));
    %     temp2 = nan(size(rdat)); % could do the same for psc
    for i = 1:size(msedat,2)
      [~,~,residual] = regress(msedat(:,i), [rdat(:,i) ones(size(rdat(:,i)))] ); % add intercept
      %       [~,~,residual2] = regress(rdat(:,i), [msedat(:,i) ones(size(rdat(:,i)))] ); % add intercept
      temp(:,i) = residual;
      %       temp2(:,i) = residual2;
    end
    temp = squeeze(reshape(temp, size(mseavg.dat{itrig}(:,6,:, :,:,4, irun ))));
    %     temp2 = squeeze(reshape(temp2, size(mseavg.dat{itrig}(:,6,:, :,:,4, irun ))));
    mseavg.dat{itrig}(:,6,:, :,:,4, irun ) = temp;
    %     mseavg.dat{itrig}(:,8,:, :,:,4, irun ) = temp;
  end
end

mseavg.fsample = mse.fsample;
mseavg.SUBJ = SUBJ;
% mseavg.toi_leg = {'Presti`m' 'Poststim' 'Post-Prestim_raw' 'Post-Prestim_psc'};
mseavg.mse_leg = {'MSEraw' 'MSEpsc' 'r_parameter' 'r_para_psc' 'MSEraw r removed' 'MSEpsc r removed' 'r MSE removed'};
mseavg.dimord = 'subj_msetype_chan_scale_time_cond_run';
mseavg.sens = critEEG_sensorselection();
mseavg.label = mseavg.sens.label;

mseavg.behav_conds = {'Conservative' 'Liberal' 'allbehavconds' 'Lib-Cons'};
mseavg.cfg = mse.cfg;
mseavg.dimordsize = size(mseavg.dat);

behavperrun.criterion(:,4,:) = behavperrun.criterion(:,2,:) - behavperrun.criterion(:,1,:);
behavperrun.dprime(:,4,:) = behavperrun.dprime(:,2,:) - behavperrun.dprime(:,1,:);
mseavg.behavperrun = behavperrun;


% put behavior in
behavin = fullfile('/Users/kloosterman/gridmaster2012/kloosterman', 'projectdata', 'critEEG', 'behavior');
% behavin = fullfile('/Users/kloosterman/gridmaster2012/kloosterman', 'projectdata', 'critEEG', 'behavior_noartfremoved');
% behavin = fullfile('/Users/kloosterman/gridmaster2012/kloosterman', 'projectdata', 'critEEG', 'behavior_noartfremoved_incsoftFAs');
load(fullfile(behavin, 'behavstruct.mat'));
mseavg.behavior = behav;

load(fullfile(behavin, 'ntrlpercond.mat'));
mseavg.behavior.ntrlpercond = ntrlpercond;

% % put reward in
% mseavg.behavior.reward = critEEG_compute_reward()

% PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/plots/mse';
% PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/plots';
mseavg.PREOUT = fullfile(PREIN, 'plots');
mkdir(mseavg.PREOUT)
