function critEEG_entropy(cfg)

ft_defaults
global ft_default
ft_default.checkpath = 'no'; % once

disp(cfg.fileload);
load(cfg.fileload); % data comes out
disp(cfg.outputfile);

fprintf('throw out trials without baseline, as in power paper\n')
cfg2 = [];
cfg2.vartrllength = 2;
cfg2.keeptrials = 'yes';
timelock = ft_timelockanalysis(cfg2, data);
tind = timelock.time <= -0.4;
cfg2 = [];
cfg2.trials = any(~isnan( timelock.trial(:,1,tind) ),3);
data = ft_selectdata(cfg2, data);
clear timelock tind

if strcmp(cfg.trialduration, 'fixeddur') % lib trials are often shorter
  cfg2 = [];
  cfg.toilim = [-0.5 1];
  cfg2.toilim = cfg.toilim ;  
  data = ft_redefinetrial(cfg2, data);
  cfg2 = [];
  cfg2.minlength = diff(cfg.toilim);
  data = ft_redefinetrial(cfg2, data);
end

if strcmp(cfg.removeERP, 'yes')  % subtract ERP per condition
  data_ERP1 = [];
  cfg3=[];
  for icond = 1:2
    cfg3.trials = data.trialinfo(:,1) == icond;
    tempdata = ft_selectdata(cfg3, data);
    data_ERP1{icond} = remove_ERP_fromdata(tempdata, cfg.removeERPmethod); % TODO add remove type to cfg
  end
  data = ft_appenddata([], data_ERP1{:});
  data.cfg.previous = data.cfg; % log detrending done in cfg to keep overview
  data.cfg = keepfields(data.cfg, 'previous');
  data.cfg.funthatwasrun = 'remove_ERP_fromdata';
  data.cfg.method = [cfg.removeERPmethod '_percondition'];
  clear data_ERP1
end

% if strcmp(cfg.robustdetrend, 'yes')
% %   order = 2; % 1 = linear fit, 2 = quadratic
%   for itrial=1:length(data.trial)
%     data.trial{itrial} = transpose(nt_detrend(data.trial{itrial}', cfg.robustdetrendorder)); 
% %   for itrial=1:20:length(data.trial)
% %     f = figure; f.Position = [        680          65        1000        1000];
% %     nt_detrend(data.trial{itrial}', cfg.robustdetrendorder); 
%   end
%   data.cfg.previous = data.cfg; % log detrending done in cfg to keep overview
%   data.cfg = keepfields(data.cfg, 'previous');
%   data.cfg.funthatwasrun = 'nt_detrend';
%   data.cfg.order = cfg.robustdetrendorder;
% end

% if strcmp(cfg.hpfilter, 'yes')
%   cfg2=[];
%   cfg2.hpfilter = 'yes';
%   cfg2.hpfreq = 1;
%   cfg2.padding = 7;
%   data = ft_preprocessing(cfg2, data);
% end

if strcmp(cfg.trigger, 'resp')
  cfg_resp=[];
  cfg_resp.offset = -data.trialinfo(:,4);
  cfg_resp.trials = find(cfg_resp.offset < 1);
  data = ft_redefinetrial(cfg_resp, data);
end

if strcmp(cfg.SDTtype, 'onlyhits')
  cfg2=[]; % select only hits
  cfg2.trials = data.trialinfo(:,2) == 1 & data.trialinfo(:,3) == 1; % stim 1 and resp 1
  data = ft_selectdata(cfg2, data);
  cfg.trials = 'all';
end

cfg.trials = data.trialinfo(:,8) == cfg.run; % or...:
if strcmp(cfg.trlcounts, 'matchedntrl')
  chronolock = [2 1 4 3 6 5 8 7 10 9 12 11 14 13 16 15 18 17];
  ctrpart = chronolock(cfg.run);
  ntrl1 = length(find(data.trialinfo(:,8) == cfg.run));
  ntrl2 = length(find(data.trialinfo(:,8) == ctrpart));
  
  cfg.trials = find(data.trialinfo(:,8) == cfg.run);
  if ntrl1 > ntrl2
    cfg.trials = cfg.trials(1:ntrl2);
  end
  data.cfg.previous = data.cfg; % log to keep overview
  data.cfg = keepfields(data.cfg, 'previous');
  data.cfg.funthatwasrun = 'matchedntrl';
end

if length(find(cfg.trials)) < 1
  warning('No trials remain for this condition, skipping')
  return
end

if cfg.hpfilter > 0
  cfg2 = [];
  cfg2.hpfilter = 'yes';
  cfg2.hpfreq = cfg.hpfilter;
  cfg2.padtype = 'mirror';
  cfg2.padding = 4;
  data = ft_preprocessing(cfg2, data);
end
if isfield(cfg, 'bpfilter') && numel(cfg.bpfilter)==2
  disp 'run bp filter on input data'
  cfg2 = [];
  cfg2.bpfilter = 'yes';
  cfg2.bpfreq = cfg.bpfilter;
  cfg2.padtype = 'mirror';
  cfg2.padding = 4;
  cfg2.bpfiltord = 5; % per JQK recommendation
  data = ft_preprocessing(cfg2, data);
end
if isfield(cfg, 'bsfilter') && numel(cfg.bsfilter)==2
  disp 'run bs filter on input data'
  cfg2 = [];
  cfg2.bsfilter = 'yes';
  cfg2.bsfreq = cfg.bsfilter;
  cfg2.padtype = 'mirror';
  cfg2.padding = 4;
  cfg2.bsfiltord = 5; % per JQK recommendation
  data = ft_preprocessing(cfg2, data);
end
% if ismac
%   cfgfreq              = [];
%   cfgfreq.output       = 'pow';
%   %   cfgfreq.channel      = 'all';
%   cfgfreq.method       = 'mtmfft';
%   cfgfreq.taper        = 'hanning';
%   cfgfreq.keeptrials   = 'no';
%   cfgfreq.foilim       = [0 min(256/2, 200)];
%   cfgfreq.pad='nextpow2';
%   tempfreq = ft_freqanalysis(cfgfreq, data)
%   figure; semilogy(tempfreq.freq, mean(tempfreq.powspctrm))
% end
if isfield(cfg, 'phasetimecourse') && strcmp(cfg.phasetimecourse, 'yes')
  disp 'get phase time course'
  cfg2           = [];
  cfg2.method    = 'mtmconvol'; % mtmfft mtmconvol
  cfg2.output    = 'fourier';
  cfg2.taper     = 'hanning';
  cfg2.tapsmofrq = 2;
  cfg2.pad = 'nextpow2';
  cfg2.channel = 'EEG';  
  cfg2.toi = 'all';
%   cfg2.foi = 4;
  cfg2.foi = 1;
  cfg2.t_ftimwin = 1./cfg2.foi*1; % 3 cycles for each freq
%   cfg2.trials = 1;
%   cfg2.t_ftimwin(1:3) = 1; % fewer cycles for 2 and 3 Hz
  freq          = ft_freqanalysis(cfg2, data);
  
  freq.fourierspctrm = angle(freq.fourierspctrm); % keep the phase
    %convert to ft raw data format
  data = ft_checkdata(freq, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'no');

end

if ismac
  % plot ERP's
  cfg2 = [];
  cfg2.vartrllength = 2;
  cfg2.keeptrials = 'yes';
  timelock = ft_timelockanalysis(cfg2, data);
  
%   cfg2 = [];
%   cfg2.layout = 'elec1010.lay';
%   cfg2.hotkeys = 'yes';
%   
%   cfg2.parameter = 'var' % var avg
%   
%   cfg2.xlim = [-1 1.25];
%   cfg2.zlim = 'maxabs';
%   cfg2.colorbar = 'yes';
%   figure
%   ft_multiplotER(cfg2, timelock)
  
  % browse data
  cfg2 = [];
  cfg2=[];
  cfg2.viewmode = 'vertical';
  cfg2.channel = [12,13,14,16,18,19,34,35,39,43,48];
  figure; ft_databrowser(cfg2, data)
end

if isfield(cfg, 'splithalf') && strcmp(cfg.splithalf, 'yes')
  alltrials = find(cfg.trials);
  [outputpath, outputfile] = fileparts(cfg.outputfile);
  
  cfg.trials = alltrials(1:2:length(alltrials));  
  cfg.outputfile = fullfile(outputpath, 'odd', [outputfile '.mat']);
  mkdir(fullfile(outputpath, 'odd'))
  mse = ft_entropyanalysis(cfg, data);

  cfg.trials = alltrials(2:2:length(alltrials));
  cfg.outputfile = fullfile(outputpath, 'even', [outputfile '.mat']);
  mkdir(fullfile(outputpath, 'even'))
  mse = ft_entropyanalysis(cfg, data);
else
  mse = ft_entropyanalysis(cfg, data);
end



