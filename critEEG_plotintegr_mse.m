% do stats in 3D, multiplot, integrated plots

try isstruct(mseavg) % nargin == 0
catch
  mseavg = critEEG_mse_merge();
end
nsub = length(mseavg.SUBJ);
tois = [];
% tois(3,:,:) = [-0.2 0.6; -0.6 0.6]; % lib-cons
tois(3,:,:) = [-0.2 0.5; -0.6 0.6]; % lib-cons
tois(4,:,:) = [-0.2 0.6; -0.6 0.6]; % for cond avg

toitind = {};
controlfordprime = 0;
stattypes = {'ttest', 'corr'};
corrwith = 'crit'; % crit driftbias rpara
clear msestat
for idat = 1%:2 % [1 5] %1:6 %  :2 %1:6 % mse_leg: {'MSEraw'  'MSEpsc'  'r_parameter'  'r_para_psc'  'MSEraw r removed' 'MSEpsc r removed'  'r MSE removed'}
  for icond = 4 %:4 %:4%:4
    for istat = 2 1:2 %:2%:2
      stattype = stattypes{istat}; % lib vs cons ttest, or corr with behavior
      for itrig = 1 %:2 % 1:2
        cfg = [];
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';
        cfg.correctm         = 'cluster';
%             cfg.correctm         = 'no';
        cfg.clusteralpha     = 0.05;
%             cfg.clusteralpha     = 0.1;
        cfg.clusterstatistic = 'maxsum';
        cfg.tail             = 0;
        cfg.clustertail      = 0;
        cfg.alpha            = 0.025;
        % prepare_neighbours determines what sensors may form clusters
        cfg0_neighb = [];
        cfg0_neighb.method    = 'template';
        cfg0_neighb.template  = 'elec1010_neighb.mat';
        cfg0_neighb.elec  = 'standard_1020.elc';
        cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);
        %     cfg.minnbchan = 1;
        
        if strcmp(stattype, 'ttest')
          cfg.numrandomization = 1000;
          design = zeros(2,2*nsub);
          for i = 1:nsub
            design(1,i) = i;
          end
          for i = 1:nsub
            design(1,nsub+i) = i;
          end
          design(2,1:nsub)        = 1;
          design(2,nsub+1:2*nsub) = 2;
          cfg.design   = design;
          cfg.uvar     = 1;
          cfg.ivar     = 2;
        elseif strcmp(stattype, 'corr')
          cfg.numrandomization = 100;
          corrtype = 'Spearman'; % used in paper
          %           corrtype = 'Pearson';
          cfg.statistic        = 'ft_statfun_partialcorrelationT';

          cfg.type = corrtype;
          cfg.uvar     = [];
          cfg.ivar     = 1;
          if strcmp(corrwith, 'crit')
            ib = 2; % 1= dprime, 2 = crit
            cfg.design = mseavg.behavior.(mseavg.behavior.behav_measures{ib})(:,4,icond)';
            % include the other behav var to control for it
            if controlfordprime == 1
              cfg.design = [cfg.design; ...
                mseavg.behavior.(mseavg.behavior.behav_measures{ mod(ib,2)+1 })(:,4,icond)'];
            end
          elseif strcmp(corrwith, 'driftbias')
            driftbias = respavg.behavior.zfix_dcfree.dc;
            cfg.design = driftbias(:,2) - driftbias(:,1);
          elseif strcmp(corrwith, 'rpara') % corr mMSE with r (idat=3)
            cfg.statistic        = 'ft_statfun_correlationT'; % edited to make faster
            cfg.design = reshape(squeeze( mseavg.dat{itrig}(:, 3, :,:,:,  icond, 11) ), 16, []);
          end
        end
        
        cfg.latency = tois(icond,itrig,:); %[-0.2 0.45];
        toitind{icond,itrig} = mseavg.time{itrig} >= cfg.latency(1) & mseavg.time{itrig} <= cfg.latency(2);
        
        freq = [];
        freq.dimord = 'subj_chan_freq_time';
        freq.label = mseavg.label;
        freq.freq = mseavg.timescales;
        freq.time = mseavg.time{itrig};

        freq.powspctrm = squeeze( mseavg.dat{itrig}(:, idat, :,:,:,  icond, 11) ); % all runs
%         freq.powspctrm = squeeze( nanmean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 1:3), 7) ); % ses 1
%         freq.powspctrm = squeeze( nanmean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 4:6), 7) ); % ses 2
%         freq.powspctrm = squeeze( nanmean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 7:9), 7) ); % ses 3

%         subj = ~isnan(freq.powspctrm(:,1,1,1)); % remove empty subj
%         freq.powspctrm = freq.powspctrm(subj,:,:,:);        
%         cfg.design = cfg.design(subj);
        
        freqzero = freq; %create zero freq to test against
        freqzero.powspctrm = zeros(size(freq.powspctrm));
        %         msestat{ itrig, icond, istim, iresp } = ft_freqstatistics(cfg, freq, freqzero);
        
        if strcmp(stattype, 'ttest')
          msestat{ idat, istat, icond, itrig } = ft_freqstatistics(cfg, freq, freqzero);
        elseif strcmp(stattype, 'corr')
          msestat{ idat, istat, icond, itrig } = ft_freqstatistics(cfg, freq);
        end
        %         cfg.avgoverrpt = 'yes';
        temp = ft_selectdata(cfg,freq);
        msestat{ idat, istat, icond, itrig }.powspctrm_subj = temp.powspctrm;
        msestat{ idat, istat, icond, itrig }.powspctrm = squeeze(mean(temp.powspctrm));
        if idat==1 && istat==2 && icond==4 && itrig==1 % mask of the neg critvsmmse corr in main paper
          mseavg.corrmask = msestat{ idat, istat, icond, itrig }.negclusterslabelmat == 1;
          mseavg_root.corrmask = mseavg.corrmask;
        end
        
        clear temp
      end
    end
  end
end

%%  plot 3D integrated cluster NEW STYLE
SAV = 1;
close all
load colormap_jetlightgray.mat
subplotind = [2 1];
clussign = {'pos', 'neg'};
param = {'powspctrm' 'rho'};

cfg=[];
cfg.layout = 'elec1010.lay';
cfg.clus2plot = 1; % to do report pvals
cfg.integratetype = 'trapz'; % mean or trapz
cfg.colormap = cmap;
cfg.subplotsize = [2 2];
f = figure;
f.Position = [ 680          75         350 350 ]; % A4 formaat
irow = 0;
for icond = 4
  for idat = 1 % 1 is raw MSE
    for istat = 2%:2 % 1= main effect, 2 =corr
      cfg.parameter = param{istat};
      if istat==1
%         cfg.titleTFR = 'MMSE Thetaphase lib-cons'; %, 'Fontsize', 12 )
% %         cfg.titleTFR = 'MMSE Theta bp lib-cons'; %, 'Fontsize', 12 )
      else
%         cfg.titleTFR = 'lib-cons crit corr'; %, 'Fontsize', 12 )
      end
      if irow+1 == 1
        cfg.titleTFR = 'lib-cons mMSE'; 
      else
        cfg.titleTFR = 'corr lib-cons mMSE vs crit'; 
      end
      for isign = 1:2 %1:2
        cfg.clussign = clussign{isign};
        for itrig = 1
          cfg.subplotind = irow*numel(subplotind) + subplotind(itrig,:);
          ft_clusterplot3D(cfg,  msestat{ idat, istat, icond, itrig }) % megdatraw.stat{idrug, imod, itrig, ifreq, idiff}
        end
        irow = irow+1;
      end
    end
    if SAV
      saveas(gcf, fullfile(mseavg.PREOUT, sprintf('dat%d_stat%d_cond%d_clus%d%s.pdf',  idat, istat, icond, cfg.clus2plot, clussign{isign})))
      saveas(gcf, fullfile(mseavg.PREOUT, sprintf('dat%d_stat%d_cond%d_clus%d%s.png',  idat, istat, icond, cfg.clus2plot, clussign{isign})))
    end    
  end
end
cd(mseavg.PREOUT)

%% plot integrated clusters: TFR and topo's As in power paper Fig 3
idat = 1; icond = 4; istat = 2; itrig=1; % 1= main effect, 2 =corr
% idat = 1; icond = 4; istat=2; itrig=1; % lib-cons mMSE or neg corr crit vs mMSE

close all
SAV = 1;
load colormap_jetlightgray.mat

trigger_leg = {'stim', 'resp'};
trap = 1;
XTICKS = -2:0.2:2;
YTICKS = {0:20:160, 40:20:150};

nrows = 4;
ncols = 4;
plorder = [2 4 1 4];  %[2 3 1 4]
% plorder = [4 1  ];

clus_leg = {'pos' 'neg'};
THR = 0.99;

% clus2plot{1} = [2 4]; % pos, stim resplocked
clus2plot{1} = [1 1]; % pos, stim resplocked
clus2plot{2} = [1 1]; % neg, stim resplocked

CLIMStfr = [-87 87; -160 160; -10 10; ];
CLIMStopo = [-250 250;    -8 8; -5 5 ];

f = figure;
% set(gcf, 'Position', [100 150 750 125*4]);
set(gcf, 'Position', [100 150 600 100*4]);
set(0, 'defaultaxesfontsize', 6);
irow = 0; hold on
panelind = 'ABCD';
senscols = {'k' 'k'}; % topo circle color
for isign = 2 %1:2 % pos neg
  iplot=0;
  irow = irow + 1; % select row
  for itrig = 1 %1:2%:2 % stim resp
    iclus = clus2plot{isign}(itrig);
    clusfield = [clus_leg{isign} 'clusters'];
    plotstat = msestat{ idat, istat, icond, itrig };
    if ~isfield(plotstat, clusfield) || isempty(plotstat.(clusfield)); disp(clusfield); disp('No clusters found'); continue; end
    
    clus = find([plotstat.(clusfield).prob] < THR);
    cluslabels = plotstat.([clusfield 'labelmat']);
    
    if isempty(clus); continue; end
    chansel = any(any(cluslabels == iclus, 2),3); % get chans before high and low are separated
    
    % plot TFR
    %     dat = squeeze( mean(mseavg.dat{itrig}(:,idat, :,:,toitind{icond, itrig},  icond, 11) ));
    if istat == 1
      dat = plotstat.powspctrm;
    elseif istat == 2
      dat = plotstat.rho;
    end
    if trap
      dat(cluslabels~=clus(iclus)) = 0; % only visualize cluster
      dat = squeeze(trapz(dat(chansel,:,:),1));
    else
      dat(cluslabels~=clus(iclus)) = NaN;
      dat = squeeze(nanmean(dat,1));%./numel(find(dat~=0));
      dat(isnan(dat))=0;
      scale = [-0.05 0.05];
    end
    
    thrcluster = squeeze(any(cluslabels == clus(iclus), 1));
    thrcluster2 = zeros(size(thrcluster)) + 0.1;
    thrcluster2(thrcluster) = 1;
    
    freq = [];
    freq.dat = shiftdim(dat,-1);
    freq.thrcluster2 = shiftdim(thrcluster2,-1);
    freq.time = plotstat.time;
    freq.freq = plotstat.freq;
    freq.dimord = 'chan_freq_time';
    freq.label = {'avg'};
    
    iplot=iplot+1;
    s = subplot(nrows,ncols, plorder(iplot) + (irow*4-4) ); % put on the correct row
    colormap(cmap); hold on
    
    cfg = [];
    cfg.parameter = 'dat';
    cfg.maskparameter = 'thrcluster2';
    cfg.zlim = 'maxabs';
    %       cfg.zlim = CLIMStfr(irow,:);
    cfg.colorbar = 'no';
    cfg.interactive = 'no';
    if itrig == 2
      cfg.title = mseavg.mse_leg{idat};
    else
      cfg.title = mseavg.behav_conds{icond};
    end
    ft_singleplotTFR(cfg, freq);
    hold on
    
    ax = gca;
    ax.Position(3) = (plotstat.time(end) - plotstat.time(1)) * ax.Position(3);
    ax.Box = 'on';
    ax.YTick = YTICKS{1};
    if itrig == 1
      ylabel('Time scale (ms)')
      CLIM = ax.CLim;
    elseif itrig == 2
      %       ax.Position(1) = 0.4475;
      ax.YTickLabel = [];
      %       if strcmp(dattype, 'dat')
      %         c.Label.String = 'Integr. M. (x100%)';
      %       elseif strcmp(dattype, 'pow')
      %         c.Label.String = 'Integrated raw power';
      %       end
    end
    c = colorbar;
    if itrig == 2; c.Label.String = []; end
    c.Position(1) = c.Position(1)+0.065;
    c.Position(2) = c.Position(2)+0.05;
    c.Position(3) = 0.005;
    c.Position(4) = 0.05;
    c.Box = 'off';
    if length(c.Ticks) == 5
      c.Ticks = c.Ticks([1 3 5]);
    end
    
    ax.CLim = CLIM;
    ax.XTick = XTICKS;
    ax.XLim = [plotstat.time(1) plotstat.time(end)];
    ax.YLim = [plotstat.freq(1) plotstat.freq(end)];
    
    plot([0,0], ax.YLim,'k',[0,0], ax.YLim,'k', 'Linewidth', 0.5);
    if strcmp(trigger_leg{itrig}, 'stim')
      stimonset = 0.16;
      plot([stimonset,stimonset], [plotstat.freq(1) plotstat.freq(end)], '--k', 'Linewidth', 0.5 );
      RT = mean(mseavg.behavior.RT(:, 4, 3, 3)); % cond avg RT
      plot(RT, plotstat.freq(end), 'vk', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
%       if irow == 1
%         text(RT, plotstat.freq(end)+10, sprintf('mean\nRT'), 'Fontsize', 6, 'HorizontalAlignment','center');
%       end
    end
%     if irow == 4
      if itrig == 1
        xlabel(sprintf('Time from trial start (s)'), 'HorizontalAlignment','center');
      else
        xlabel(sprintf('Time from response (s)'), 'HorizontalAlignment','center');
      end
%     end
  end
  
  % plot topo
  for itrig = 1%:2
    plotstat = msestat{idat, istat, icond, itrig};
    iclus = clus2plot{isign}(itrig);
    
    if ~isfield(plotstat, clusfield) || isempty(plotstat.(clusfield)); disp(clusfield); disp('No clusters found'); continue; end
    
    clus = find([plotstat.(clusfield).prob] < THR);
    if isempty(clus); continue; end
    
    cluslabels = plotstat.([clusfield 'labelmat']);
    
    iplot=iplot+1;
    subplot(nrows, ncols, plorder(iplot) + (irow*4-4) ); colormap(cmap); hold on
    cfg = [];
    cfg.layout = 'elec1010.lay';
    cfg.comment = 'no';
    cfg.marker = 'off';
    cfg.shading = 'flat';
    cfg.style = 'straight_imsat'; %both  straight
    cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
    cfg.markersize = 3;
    cfg.highlight = 'off';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 20;
    cfg.zlim = 'maxabs'; %  [-150 150]
%     cfg.zlim = [-100 100]; % maxabs [-150 150]
%     cfg.zlim = [-175 175]; % maxabs [-150 150]
%           cfg.zlim = CLIMStopo(irow,:); % ax.CLim;
    
    %                             cfg.highlightchannel = mseavg.label(any(any(plotstat.([clusfield 'labelmat']) == iclus, 2),3));
    cfg.highlightchannel = mseavg.label(any(any(cluslabels == iclus, 2),3));
    cfg.parameter = 'powspctrm';
    cfg.interactive = 'no';
    
    %         dat = squeeze( mean(mseavg.(dattype)(:,1:48, 1:length(mseavg.freq{iband}), tind{itrig}, ...
    %           iband, itrig, 4,icond, istim, iresp) ));
    %     dat = squeeze( mean(mseavg.dat{itrig}(:,4, 2, :,:,:,  icond, istim, iresp) ));
%     dat = squeeze( mean(mseavg.dat{itrig}(:, idat, :,:,toitind{icond,itrig},  icond, 11) ));
    % plot TFR
    %     dat = squeeze( mean(mseavg.dat{itrig}(:,idat, :,:,toitind{icond, itrig},  icond, 11) ));
    if istat == 1
      dat = plotstat.powspctrm;
    elseif istat == 2
      dat = plotstat.rho;
    end
    
    if trap
      dat(cluslabels~=clus(iclus)) = 0;
      if any(isnan(dat))
        warning('nans present!')
        dat(isnan(dat)) = 0;
      end
      dat = squeeze(trapz(dat,2));
      dat = squeeze(trapz(dat,2));
    else
      dat(cluslabels~=clus(iclus)) = NaN;
      dat = squeeze(nanmean(dat,2));
      dat = squeeze(nanmean(dat,2));
      dat(isnan(dat)) = 0;
      %           cfg.zlim = [-0.25 0.25];
    end
    if any(isnan(dat))
      warning('nans present!')
      continue
    end
    
    freq = [];
    freq.label = plotstat.label; % chlabel(LR_subtract_mat(:,2));
    freq.dimord = 'chan';
    freq.powspctrm = dat;
    
    cfgout = ft_topoplotTFR(cfg, freq);
    
    ax=gca;
    if itrig == 1
    else
      ax.Position(1) = 0.65;
    end
    
    %plot sensors with size determined by weight
    weights = double(cluslabels == clus(iclus));
    weights = sum(weights(:,:),2);
    
%     weights = weights / max(weights); % used for plotting
    %         weights = weights / mean(weights(weights > 0)); % weights centered around 1
            weights = weights / sum(weights); % make weights sum to 1
    hold on
    for isens = find(weights > 0)'
      plot(mseavg.sens.pos(isens,1), mseavg.sens.pos(isens,2), 'marker', 'o', 'color', senscols{isign}, 'markersize', 4 * weights(isens), 'linestyle', 'none', 'LineWidth', 0.75)
      %                                 plot(senspos(isens,1), senspos(isens,2), 'marker', 'x', 'color', 'k', 'markersize', 10, 'linestyle', 'none')
    end
    %                             title(sprintf('%s %s p = %g', mseavg.sdt_conds{istim,iresp}, ...
    %                                 mseavg.behav_conds{icond}, plotstat.(clusfield)(iclus).prob))
    %                 title(sprintf('%s p = %g', mseavg.behav_conds{icond}, plotstat.(clusfield)(iclus).prob))
    title(sprintf('p = %1.3f', plotstat.(clusfield)(iclus).prob), 'FontWeight', 'normal')
    
    %         if itrig == 2 || irow == 4
    c = colorbar;
    %     if strcmp(dattype, 'dat')
    %       c.Label.String = 'Integr. M. (x100%)';
    %       if itrig == 1; c.Label.String = []; end
    %     elseif strcmp(dattype, 'pow')
    %       c.Label.String = 'Integrated raw power';
    %     end
    c.Position(1) = c.Position(1)+0.04;
    c.Position(2) = c.Position(2)+0.01;
    c.Position(3) = 0.005;
    c.Position(4) = 0.05;
    c.Box = 'off';
    if length(c.Ticks) == 5
      c.Ticks = c.Ticks([1 3 5]);
    end
    %         end
    if plorder(iplot) == 1
      t = text(-0.6 , 0.7 , panelind(irow), 'Fontsize', 12, 'FontWeight', 'bold');
    end
  end
end


if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('mse_sigclusters_%s_%s_%s_clus%d', mseavg.mse_leg{idat}, mseavg.behav_conds{icond}, stattypes{istat}, clus2plot{isign}(itrig) ));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
  print('-depsc2', outfile)
  print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end
set(0, 'defaultaxesfontsize', 12);

%% plot multiplot for condavg
% close all
SAV = 0;

idat = 4; icond = 3; istat=1; itrig=1; % lib-cons mMSE or neg corr crit vs mMSE
stat = msestat{ idat, istat, icond, itrig };
% stat.mask = stat.([sign 'clusterslabelmat']) == iclus;
% stat.mask = stat.posclusterslabelmat < 57 & stat.posclusterslabelmat > 0;
% stat.mask = stat.negclusterslabelmat == 1;
% stat.mask = stat.posclusterslabelmat == 1;
% stat.mask = stat.negclusterslabelmat == 1 | stat.posclusterslabelmat == 2;
% stat.mask = stat.negclusterslabelmat > 0 | stat.posclusterslabelmat > 0;
%   stat.mask = stat.posclusterslabelmat == 3 | stat.posclusterslabelmat == 1;
% stat.mask = stat.negclusterslabelmat < 4 & stat.negclusterslabelmat > 0;

cfg = [];
if istat == 1;  cfg.parameter = 'powspctrm'; else  cfg.parameter = 'rho'; end
% cfg.maskparameter = 'mask';
% cfg.maskalpha = 0.25;
cfg.layout = 'biosemi64incI1I2.lay';  % biosemi64  elec1010
% cfg.layout = 'biosemi64.lay';  % biosemi64  elec1010
% cfg.layout = 'elec1010.lay';  % biosemi64  elec1010
cfg.layout = ft_prepare_layout(cfg);
cfg.layout.width(:) = 0.075;
cfg.layout.height(:) = 0.075;

% cfg.xlim = [-0.2 0.6];
%                 cfg.zlim = [-0.1 0.1];
cfg.zlim = 'maxabs';
% cfg.zlim = 'maxmin';
% cfg.zlim = [0.2 1.2];
load( 'colormap_jetlightgray.mat')
% cfg.colormap = cmap(129:end, :);
cfg.colormap = cmap;
% cfg.colormap = jet;
% cfg.colormap = parula;
cfg.hotkeys = 'yes';
cfg.fontsize = 12;
cfg.colorbar = 'yes';
cfg.box = 'on';
cfg.showoutline = 'yes';

f = figure;
f.Position = [ 680   678   300   300];

ft_multiplotTFR(cfg, stat);

% % try
% %   title(sprintf('%s %s p = %1.4f', mseavg.mse_leg{idat}, mseavg.behav_conds{icond}, stat.([sign 'clusters'])(iclus).prob))
% % catch
% %   title(sprintf('%s %s', mseavg.mse_leg{idat}, mseavg.behav_conds{icond}))
% % end
% ax=gca;
% ax.FontSize = 8;
% ax.CLim = [0.95 1.25];
% c = colorbar;
% c.Position(1) = 0.9;
% c.Position(2) = 0.5;
% c.Position(3) = c.Position(3) * 0.5;
% c.Position(4) = c.Position(4) * 0.25;
% c.Box = 'off';

if SAV 
  outfile = fullfile(mseavg.PREOUT, sprintf('multiplot_%s_%s', mseavg.mse_leg{idat}, mseavg.behav_conds{icond} ));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
  print('-depsc2', outfile)
  print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end

%% plot scatter significant corr cluster, control for power bands
% % compute traditional bandpower: theta (3-7 Hz), % alpha (8-12 Hz)% beta (13-30 Hz)% gamma (31 - 100 Hz)
% % to control for in the subject MSE measures
% % clus pos: 0 - 0.5 s, lateral frontal chans
% % clus neg: 0.2 - 1 s, occipital chans
% tims = [0 0.5; 0.2 1]; % pos, neg clusters

tims = [-0.2 0.6;-0.2 0.6;-0.2 0.6 ]; % pre- and poststim
% freqbands = [3 7; 8 12; 13 30; 31 100];
freqbands = [1 2; 3 7; 8 12; 13 30; 60 100];
chanroi = [11 10 18]; % ind in sens field
bandpow = nan(length(respavg.SUBJ), 4,length(chanroi)); % SUBJ freqband clus

% chan weights based on MFC mse crit corr
useweights = 1;
if useweights
%   weights = [1,0.144886363636364,0.00852272727272727,0.0880681818181818,0.238636363636364,0.0738636363636364,0,0.204545454545455,0.332386363636364,0.139204545454545,0,0.201704545454545,0,0.0454545454545455,0,0,0.213068181818182,0.392045454545455,0.272727272727273,0.0511363636363636,0.150568181818182,0.267045454545455,0.0369318181818182,0,0,0.0284090909090909,0,0,0,0,0.840909090909091,0.0255681818181818,0.125000000000000,0.241477272727273,0,0.144886363636364,0,0,0.0681818181818182,0,0,0,0.0482954545454545,0.0426136363636364,0.0284090909090909,0.0170454545454545,0.213068181818182,0.323863636363636];

% from r controlled mMSE critcorr chch corr I THINK
%   weights = [0.166430260047281,0.0241134751773050,0.00141843971631206,0.0146572104018913,0.0397163120567376,0.0122931442080378,0,0.0340425531914894,0.0553191489361702,0.0231678486997636,0,0.0335697399527187,0,0.00756501182033097,0,0,0.0354609929078014,0.0652482269503546,0.0453900709219858,0.00851063829787234,0.0250591016548463,0.0444444444444444,0.00614657210401891,0,0,0.00472813238770686,0,0,0,0,0.139952718676123,0.00425531914893617,0.0208037825059102,0.0401891252955083,0,0.0241134751773050,0,0,0.0113475177304965,0,0,0,0.00803782505910165,0.00709219858156028,0.00472813238770686,0.00283687943262411,0.0354609929078014,0.0539007092198582];
  
% % from raw mMSE critcorr chch corr: idat = 1; istat = 2; icond = 4; itrig =
% % 1; CONTROLLED for D'
% weights = [0,0,0,0,0.0256157635467980,0.183251231527094,0,0.00344827586206897,0.119211822660099,0,0,0.00443349753694581,0,0.0147783251231527,0,0,0.00591133004926108,0.0620689655172414,0.0118226600985222,0.00541871921182266,0.0172413793103448,0.112315270935961,0.00689655172413793,0,0,0.0536945812807882,0,0,0,0,0.286206896551724,0.0108374384236453,0.0241379310344828,0,0,0,0,0,0,0,0,0,0,0,0.00197044334975369,0,0.0206896551724138,0.0300492610837438];

% % from raw mMSE critcorr chch corr: idat = 1; istat = 2; icond = 4; itrig =
% % 1; NOT CONTROLLED for D', spearman
% weights = [0,0,0,0.000476190476190476,0.0357142857142857,0.190476190476190,0,0.00285714285714286,0.0595238095238095,0.0157142857142857,0,0,0,0,0,0,0.0138095238095238,0.0542857142857143,0.0180952380952381,0.00619047619047619,0.0542857142857143,0.102857142857143,0.0123809523809524,0.0271428571428571,0.0138095238095238,0.0371428571428571,0,0,0,0,0.234285714285714,0.0109523809523810,0.0409523809523810,0,0,0.000952380952380952,0,0,0.000476190476190476,0,0,0,0,0,0.00380952380952381,0,0.0142857142857143,0.0495238095238095];
% time -0.5 to 1
% weights = [0,0,0,0.000342348510783978,0.0311537144813420,0.159192057514550,0,0.00205409106470387,0.0427935638479973,0.0112975008558713,0,0,0.00616227319411161,0,0,0,0.0116398493666553,0.0578568983224923,0.0215679561793906,0.00445053064019172,0.0725778842862034,0.129407737076344,0.0130092434097912,0.0195138651146868,0.00992810681273536,0.0506675795960288,0,0,0,0,0.200616227319411,0.0157480314960630,0.0438206093803492,0,0,0.000684697021567956,0,0,0.000342348510783978,0,0,0,0,0,0.0356042451215337,0.00513522766175967,0.0126668948990072,0.0417665183156453];

% eLife revision: ERP not removed:
weights = [0,0,0,0.00123001230012300,0.0278802788027880,0.162771627716277,0,0.00492004920049201,0.0742107421074211,0.0155801558015580,0.00246002460024600,0,0,0,0,0,0.0110701107011070,0.0459204592045921,0.00820008200082001,0.0135301353013530,0.0631406314063141,0.151291512915129,0.0225502255022550,0.0196801968019680,0.00861008610086101,0.0992209922099221,0,0,0,0,0.148421484214842,0.00492004920049201,0.0590405904059041,0,0,0.00369003690036900,0,0,0,0,0,0,0,0,0.000820008200082001,0,0.00615006150061501,0.0446904469044690];

% % from r-controlled mMSE critcorr chch corr: idat = 5; istat = 2; icond = 4; itrig =
% % 1; NOT CONTROLLED for D', spearman
% weights = [0,0,0,0.000907441016333938,0.0517241379310345,0.000907441016333938,0.0970961887477314,0.0254083484573503,0.0381125226860254,0.0526315789473684,0,0.0362976406533575,0,0,0,0,0.0317604355716878,0.0980036297640653,0.0499092558983666,0.00816696914700544,0.108892921960073,0.0417422867513612,0.00816696914700544,0.0480943738656987,0,0,0,0,0,0,0.0299455535390200,0.00181488203266788,0.0626134301270417,0.00181488203266788,0,0.00816696914700544,0,0,0.0299455535390200,0,0,0,0,0,0.0190562613430127,0.00453720508166969,0.0317604355716878,0.112522686025408];

% % from raw mMSE critcorr chch corr: idat = 1; istat = 2; icond = 4; itrig =
% % 1; NOT CONTROLLED for D', pearson!
% weights = [0,0,0,0,0.0760697305863708,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00158478605388273,0,0,0.00158478605388273,0.00316957210776545,0,0,0.131537242472266,0,0,0,0,0.749603803486529,0,0.0364500792393027,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

%   weights = weights > 0;
  for iband = 1:5
    for iclus = 1%:length(chanroi)
      if iband == 1 % ultralowfreq
        froi = respavg.frequltralow{4}.freq >= freqbands(iband,1) & respavg.frequltralow{4}.freq <= freqbands(iband,2);
        toi = respavg.frequltralow{4}.time >= tims(iclus,1) & respavg.frequltralow{4}.time <= tims(iclus,2);  
        temp = respavg.frequltralow{4}.powspctrm(:,:,froi,toi);
        temp = mean(temp(:,:,:), 3);
        temp = temp .* weights;
        temp = sum(temp,2);
        bandpow(:,iband, iclus) = squeeze(temp);
        
      elseif iband == 5
        froi = respavg.freq{2} >= freqbands(iband,1) & respavg.freq{2} <= freqbands(iband,2);
        toi = respavg.time{1} >= tims(iclus,1) & respavg.time{1} <= tims(iclus,2);
        temp = respavg.pow(:,1:48,froi,toi,2,1, 4,4,1,1);
        temp = temp(:,:,:);
        temp = mean(temp,3);
        temp = temp .* weights;
        temp = sum(temp,2);
        bandpow(:,iband, iclus) = squeeze(temp);
      else
        froi = respavg.freq{1} >= freqbands(iband,1) & respavg.freq{1} <= freqbands(iband,2);
        toi = respavg.time{1} >= tims(iclus,1) & respavg.time{1} <= tims(iclus,2);
        if iband < 4
          temp = respavg.pow(:,1:48,froi,toi,1,1, 4,4,1,1);
        else
          temp = respavg.pow(:,1:48,froi,toi,2,1, 4,4,1,1);
        end
        temp = temp(:,:,:);
        temp = mean(temp,3);
        temp = temp .* weights;
        temp = sum(temp,2);
        bandpow(:,iband, iclus) = squeeze(temp);
      end
    end
  end
else
  for iband = 1:4
    for iclus = 1:length(chanroi)
      choi = respavg.sens.ind{chanroi(iclus)};
      froi = respavg.freq{1} >= freqbands(iband,1) & respavg.freq{1} <= freqbands(iband,2);
      toi = respavg.time{1} >= tims(iclus,1) & respavg.time{1} <= tims(iclus,2);
      if iband < 4
        temp = respavg.pow(:,choi,froi,toi,1,1, 4,4,1,1);
      else
        temp = respavg.pow(:,choi,froi,toi,2,1, 4,4,1,1);
      end
      temp = mean(temp,2);
      temp = mean(temp,3);
      temp = mean(temp,4);
      bandpow(:,iband, iclus) = squeeze(temp);
    end
  end
end

bandpow = bandpow*1e8;

%% scatter sig cluster mMSE vs crit
close all
SAV=1;
f = figure;
f.Position = [680   678   300   250]; % 450

dotsize = 35;
% idat = 5;  istat = 2; icond = 4; itrig = 1; % neg lib-cons crit corr cluster
idat = 1;  istat = 2; icond = 4; itrig = 1; % neg lib-cons crit corr cluster
% mask =  msestat{ idat, istat, icond, itrig }.mask;
mask  =  msestat{ idat, istat, icond, itrig }.negclusterslabelmat == 1;
msedat = msestat{ idat, istat, icond, itrig }.powspctrm_subj(:,mask);
msedat = mean(msedat,2);
% rdat = squeeze(mseavg.dat{itrig}(:,3, :,:,5:21,  icond, 11)); % 5:21 for -0.2 to 0.6 s
rdat = squeeze(mseavg.dat{itrig}(:,3, :,:,:,  icond, 11)); % 5:21 for -0.2 to 0.6 s
rdat = mean(rdat(:,mask),2);
critdat = mseavg.behavior.criterion(:,4,4);
% critdat = mseavg.behavior.dprime(:,4,4);

% % plot raw mMSE corr with r in cluster
% subplot(3,3,1); hold on; axis square; box on; axis tight
% ax=gca;
% ax.FontSize = 7;
% [r,p] = corr(msedat, rdat, 'type', 'Pearson');
% [rho,prho] = corr(msedat, rdat, 'type', 'Spearman');
% scatter(msedat, rdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% % title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
% title(sprintf('mMSE\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
% xlabel(sprintf('Lib-cons mMSE'));
% ylabel('Lib-cons rdat');
% if p < 0.05; lsline; end;
% ax=gca; ax.XTick = -0.1:0.05:0.1; ax.XLim = [-0.075 0.075];

% plot raw mMSE corr in cluster
subplot(2,3,1); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[r,p] = corr(msedat, critdat, 'type', 'Pearson');
[rho,prho] = corr(msedat, critdat, 'type', 'Spearman');
scatter(msedat, critdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
title(sprintf('mMSE vs crit\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
xlabel(sprintf('Lib-cons mMSE'));
ylabel('Lib-cons crit');
if p < 0.05; lsline; end;
% ax=gca; ax.XTick = -0.1:0.05:0.1; ax.XLim = [-0.075 0.075];


%% plot scatter mMSE vs r for elife revision
close all
SAV=1;
f = figure;
f.Position = [680   678   500   250]; % 450

subplot(2,4,1); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[r,p] = corr(msedat, rdat, 'type', 'Pearson');
[rho,prho] = corr(msedat, rdat, 'type', 'Spearman');
scatter(msedat, rdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
title(sprintf('mMSE vs SD\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
xlabel(sprintf('Lib-cons mMSE'));
ylabel('Lib-cons SD');
if p < 0.05; lsline; end;

% plot scatter delta vs mmse
whichbands = 1;
subplot(2,4,2); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[r,p] = corr(bandpow(:,whichbands,1), msedat, 'type', 'Pearson');
[rho,prho] = corr(bandpow(:,whichbands,1), msedat, 'type', 'Spearman');
scatter(bandpow(:,whichbands,1), msedat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
title(sprintf('delta vs mMSE\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
xlabel(sprintf('Lib-cons delta'));
ylabel('Lib-cons mMSE');
if p < 0.05; lsline; end;

% plot scatter delta vs r
subplot(2,4,3); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[r,p] = corr(bandpow(:,whichbands,1), rdat, 'type', 'Pearson');
[rho,prho] = corr(bandpow(:,whichbands,1), rdat, 'type', 'Spearman');
scatter(bandpow(:,whichbands,1), rdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
title(sprintf('delta vs rdat\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
xlabel(sprintf('Lib-cons delta'));
ylabel('Lib-cons rdat');
if p < 0.05; lsline; end;

% plot scatter crit vs delta
iband=1;
% plot scatter crit vs theta
% iband=2;
subplot(2,4,4); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[r,p] = corr( bandpow(:,iband,1), critdat, 'type', 'Pearson');
[rho,prho] = corr( bandpow(:,iband,1), critdat, 'type', 'Spearman');
scatter( bandpow(:,iband,1), critdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize );
% title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
title(sprintf('delta vs crit\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
xlabel('Lib-cons delta');
ylabel(sprintf('Lib-cons crit'));
if p < 0.05; lsline; end;

% plot scatter crit vs SD
subplot(2,4,5); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[r,p] = corr(critdat, rdat, 'type', 'Pearson');
[rho,prho] = corr(critdat, rdat, 'type', 'Spearman');
scatter(rdat, critdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize );
% title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
title(sprintf('SD vs criterion\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
xlabel(sprintf('Lib-cons SD'));
ylabel('Lib-cons crit');
if p < 0.05; lsline; end;

% plot corr controlled for SD
subplot(2,4,6); hold on; axis square; box on; axis tight
ax=gca;
ax.FontSize = 7;
[~,~,residualmse] = regress(msedat, [rdat ones(16,1)] );
[~,~,residualcrit] = regress(critdat, [rdat  ones(16,1)] );
% [r,p] = corr(residualmse, residualcrit, 'type', corrtype);
[r,p] = partialcorr([msedat critdat], rdat, 'type', 'Pearson'); r = r(2); p = p(2);
[rho,prho] = partialcorr([msedat critdat], rdat, 'type', 'Spearman'); rho = rho(2); prho = prho(2);
scatter(residualmse, residualcrit, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
title(sprintf('controlled for SD\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
% if p < 0.05;    lsline;   end
xlabel(sprintf('Lib-cons mMSE'));
ylabel('Lib-cons criterion');
ax.XLim = [-0.04 0.035];
if p<0.05; lsline; end;

% plot scatter controlled for major pow bands
whichbands = 1 %1:5;
dotsize = 35
subplot(2,4,7);
hold on; axis square; box on; axis tight
ax=gca; ax.FontSize = 7;
[~,~,residualmse2] = regress(msedat, [bandpow(:,whichbands,1) ones(16,1)] );
[~,~,residualcrit2] = regress(critdat, [bandpow(:,whichbands,1) ones(16,1)] );
% [r,p] = corr(residualmse2, residualcrit2, 'type', corrtype);
[r,p] = partialcorr([msedat critdat],  bandpow(:,whichbands,1), 'type', 'Pearson'); r = r(2); p = p(2);
[rho,prho] = partialcorr([msedat critdat],  bandpow(:,whichbands,1), 'type', 'Spearman'); rho = rho(2); prho = prho(2);
scatter(residualmse2, residualcrit2, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('controlled for powbands\nrho = %1.2f\np = %1.3f', r(2), p(2)))
% title(sprintf('controlled for powbands\nrho = %1.2f\np = %1.3f', r(1), p(1)))
title(sprintf('controlled for all bands\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
% if p < 0.05;    lsline;   end
xlabel(sprintf('Lib-cons mMSE'));
ylabel('Lib-cons criterion');
ax.XLim = [-0.04 0.035];
if p<0.05; lsline; end;

% plot scatter rdat vs crit controlled for mmse
subplot(2,4,8); hold on; axis square; box on; %axis tight
ax=gca; ax.FontSize = 7;
[~,~,residualmse2] = regress(rdat, [msedat ones(16,1)] );
[~,~,residualcrit2] = regress(critdat, [msedat ones(16,1)] );
% [r,p] = corr(residualmse2, residualcrit2, 'type', corrtype);
[r,p] = partialcorr([rdat critdat], msedat, 'type', 'Pearson'); r = r(2); p = p(2);
[rho,prho] = partialcorr([rdat critdat],  msedat, 'type', 'Spearman'); rho = rho(2); prho = prho(2);
scatter(residualmse2, residualcrit2, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('controlled for powbands\nrho = %1.2f\np = %1.3f', r(2), p(2)))
% title(sprintf('controlled for powbands\nrho = %1.2f\np = %1.3f', r(1), p(1)))
title(sprintf('controlled for mMSE\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
% if p < 0.05;    lsline;   end
xlabel(sprintf('Lib-cons SD'));
ylabel('Lib-cons criterion');
% ax.XLim = [-0.04 0.035];
if p<0.05; lsline; end;

% % plot raw mMSE corr in cluster controlled for d'
% dprime = mseavg.behavior.dprime(:,4,icond);
% subplot(3,3,2); hold on; axis square; box on; %axis tight
% ax=gca;
% ax.FontSize = 7;
% [~,~,residualmse] = regress(msedat, [dprime ones(16,1)] );
% [~,~,residualcrit] = regress(critdat, [dprime  ones(16,1)] );
% % [r,p] = corr(residualmse, residualcrit, 'type', corrtype);
% [r,p] = partialcorr([msedat critdat], dprime, 'type', 'Pearson'); r = r(2); p = p(2);
% [rho,prho] = partialcorr([msedat critdat], dprime, 'type', 'Spearman'); rho = rho(2); prho = prho(2);
% scatter(residualmse, residualcrit, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% title(sprintf('controlled for dprime\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
% % if p < 0.05;    lsline;   end
% xlabel(sprintf('Lib-cons mMSE'));
% ylabel('Lib-cons criterion');
% ax.XLim = [-0.04 0.035];
% if p<0.05; lsline; end;
% 
% % plot scatter delta vs crit
% subplot(3,3,5); hold on; axis square; box on; %axis tight
% ax=gca;
% ax.FontSize = 7;
% [r,p] = corr(bandpow(:,whichbands,1), critdat, 'type', 'Pearson');
% [rho,prho] = corr(bandpow(:,whichbands,1), critdat, 'type', 'Spearman');
% scatter(bandpow(:,whichbands,1), critdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
% % title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
% title(sprintf('delta vs crit\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', r, p, rho, prho))
% xlabel(sprintf('Lib-cons delta'));
% ylabel('Lib-cons criterion');
% if p < 0.05; lsline; end;

if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('corr_MSEvs_crit_controls'));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
    print('-depsc2', outfile)
    print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end

%% plot corr separately for lib and cons, per revision 2 request
close all
SAV=1;
f = figure;
f.Position = [680   678   500   250]; % 450

dotsize = 45;
idat = 1;  istat = 2; icond = 4; itrig = 1; % neg lib-cons crit corr cluster
mask  =  msestat{ idat, istat, icond, itrig }.negclusterslabelmat == 1;

msecond = []; behav=[];
for icond = 1:2
  msecond(icond).label = mseavg.label;
  msecond(icond).time = mseavg.time{itrig};
  msecond(icond).freq = mseavg.timescales;
  msecond(icond).powspctrm = squeeze( mseavg.dat{itrig}(:, idat, :,:,:,  icond, 11) ); % all runs
  msecond(icond).dimord = 'rpt_chan_freq_time';
  behav(:,icond) = squeeze(mseavg.behavior.criterion(:,4,icond));
end
cfg = [];
cfg.latency = [-0.2 0.6];
msecond = arrayfun(@(x) ft_selectdata(cfg, x), msecond);
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
msecond(4) = ft_math(cfg, msecond(1), msecond(2));

% plotting
for icond= 1:2 [1,2,4]
  msedat = mean(msecond(icond).powspctrm(:,mask),2);
  behavdat = behav(:,icond);
  
  subplot(2,4,icond); hold on; axis square; box on; axis tight
  ax=gca;
  ax.FontSize = 8;
  [r,p] = corr(msedat, behavdat, 'type', 'Pearson');
  [rho,prho] = corr(msedat, behavdat, 'type', 'Spearman');
  scatter(msedat, behavdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
  % title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
  title(sprintf('%s mMSE vs crit\nr = %1.2f,p = %1.3f\nrho = %1.2f,p = %1.3f', mseavg.behav_conds{icond}, r, p, rho, prho))
  xlabel(sprintf('mMSE'));
  ylabel('criterion');
  if p < 0.05; lsline; end;
end

msedat_cons = mean(msecond(1).powspctrm(:,mask),2);
msedat_lib = mean(msecond(2).powspctrm(:,mask),2);
deltarho = randtest_corr([msedat_cons behav(:,1)], [msedat_lib behav(:,2)] ,0,10000, 'Spearman')
deltar = randtest_corr([msedat_cons behav(:,1)], [msedat_lib behav(:,2)] ,0,10000, 'Pearson')

if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('corr_MSEsinglecondvs_crit'));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
    print('-depsc2', outfile)
    print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end
%% plot mMSE separately for lib and cons in corr mask, per revision 2 request
% integrate over chans, timescales and time
close all
SAV=1;

idat = 1;  istat = 2; icond = 4; itrig = 1; % neg lib-cons crit corr cluster
mask  =  msestat{ idat, istat, icond, itrig }.negclusterslabelmat == 1;

msecond = []; behav=[];
for icond = 1:2
  msecond(icond).label = mseavg.label;
  msecond(icond).time = mseavg.time{itrig};
  msecond(icond).freq = mseavg.timescales;
  msecond(icond).powspctrm = squeeze( mean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 11))); % all runs
  msecond(icond).dimord = 'chan_freq_time';
end
cfg = [];
cfg.latency = [-0.2 0.6];
msecond = arrayfun(@(x) ft_selectdata(cfg, x), msecond);

msecond(3) = ft_freqgrandaverage([],  msecond(2), msecond(1));
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
msecond(4) = ft_math(cfg, msecond(2), msecond(1));

% plotting 
close all
cfg=[];
cfg.layout = 'elec1010.lay';
cfg.clussign = 'neg';
cfg.clus2plot = 1; % to do report pvals
cfg.integratetype = 'mean'; % mean or trapz
cfg.subplotsize = [4 2];
cfg.parameter ='powspctrm';

f = figure;
f.Position = [ 680          75         350 700 ]; % A4 formaat
irow = 0;
for icond = [4 4] 1:2 %[1,2,4]; % 1:3 %
  if icond < 4
    cfg.colormap = cmap(129:end,:); % cmap; %
%     cfg.zlim = 'maxabs' %; [1.05 1.15]; %
%     cfg.zlim = [1 1.2]; %
    cfg.zlimTFR = [1 1.2]; %
    cfg.zlimtopo = [0.1 0.5]; %
  else
    cfg.colormap = cmap;
    %     cfg.zlim = 'maxabs' %;
%     cfg.zlim = [-0.1 0.1];
    cfg.zlimTFR = [-0.1 0.1]; %
    cfg.zlimtopo = [-0.025 0.025]; %
  end
  cfg.subplotind = irow*numel(subplotind) + subplotind(itrig,:);
  stat = msestat{ 1, 2, 4, 1 }; % idat, istat, icond, itrig % take corr stat, but plot powspctrm
  stat.mask = mask; % always give clusterplot3D a mask
  stat.powspctrm = squeeze(msecond(icond).powspctrm);
  %   stat.powspctrm = squeeze(mean(msecond(2).powspctrm)) - squeeze(mean(msecond(1).powspctrm));
%   cfg.maskparameter = 'mask';
  ft_clusterplot3D(cfg,  stat) % megdatraw.stat{idrug, imod, itrig, ifreq, idiff}
  irow = irow+1;
end
  if SAV
    saveas(gcf, fullfile(mseavg.PREOUT, sprintf('dat%d_stat%d_cond%d_clus%d%s.pdf',  idat, istat, icond, cfg.clus2plot, clussign{isign})))
    saveas(gcf, fullfile(mseavg.PREOUT, sprintf('dat%d_stat%d_cond%d_clus%d%s.png',  idat, istat, icond, cfg.clus2plot, clussign{isign})))
  end
cd(mseavg.PREOUT)

%% plot time courses lib and cons
close all
linecol = {'r' 'b'};
f = figure;
f.Position = [ 680          75         350 700 ]; % A4 formaat
irow = 0;

for icond = 1:2
  subplot(4,2,1); hold on
  MSE = squeeze( mseavg.dat{itrig}(:, idat, :,:,:,  icond, 11)); % all runs
%   MSE(~mask) = 0;
%   MSE = squeeze(trapz(MSE));
%   MSE = squeeze(trapz(MSE));
  MSE(:, ~mask) = NaN;
  MSE = squeeze(nanmean(MSE,2));
  MSE = squeeze(nanmean(MSE,2)); 
  
  shadedErrorBar(mseavg.time{itrig}, mean(MSE), std(MSE)/sqrt(size(MSE,1)), linecol{icond}, 1)
  xlim([-0.2 0.6])
  ylim([1.09 1.14])
end
for icond = 1:2
  subplot(4,2,2); hold on
  MSE = squeeze( mseavg.dat{itrig}(:, idat, :,:,:,  icond, 11)); % all runs
  MSE(~mask) = 0;
  MSE = squeeze(trapz(MSE,2));
  MSE = squeeze(trapz(MSE,2));
%   MSE(:, ~mask) = NaN;
%   MSE = squeeze(nanmean(MSE,2));
%   MSE = squeeze(nanmean(MSE,2)); 
  
  shadedErrorBar(mseavg.time{itrig}, mean(MSE), std(MSE)/sqrt(size(MSE,1)), linecol{icond}, 1)
  xlim([-0.2 1])
  ylim([2100 2200])
end
  




%% venn diagram with delta, mMSE and crit
close all
deltadat = bandpow(:,1,1);
corrtype = 'Spearman';
% corrtype = 'Pearson';
SAV=1;
f = figure;
f.Position = [680   678   300   300];
r = [];
r(1) = corr(msedat, critdat, 'type', corrtype);
r(2) = corr(msedat, deltadat, 'type', corrtype);
r(3) = corr(critdat, deltadat, 'type', corrtype);
r(4) = r(1) - r(2); % TODO overlap of the 3
% regress
r = r.^2;

% cov([msedat, critdat, deltadat])
% r= corrcoef([msedat, critdat, deltadat])
% r = r.^2;

areas = [];
areas(1) = std(msedat);
areas(2) = std(critdat);
areas(3) = std(deltadat);

% test = venn(areas, r); % mMSE, crit, delta
test = venn([2 2 2], r); % mMSE, crit, delta

legend({ 'mMSE', 'crit', 'delta'})

% test = venn([2 2 2], [0.5 0.75 0.5 0.01 ]); % mMSE, crit, delta

%% Scatters of mMSE vs 5 pow bands
close all
SAV=1;
f = figure;
f.Position = [680   678   300   300];
corrtype = 'Spearman';
% corrtype = 'Pearson';

% idat = 5;  istat = 2; icond = 4; itrig = 1; % neg lib-cons crit corr cluster
idat = 1;  istat = 2; icond = 4; itrig = 1; % neg lib-cons crit corr cluster
% mask =  msestat{ idat, istat, icond, itrig }.mask;
mask  =  msestat{ idat, istat, icond, itrig }.negclusterslabelmat == 1;
msedat = msestat{ idat, istat, icond, itrig }.powspctrm_subj(:,mask);
msedat = mean(msedat,2);
rdat = squeeze(mseavg.dat{itrig}(:,3, :,:,5:21,  icond, 11)); % 5:21 for -0.2 to 0.6 s
rdat = mean(rdat(:,mask),2);
critdat = mseavg.behavior.criterion(:,4,4);

% plot raw mMSE corr in cluster
for iband = 1:5
  subplot(2,3,iband); hold on; axis square; box on; %axis tight
  ax=gca;
  ax.FontSize = 7;
  [r,p] = corr(msedat, bandpow(:,iband,1), 'type', corrtype);
  scatter(msedat,  bandpow(:,iband,1), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', dotsize);
  title(sprintf('raw mMSE\nrho = %1.2f\np = %1.3f', r, p))
  xlabel(sprintf('Lib-cons mMSE'));
  ylabel(sprintf('Lib-cons power band %d', iband));
  if p<0.05; lsline; end;
end
if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('corr_MSEvs_5bands'));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
    print('-depsc2', outfile)
    print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end



%% export miniblock data for rm correlation Lib - cons
close all
SAV = 1;
idat = 1; icond = 4; itrig = 1;
acrossses = 0; % 0 = across runs

freq = []; critdat =[];
freq.dimord = 'subj_chan_freq_time';
freq.label = mseavg.label;
freq.freq = mseavg.timescales;
freq.time = mseavg.time{itrig};
if acrossses
  freq.powspctrm(:,:,:,:,1) = squeeze( nanmean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 1:3), 7 ));
  freq.powspctrm(:,:,:,:,2) = squeeze( nanmean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 4:6), 7 ));
  freq.powspctrm(:,:,:,:,3) = squeeze( nanmean(mseavg.dat{itrig}(:, idat, :,:,:,  icond, 7:9), 7 ));
  critdat(:,1) = squeeze( nanmean(mseavg.behavperrun.criterion(:,4,1:3), 3));
  critdat(:,2) = squeeze( nanmean(mseavg.behavperrun.criterion(:,4,4:6), 3));
  critdat(:,3) = squeeze( nanmean(mseavg.behavperrun.criterion(:,4,7:9), 3));
else
  freq.powspctrm = squeeze( mseavg.dat{itrig}(:, idat, :,:,:,  icond, 1:9) );
  critdat = squeeze(mseavg.behavperrun.criterion(:,4,1:9));
end

cfg = [];
cfg.latency = [-0.2 0.6]; %[0 0.8];
% cfg.latency = [0 0.6]; %[0 0.8];
freq = ft_selectdata(cfg, freq);

msedat = permute(freq.powspctrm, [1 5 2 3 4]);
msedat = msedat(:,:,mask);
msedat = nanmean(msedat,3);

% plot single subj
f = figure;
f.Position = [   680   303   900   800];

incsubj = [1:11, 13:16];
% incsubj = 1:16;

pcoeff= []; rho = [];
for isub = incsubj %1:size(msedat,1)
  subplot(4,4,isub); hold on; axis square; box on
  scatter(msedat(isub,:), critdat(isub,:));
  pcoeffdat = [msedat(isub,:)', critdat(isub,:)'];
  pcoeffdat = pcoeffdat(~isnan(pcoeffdat(:,1)),:);
  pcoeffdat = pcoeffdat(~isnan(pcoeffdat(:,2)),:);
  pcoeff(isub,:) = polyfit(pcoeffdat(:,1), pcoeffdat(:,2), 1);
  lsline
  rho(isub,:) = corr(pcoeffdat(:,1), pcoeffdat(:,2), 'type', 'Pearson');
end
rho = rho(incsubj);

% critdat(critdat > 0.03) = NaN
% validsubj = all(~isnan(msedat),2);
% msedat = msedat(validsubj,:);
% critdat = critdat(validsubj,:);

% fixed effects scatter
f = figure;
f.Position = [   680   303   900   800];
scatter(msedat(:), critdat(:));
lsline

% demean each subj, do fixed efffects scatter (same as rmcorr)
msedat2 = msedat - nanmean(msedat,2);
critdat2 = critdat - nanmean(critdat,2);
corrdat = [msedat2(:) critdat2(:)];
corrdat = corrdat(~isnan(corrdat(:,1)),:);
corrdat = corrdat(~isnan(corrdat(:,2)),:);

f = figure;
f.Position = [   680   303   450   400];
scatter(corrdat(:,1), corrdat(:,2));
rho = corr(corrdat(:,1), corrdat(:,2), 'type', 'Pearson');
lsline
title(sprintf('rho = %g', rho))
xlabel('lib-cons mMSE'); xlabel('lib-cons crit')

% lslines
cmap = jet(16);
pcoeff(:,2) = 0; % don't care about intercept
pcoeff = pcoeff(incsubj,:); % drop 12 (Rinske), only 5 observations
f = figure; hold on; axis square; box on
f.Position = [   680   303   150   125];
for isub = 1:size(pcoeff,1)
  lsvals = polyval(pcoeff(isub,:), [-0.1 0.1]);
  plot([-1 1], lsvals, 'Color', cmap(isub,:))
end
lsvals = polyval(mean(pcoeff), [-0.1 0.1]);
plot([-1 1], lsvals, 'k', 'Linewidth', 3)
[~,tpval]=ttest(pcoeff(:,1));
title(sprintf('%d negative, %d positive, p = %g', sum(pcoeff(:,1) < 0), sum(pcoeff(:,1) > 0), tpval ))
ylim([-1.2 1.2])

if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('lslines_corr_MSEvs_crit'));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
  print('-depsc2', outfile)
  print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end
ax=gca;
ax.FontSize = 7;

% export_for_rmcorr(msedat, critdat, 'mse', 'crit', mseavg.PREOUT);
export_for_rmcorr(msedat(incsubj,:), critdat(incsubj,:), 'mse', 'crit', mseavg.PREOUT);

%% correlate MSE psc vs r psc, use MSE psc as mask
istat=1; icond=3; itrig =1; %1:6 % mse_leg: {'MSEraw'  'MSEpsc'  'r_parameter'  'r_para_psc'  'MSEraw r removed' 'MSEpsc r removed'}
corrtype = 'Spearman'; % Spearman Pearson
close all
f = figure;
f.Position =[   680   678   400   400];
SAV = 1;
ctr=0;
maskleg = {'mse mask' 'r para mask'};
imask = 2; %[2,4]% use mse dat and r dat as mask
dattypes = {'raw' 'psc'};
signtypes = {'pos', 'neg'};
for idat = 2:-1:1%[2,4]% use mse dat and r dat as mask
  for isign = 1:2
    ctr = ctr+1;
    
    %   mse = msestat{ 2, istat, icond, itrig };
    %   rpara = msestat{ 4, istat, icond, itrig };
    mse = msestat{ idat, istat, icond, itrig };
    rpara = msestat{ idat+2, istat, icond, itrig };
    
    % neg cluster
    mask = msestat{ imask, istat, icond, itrig }.([signtypes{isign} 'clusterslabelmat']) == 1;
    msedat = mse.powspctrm(:,mask);
    msedat = mean(msedat,2);
    
    rdat = rpara.powspctrm(:,mask);
    rdat = mean(rdat,2);
    
    subplot(2,2,ctr); hold on; axis square; box on; axis tight
    scatter(msedat, rdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 100 );
    [r,p] = corr(msedat, rdat, 'type', corrtype);
    title(sprintf('%s %s cluster:\nr = %1.2f, p = %1.4f', maskleg{imask/2}, signtypes{isign}, r, p))
    if p < 0.05;    lsline;   end
    xlabel(sprintf('MSE (%s)', dattypes{idat}));
    ylabel(sprintf('r para (%s)', dattypes{idat}));
    
  end
  %   ctr = ctr+1;
  %   % pos cluster
  %   mask = msestat{ imask, istat, icond, itrig }.posclusterslabelmat == 1;
  %   msedat = mse.powspctrm(:,mask);
  %   msedat = mean(msedat,2);
  %
  %   rdat = rpara.powspctrm(:,mask);
  %   rdat = mean(rdat,2);
  %
  %   subplot(2,2,ctr); hold on; axis square; box on; axis tight
  %   scatter(msedat, rdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 100 );
  %   [r,p] = corr(msedat, rdat, 'type', corrtype);
  %   title(sprintf('%s pos cluster:\nr = %1.2f, p = %1.4f', maskleg{imask/2}, r, p))
  %   xlabel(sprintf('MSE (%s)', dattypes{idat}));
  %   ylabel(sprintf('r para (%s)', dattypes{idat}));
end

if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('corr_MSEvs_r'));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
  %   print('-depsc2', outfile)
  %   print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end





%% correlate MSE psc vs power psc, use MSE psc as mask
istat=1; icond=3; itrig =1; %1:6 % mse_leg: {'MSEraw'  'MSEpsc'  'r_parameter'  'r_para_psc'  'MSEraw r removed' 'MSEpsc r removed'}
corrtype = 'Spearman'; % Spearman Pearson
close all
f = figure;
f.Position =[   680   678   400   400];
SAV = 1;
ctr=0;
maskleg = {'mse mask' 'r para mask'};
imask = 2; %[2,4]% use mse dat and r dat as mask
dattypes = {'raw' 'psc'};
signtypes = {'pos', 'neg'};
freqoi = {[3 7] [3 7]; [12 30] [12 30]};
% freqoi = {[3 7] [3 7]; [8 12] [8 12]};
soi = [18 1];
for idat = 2% :-1:1%[2,4]% use mse dat and r dat as mask
  for isign =  1:2
    for ifreq =  1:2
      ctr = ctr+1;
      
      %   mse = msestat{ 2, istat, icond, itrig };
      %   rpara = msestat{ 4, istat, icond, itrig };
      mse = msestat{ idat, istat, icond, itrig };
      %     rpara = msestat{ idat+2, istat, icond, itrig };
      
      % get alpha power psc
      freq = [];
      freq.label = respavg.label(1:48);
      freq.freq = respavg.freq{1};
      freq.time = respavg.time{itrig};
      freq.dimord = 'subj_chan_freq_time';
      freq.powspctrm = respavg.dat(:,1:48,:,:,1,itrig, 4,icond,3,3);
      cfg = [];
      cfg.channel = respavg.sens.ind{soi(isign)};
      cfg.frequency = freqoi{ifreq, isign};
      cfg.latency = [0.1 0.5];
      cfg.avgovertime = 'yes';
      cfg.avgoverfreq = 'yes';
      cfg.avgoverchan = 'yes';
      freq = ft_selectdata(cfg, freq)
      
      % neg cluster
      mask = msestat{ imask, istat, icond, itrig }.([signtypes{isign} 'clusterslabelmat']) == 1;
      msedat = mse.powspctrm_subj(:,mask);
      msedat = mean(msedat,2);
      
      
      
      subplot(2,2,ctr); hold on; axis square; box on; %axis tight
      scatter(msedat, freq.powspctrm, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 100 );
      [r,p] = corr(msedat, freq.powspctrm, 'type', corrtype);
      title(sprintf('%s %s cluster:\nr = %1.2f, p = %1.4f', maskleg{imask/2}, signtypes{isign}, r, p))
      %     if p < 0.05;    lsline;   end
      xlabel(sprintf('%s MSE (%s)', respavg.sens.leg{soi(isign)}, dattypes{idat}));
      ylabel(sprintf('%d-%d Hz power (%s)', cfg.frequency, dattypes{idat}));
      ax=gca;
      ax.XTick = -5:2.5:5;
      if ctr > 2
        ax.XLim = [-5 1];
      end
    end
  end
  %   ctr = ctr+1;
  %   % pos cluster
  %   mask = msestat{ imask, istat, icond, itrig }.posclusterslabelmat == 1;
  %   msedat = mse.powspctrm(:,mask);
  %   msedat = mean(msedat,2);
  %
  %   rdat = rpara.powspctrm(:,mask);
  %   rdat = mean(rdat,2);
  %
  %   subplot(2,2,ctr); hold on; axis square; box on; axis tight
  %   scatter(msedat, rdat, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 100 );
  %   [r,p] = corr(msedat, rdat, 'type', corrtype);
  %   title(sprintf('%s pos cluster:\nr = %1.2f, p = %1.4f', maskleg{imask/2}, r, p))
  %   xlabel(sprintf('MSE (%s)', dattypes{idat}));
  %   ylabel(sprintf('r para (%s)', dattypes{idat}));
end

if SAV %&& ~isempty(clus) && pvalues(iclus) < 1
  outfile = fullfile(mseavg.PREOUT, sprintf('corr_MSEvs_pow_%scluster', signtypes{isign}));
  disp(outfile)
  mkdir(fullfile(mseavg.PREOUT))
  f.Renderer = 'Painters'; % Painters  OpenGL
  fprintf(outfile)
  %   print('-depsc2', outfile)
    print('-dpdf', outfile)
  print('-dpng', outfile)
  %                                         export_fig( outfile, '-deps')
  %                             export_fig( outfile, '-png')
  cd(fullfile(mseavg.PREOUT))
end

%% correlate posterior quenching with dprime
itrig = 1; idat = 2; %MSEpsc 
icond = 4;

freq = [];
freq.dimord = 'subj_chan_freq_time';
freq.label = mseavg.label;
freq.freq = mseavg.timescales;
freq.time = mseavg.time{itrig};
%         freq.powspctrm = squeeze( mseavg.dat{itrig}(:,4, 2, :,:,:,  icond, istim, iresp) );
freq.powspctrm = squeeze( mseavg.dat{itrig}(:, idat, :,:,:,  icond, 11) );

sens = respavg.sens.ind{1};
scoi = 1:10; % up to 40 ms 
tim = respavg.time{1} >= 0 & respavg.time{1} <= 0.6;

cfg=[];
% quenching
cfg.latency = [0 0.6];
cfg.channel = respavg.sens.ind{1};
cfg.frequency = [3 40];
cfg.avgoverfreq = 'yes'; cfg.avgovertime = 'yes'; cfg.avgoverchan = 'yes';
freqsel = ft_selectdata(cfg, freq);

behav = mseavg.behavior.dprime(:,4,icond);
% behav = mseavg.behavior.criterion(:,4,icond);

% corr(freqsel.powspctrm, dprime, 'type', 'Pearson')
[r,p] = corr(freqsel.powspctrm, behav, 'type', 'Spearman')
figure; scatter(freqsel.powspctrm, behav)





