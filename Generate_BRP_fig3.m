% GENERATE_BRP_FIG3 reproduces Figure 3 of
%   Schuma, Landau, Diester, Anderson & Karvat: Burst-related potentials as a
%   temporal anchor for cognition, Comm. Biol. 2026
%
% EEG data from one participant (channel O2) in a visual task:
%   Fig. 3A - rhythmicity (LAVI) spectrum, defining the frequency bands
%   Fig. 3B - burst-related potentials per band, aligned to the nearest
%             trough, the nearest minimum, and the power maximum (mean +- SD)
%   Fig. 3C - (i) ERP, (ii) baseline-normalized TFR, (iii-iv) BRPs of the
%             first burst after stimulus onset (black) vs. all bursts (colour),
%             with insets showing the latency distribution of the first burst
%             (proportion of all trials)
%
% USAGE
%   Download the functions and data from https://github.com/laaanchic/BRP into
%   one folder, and enter that folder under BRP_folder below (line 31).
%
% REQUIREMENTS
%   Tested on MATLAB R2020b
%   FieldTrip (tested on 20241211; https://www.fieldtriptoolbox.org/download/)
%
% INPUT DATA (BRP_fig3_input.mat)
%   OneChan - FieldTrip raw structure with the O2 channel
%   bst     - burst matrix from Bursts_detection_BRP (see its header for the columns)
%   lavi    - LAVI spectrum of the participant;  foi - its frequencies (Hz)
%   TRL     - trial definition (table)
%
% Author:  Golan Karvat, 2026

clearvars
BRP_folder = 'enter the folder with the code and data from GitHub';
BRP_folder = pwd;
addpath(genpath(BRP_folder));
if ~exist('ft_defaults', 'file'), error('FieldTrip not found. Add it to the path.'); end
ft_notice off; ft_info off
load(fullfile(BRP_folder, 'BRP_fig3_input.mat'), 'OneChan', 'bst', 'lavi', 'foi', 'TRL');
fs = OneChan.fsample;

%% Parameters
boi    = -1:2;                      % bands of interest (band codes from ABBA)
bnames = {'q/a', 'a', 'b1', 'b2'};
LOCKTO = [5, 15, 2];                % bst columns to align to: 5 = nearest trough,
                                    % 15 = nearest minimum, 2 = power maximum
preB   = 3;                         % cycles before the alignment point
postB  = 3;                         % cycles after the alignment point
ylm    = [-25 25];                  % y-limits of the BRP/ERP plots
xlm    = [-0.25 0.75];              % x-limits (s) of the task-related plots
blWin  = [-0.5 0];                  % baseline window (s)
c2     = [76,134,12; 61,121,187]/255;  % band colours, as in the LAVI plot
xlbl   = arrayfun(@num2str, -preB:postB, 'UniformOutput', false);

%% Fig. 3A: LAVI spectrum and frequency bands
abba = ABBA(lavi, foi); abba = abba{1};   % band limits: columns 4-5 (Hz), band code: column 10
figure(1); clf; hold on
abba2 = plotLaviInterp(foi, lavi);        % column 9: band type, used for the colour
xlabel('Frequency (Hz)');
ylabel('Rhythmicity (LAVI)');
set(gca, 'fontsize', 14);

%% Burst minima and TFR
bst = findBurstMinima(OneChan, bst, fs);  % adds the nearest minimum as column 15

cfg        = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foi    = 3:35;
cfg.toi    = 'all';
cfg.width  = 5;
cfg.pad    = 'nextpow2';
TFR = ft_freqanalysis(cfg, OneChan);

%% Fig. 3B: BRPs of all bursts, per band and alignment point
figure(2); clf;
rows = length(LOCKTO);
cols = length(boi);
set(gcf, 'position', [680 18 350*cols 740]);
BDAT    = cell(length(boi), length(LOCKTO));   % burst-locked data, band x alignment
for bi = 1:length(boi)
    band  = abba(abba(:,10) == boi(bi), :);
    f_lim = [floor(band(4)), ceil(band(5))];   % rounded, as bursts were detected in 1 Hz steps
    cyc   = fs/mean(f_lim);                    % samples per cycle
    preBi  = round(cyc*preB);
    postBi = round(cyc*postB);
    bsti  = bst(bst(:,1) >= f_lim(1) & bst(:,1) <= f_lim(2), :);
    c     = c2(2 - (abba2(bi,9) == 1), :);

    for li = 1:length(LOCKTO)
        subplot(rows, cols, (li-1)*cols + bi); cla; hold on
        trlB = [bsti(:,LOCKTO(li)) - preBi, bsti(:,LOCKTO(li)) + postBi, zeros(size(bsti,1),1) - preBi];
        trlB(trlB(:,1) < OneChan.sampleinfo(1) | trlB(:,2) > OneChan.sampleinfo(end), :) = []; % drop bursts at the edges
        cfg     = [];
        cfg.trl = trlB;
        BDAT{bi,li}   = ft_redefinetrial(cfg, OneChan);

        if li == 1   % the same example bursts are highlighted in all rows
            [~, inds] = plotAvgErr_v01(BDAT{bi,li}, c, 1, ylm, 1, 1, 1);
        else
            plotAvgErr_v01(BDAT{bi,li}, c, 1, ylm, 1, 1, 1, inds);
        end
        xtck = (-cyc*preB : cyc : cyc*postB)/fs;
        set(gca, 'xtick', xtck, 'xticklabel', xlbl, 'xlim', xtck([1 end]));
        if li == 1, title(sprintf('%i-%i Hz (%s)', f_lim(1), f_lim(2), bnames{bi})); end
        if li == length(LOCKTO), xlabel('Cycles'); else, set(gca, 'xticklabel', {}); end
        if bi == 1, ylabel('\muV'); else, set(gca, 'yticklabel', {}); end
        set(gca, 'fontsize', 14, 'tickdir', 'out');
    end
end

%% Fig. 3C: task-related activity
figure(3); clf;
set(gcf, 'position', [940 230 642 408]);
lockto = 5;   % align first bursts to their nearest trough

% Trials: 1.5 s epochs, starting 0.25 s after the start sample in TRL, with t = 0 at 0.75 s into the epoch
trl = str2double(table2array(TRL));
trl(:,1) = trl(:,1) + 0.25*fs;
trl(:,2) = trl(:,1) + 1.5*fs;
trl(:,3) = -0.75*fs;
nTrials  = size(trl, 1);

% (i) ERP
cfg     = [];
cfg.trl = trl;
dat0 = ft_redefinetrial(cfg, OneChan);
cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = blWin;
dat0 = ft_preprocessing(cfg, dat0);

subplot(2, 2, 1); hold on
plotAvgErr_v01(dat0, [0 0 0], 1, ylm, 0);
xlim(xlm);
ylabel('\muV');
set(gca, 'fontsize', 14, 'tickdir', 'out');

% (ii) TFR, relative to baseline
cfg     = [];
cfg.fs  = fs;
cfg.trl = trl;
cfg.baselinewindow = blWin;
cfg.blmethod = 'relative';
tfr0 = ft_GK_redefinefreq_v04(cfg, TFR);
pwr  = squeeze(mean(tfr0.powspctrm, 1, 'omitnan'));

subplot(2, 2, 3); hold on
imagesc(tfr0.time, tfr0.freq, pwr);
colormap(flip(ft_colormap('RdBu')));
cb = colorbar('Location', 'north'); cb.AxisLocation = 'out';
cb.TickDirection = 'out'; cb.Ticks = 1:5;
pos = get(gca, 'position'); cpos = cb.Position;
cpos(1) = pos(1); cpos(2) = pos(2) + pos(4)*1.02; cpos(3) = pos(3); cb.Position = cpos;
xlabel('Seconds');
ylabel('Hz');
set(gca, 'xlim', xlm, 'ylim', tfr0.freq([1 end]), 'tickdir', 'out', 'fontsize', 14);

% (iii-iv) BRPs of the first burst after stimulus onset, for the first two bands
for bi = 1:2
    band  = abba(abba(:,10) == boi(bi), :);
    f_lim = [floor(band(4)), ceil(band(5))];
    cyc   = fs/mean(f_lim);
    preBi  = round(cyc*preB);
    postBi = round(cyc*postB);
    bsti  = bst(bst(:,1) >= f_lim(1) & bst(:,1) <= f_lim(2), :);
    c     = c2(2 - (abba2(bi,9) == 1), :);
    [first_ind, first_lat] = findFirstBurst(trl, bsti, fs, lockto);  % sample and latency of each trial's first burst

    subplot(2, 2, 2*bi); cla; hold on
    cfg     = [];
    cfg.trl = [first_ind - preBi, first_ind + postBi, zeros(size(first_ind,1),1) - preBi];
    bdat = ft_redefinetrial(cfg, OneChan);
    plotAvg_line_v01(BDAT{bi,1}, c, 1, ylm);          % all bursts (colour)
    plotAvgErr_2lines_v01(bdat, [0 0 0], 1, ylm);     % first bursts (black)
    set(gca, 'xtick', (-cyc*preB : cyc : cyc*postB)/fs, 'xticklabel', xlbl);
    if bi > 1, xlabel('Cycles'); end
    ylabel('\muV');
    set(gca, 'fontsize', 14, 'tickdir', 'out');

    % inset: latency of the first burst, as a proportion of ALL trials
    % (trials without a burst in the window are counted in the denominator)
    pos = get(gca, 'position');
    ax  = axes;
    ax.Position = [pos(1) + pos(3)*0.65, pos(2) + pos(4)*0.275, pos(3)*0.35, pos(4)*0.2];
    [n, e] = histcounts(first_lat/fs, 0 : diff(xlm)/14 : xlm(end));
    bar(e(1:end-1), n/nTrials, 'facecolor', 'k', 'linewidth', 0.5, 'edgecolor', 'k', 'facealpha', 0.5);
    ylim([0 0.2]);
    set(gca, 'tickdir', 'out');
end