% GENERATE_BRP_FIG2 reproduces Figure 2 of
%   Schuma, Landau, Diester, Anderson & Karvat: Burst-related potentials as a
%   temporal anchor for cognition, Comm. Biol. 2026
%
% The figure shows a single burst recorded in rat cortex:
%   Fig. 2A top    - time-frequency power around the burst (fold of median)
%   Fig. 2A bottom - raw LFP (black), burst duration (green), and an ideal
%                    cosine at the burst frequency (grey)
%   Fig. 2B top    - spike raster of all well-isolated units
%   Fig. 2B bottom - population firing rate (mean of z-scored rates)
% Red lines/shading mark troughs (+-1/4 cycle); blue lines/shading mark peaks.
%
% USAGE
%   Download the functions and data from https://github.com/laaanchic/BRP into
%   one folder, and enter that folder under BRP_folder below (line 34).
%
% REQUIREMENTS
%   Tested on MATLAB R2020b
%   FieldTrip (tested on version 20241211; https://www.fieldtriptoolbox.org/download/)
%   Image Processing Toolbox (imregionalmax, used in burst detection)
%   Statistics and Machine Learning Toolbox (zscore; prctile in burst detection)
%
% INPUT DATA (BRP_fig2_input.mat)
%   combined - FieldTrip raw structure, whole session as one trial. 
%              Channel 1 is the LFP (uV); channels 2:end are binary spike 
%              trains of sorted units, whose labels end with a quality flag 
%              (0 = noise, excluded here).
%
% Author:  Golan Karvat, 2026

clearvars;

%% Paths and data
BRP_folder = 'enter the folder with the code and data from GitHub';
addpath(genpath(BRP_folder));
if ~exist('ft_defaults', 'file'), error('FieldTrip not found. Add it to the path.'); end

load(fullfile(BRP_folder, 'BRP_fig2_input.mat'), 'combined'); 
fs  = combined.fsample;
lfp = combined.trial{1}(1,:);
spk = combined.trial{1}(2:end,:);
unit_labels = combined.label(2:end);
good_units  = cellfun(@(s) s(end), unit_labels) ~= '0';   % exclude noise clusters
spk = spk(good_units,:);
clear combined unit_labels good_units

%% Population firing rate
width_ms = 25;                       % FWHM of the Gaussian kernel, in ms
gw  = makeGaussWin(fs, width_ms);
FR  = conv2(1, gw, spk, 'same');     % smooth each unit along time
FRz = zscore(FR, [], 2);             % z-score each unit

%% LFP as a FieldTrip structure
raw            = [];
raw.trial{1}   = lfp;
raw.time{1}    = (1:length(lfp))/fs;
raw.fsample    = fs;
raw.label{1}   = 'LFP';
raw.sampleinfo = [1 length(lfp)];

%% Burst detection
pmtr              = [];
pmtr.foi          = 3:30;                         % frequencies of interest (Hz)
pmtr.wave_width   = 5;                            % wavelet width (cycles)
pmtr.min_dur      = 1;                            % minimal burst duration (cycles)
pmtr.min_f_gap    = max(ceil(pmtr.foi./4), 4);    % min. frequency gap between overlapping bursts (Hz)
pmtr.output_power = 1;                            % return the TFR for plotting

[Output, pmtr] = Bursts_detection_BRP(raw, [], pmtr);
bst = Output.burstDists;   % bursts x 18; see Output.col_names for the columns
tfr = Output.pwr;          % TFR, fold of the median power per frequency

%% Select the burst shown in the figure
bi = 278;   % the burst shown in the figure (specific to this dataset)

ntr     = 2.5;                                  % cycles to plot on each side of the trough
ind0    = bst(bi,5);                            % nearest trough (sample)
f       = bst(bi,1);                            % burst frequency (Hz)
cyc     = fs/f;                                 % samples per cycle
inds    = bst(bi,3):bst(bi,4);                  % burst duration (samples)
inds2   = round(ind0-ntr*cyc) : round(ind0+ntr*cyc);   % plotted window (samples)
t       = (inds  - ind0)/fs;                    % burst time axis (s, trough = 0)
t2      = (inds2 - ind0)/fs;                    % window time axis (s, trough = 0)
cyc_num = ceil(-ntr) : floor(ntr);              % cycle numbers relative to the trough
troughs = cyc_num*cyc/fs;                       % trough times (s)
peaks   = troughs + 0.5*cyc/fs;                 % peak times (s)

wv   = cos(2*pi*f*t2) * min(lfp(inds2))/2;      % ideal cosine with a trough at t = 0
spki = spk(:,inds2);                            % spikes in the window
C    = mean(FRz(:,inds2));                      % population firing rate in the window

%% Plotting
figure(1); clf;
set(gcf, 'position', [700 50 530 700]);
rows = 16;

% TFR (Fig. 2A top)
subplot(rows, 1, 1:6); hold on
imagesc(t2, pmtr.foi, tfr(:,inds2));
scatter(0, f, 50, 'w', 'filled');               % burst peak
scatter(0, f, 20, 'k', 'filled');
axis xy
xlim(t2([1 end])); ylim([6 30]);
colormap(flip(ft_colormap('RdBu')));
cb = colorbar; cb.Location = 'northoutside';
set(gca, 'xticklabel', {});

% LFP (Fig. 2A bottom)
subplot(rows, 1, 7:8); hold on
plot(t2, wv*1e-3,          'color', [.5 .5 .5],          'linewidth', .75);   % uV -> mV
plot(t2, lfp(inds2)*1e-3,  'k',                          'linewidth', .75);
plot(t,  lfp(inds)*1e-3,   'color', [4 181 13]/255,      'linewidth', .75);
axis tight
for k = 1:length(troughs), xline(troughs(k), 'r'); xline(peaks(k), 'b'); end   % vector xline needs R2021a+
xlim(t2([1 end]));
xlabel('Seconds');

% Spike raster (Fig. 2B top)
subplot(rows, 1, 11:rows-2); hold on
ys = [1, size(spk,1)];
for k = 1:length(troughs)
    xs = troughs(k) + [-1 1]*0.25*cyc/fs;       % +-1/4 cycle around each trough
    patch(xs([1 2 2 1]),            ys([1 1 2 2]), 'r', 'facealpha', 0.2, 'edgecolor', 'none');
    patch(xs([1 2 2 1]) + cyc/fs/2, ys([1 1 2 2]), 'b', 'facealpha', 0.1, 'edgecolor', 'none');
    if k == 1   % also shade the peak before the first trough
        patch(xs([1 2 2 1]) - cyc/fs/2, ys([1 1 2 2]), 'b', 'facealpha', 0.1, 'edgecolor', 'none');
    end
end
for k = 1:length(troughs), xline(troughs(k), 'r'); xline(peaks(k), 'b'); end   % vector xline needs R2021a+
[unit_i, samp_i] = find(spki);
scatter(t2(samp_i), unit_i, 10, 'ok', 'filled');
xlim(t2([1 end])); ylim(ys);
set(gca, 'xtick', []);

% Population firing rate (Fig. 2B bottom)
subplot(rows, 1, rows + [-1 0]); hold on
plot(t2, C + 0.5, 'k', 'linewidth', .75);       % +0.5: vertical offset for display only
for k = 1:length(troughs), xline(troughs(k), 'r'); xline(peaks(k), 'b'); end   % vector xline needs R2021a+
xlim(t2([1 end]));
set(gca, 'xtick', troughs, 'xticklabel', arrayfun(@num2str, cyc_num, 'UniformOutput', false));
xlabel('Cycles');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
set(findall(gcf, 'type', 'axes'), 'TickDir', 'out');