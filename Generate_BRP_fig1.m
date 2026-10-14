%% Generate_BRP_fig1.m
% Simulation code for the schematic figure (Fig. 1) in:
% Schuma, Landau, Diester, Anderson, & Karvat: Burst-related potentials as 
% a temporal anchor for cognition, Comm. Biol. 2026
%
% Simulates T trials of aperiodic (1/f) background activity, each containing
% one oscillatory burst whose onset is jittered relative to a "stimulus" at
% t = 0. The same data are then averaged in two ways:
%   left column:  stimulus-locked (ERP)
%   right column: burst-locked (BRP), each trial aligned to its deepest trough
% Top panels show the first ttp (trials to plot) trials concatenated as 
% continuous data, with bursts shaded and alignment events marked. Bottom 
% panels show these trials individually (grey) and the mean +/- SEM across 
% all T trials.
%
% The MATLAB output was subsequently edited in Adobe Illustrator for the
% published figure.
%
% Requirements:
%   - MATLAB R2018b or later 
%   - Signal Processing Toolbox
%
% rng(42) reproduces the published figure. Changing parameters that alter
% the number of random draws (e.g. T, or the trial window) changes all
% subsequent random values and therefore the figure.
%
% Author:  Golan Karvat, 2026

clearvars;
rng(42);

%% Parameters
T            = 10;          % total number of trials
ttp          = 3;           % number of trials shown individually
fs           = 1000;        % sampling frequency (Hz)
tlim         = [-0.3 1.5];  % trial window relative to stimulus (s)
blim         = [-0.4 0.8];  % window around each burst trough for the BRP (s)
aperExp      = -1.5;        % aperiodic exponent of the power spectrum (-1 pink, -2 brown;
                            % electrophysiological data typically -0.5 to -2)
bf           = 9;           % burst frequency (Hz)
bA           = 1.5;         % burst amplitude (relative to the normalised background)
bCyc         = 3.5;         % burst duration (cycles)
bCycVar      = 0;           % random variability in burst duration (cycles); 0 = fixed
jitter_max   = 0.25;        % maximal stimulus-to-burst delay for random trials (s),
trialSpacing = 2;           % vertical spacing between plotted trials (units of max amplitude)
bpfreq       = [0.95 30];   % band-pass applied to the background (Hz)
mod_freq     = 5;           % low-pass cut-off for amplitude modulation within bursts (Hz)

% Colours
col_burst = [0.20 0.45 0.75];   % blue, burst-locked
col_stim  = [0.85 0.33 0.10];   % orange-red, stimulus-locked
col_grey  = [0.50 0.50 0.50];   % individual trials

%% Time axes
t  = tlim(1):1/fs:tlim(2);      % trial time vector (s)
nt = numel(t);                  % samples per trial
N  = nt * T;                    % total number of samples
t0 = round(-tlim(1)*fs) + 1;    % sample index of the stimulus (t = 0)

%% Aperiodic background
% Approximate 1/f noise: power-law amplitude spectrum with random phases,
% inverse FFT, real part taken. Adequate for illustration; this is not an
% exact Hermitian-symmetric construction.
assert(mod(N,2) == 0, 'numel(t)*T must be even.');
f      = 1:N/2;
S      = f.^(aperExp/2);                        % amplitude spectrum (power ~ f^aperExp)
phases = exp(1i*2*pi*rand(1,N/2));
spec   = [S.*phases, fliplr(conj(S.*phases))];
noise  = real(ifft(spec));
[bBP, aBP] = butter(4, bpfreq/(fs/2));          % band-pass Butterworth
noise  = filtfilt(bBP, aBP, noise);
noise  = noise / max(abs(noise));
data   = reshape(noise, nt, T)';                % trials x samples

%% Bursts
% Burst duration in samples. The rand call is kept even when bCycVar = 0 so
% that the random sequence, and hence the figure, is unchanged.
bdur = round((bCyc + bCycVar*rand(1,T)) / bf * fs);

% Burst onset delay relative to the stimulus (samples). Trials 4..T are
% random; the first three (the ones plotted individually) are set by hand
% for visual clarity.
J    = rand(1,T) * jitter_max * fs;
J(1) = fs/bf;            % one cycle after the stimulus
J(2) = J(1)*4.2;
J(3) = J(1)*2.1;
J    = round(J);
assert(all(t0 + J + bdur <= nt), 'A burst extends beyond the trial window.');

[bLP, aLP] = butter(6, mod_freq/(fs/2));        % low-pass Butterworth for the modulation

ERP = data;
MSK = zeros(size(data));                        % burst envelopes, for shading
for tr = 1:T
    g     = gausswin(bdur(tr))';  g = g/max(g); % Gaussian burst envelope, peak = 1
    moder = filtfilt(bLP, aLP, rand(size(g)));  % slow amplitude modulation
    env   = g + moder;  env = env/max(env);
    inds  = t0 + J(tr) + (1:bdur(tr));
    burst = bA * sin(2*pi*bf*(1:bdur(tr))/fs);
    ERP(tr,inds) = ERP(tr,inds) + burst .* env;
    MSK(tr,inds) = g;
end
normer = max(ERP(:));                           % scaling used for plotting

%% Shared plotting variables
r    = 6;  c = 2;                               % subplot grid
smr  = 1/3;                                     % fraction of rows for the continuous panel
nTop = floor(r*smr);                            % rows for the continuous panel

cnt  = reshape(ERP(1:ttp,:)', 1, []);              % first ttp trials, concatenated
cntG = reshape(MSK(1:ttp,:)', 1, []) * normer*1.2; % burst envelopes, scaled for shading
tt   = (0:numel(cnt)-1)/fs + t(1);                 % time of the concatenated trace
yoff = normer * (2 + trialSpacing*(0:ttp-1))';     % vertical offsets of single trials
mrkY = repmat([0; 1.5]*normer*1.2, 1, ttp);        % y-extent of event markers

figure(5354); clf;
set(gcf, 'Units','pixels', 'Position',[100 100 1140 390]);

%% Stimulus-locked (left column)
C   = col_stim;
sps = 1:2:r*c;
stm = t0 + nt*(0:ttp-1);                        % stimulus samples in the concatenated trace

subplot(r, c, sps(1:nTop)); hold on
patch([tt fliplr(tt)], [cntG fliplr(-cntG)], C, 'EdgeColor','none', 'FaceAlpha',.75);
plot(tt, cnt, 'k', 'LineWidth',1);
plot(tt([stm; stm]), mrkY, 'Color',C, 'LineWidth',1);
style_raw_axes();

subplot(r, c, sps(nTop+1:end)); hold on
erp = flipud(ERP(1:ttp,:));                     % first trial plotted on top
plot(t, (erp + yoff)', 'Color',col_grey, 'LineWidth',1);
plot_mean_sem(t, ERP, C);
xline(0, 'Color',C, 'LineWidth',1);
xlabel('Time from stimulus');
style_mean_axes(C);

%% Burst-locked (right column)
% Each trial is aligned to its deepest trough (global minimum of the trial)
C    = col_burst;
sps  = 2:2:r*c;
bind = round(blim(1)*fs) : round(blim(2)*fs);   % samples around the trough
bt   = bind / fs;                               % time from burst (s)

BRP  = nan(T, numel(bind));
stmB = zeros(1, T);
for tr = 1:T
    [~, t0b] = min(ERP(tr,:));
    assert(t0b + bind(1) >= 1 && t0b + bind(end) <= nt, ...
        'Trial %d: BRP window exceeds the trial; adjust blim, tlim or jitter.', tr);
    BRP(tr,:) = ERP(tr, t0b + bind);
    stmB(tr)  = t0b + nt*(tr-1);                % trough sample in the concatenated trace
end
stmB = stmB(1:ttp);

% Same continuous data as the left column, now marked at the burst troughs
subplot(r, c, sps(1:nTop)); hold on
patch([tt fliplr(tt)], [cntG fliplr(-cntG)], C, 'EdgeColor','none', 'FaceAlpha',.75);
plot(tt, cnt, 'k', 'LineWidth',1);
plot(tt([stmB; stmB]), mrkY, 'Color',C, 'LineWidth',1);
style_raw_axes();

subplot(r, c, sps(nTop+1:end)); hold on
brp = flipud(BRP(1:ttp,:));                     % first trial plotted on top
plot(bt, (brp + yoff)', 'Color',col_grey, 'LineWidth',1);
plot_mean_sem(bt, BRP, C);
xline(0, 'Color',C, 'LineWidth',1);
xlabel('Time from burst');
style_mean_axes(C);

%% Local functions
function plot_mean_sem(x, X, C)
% Mean across rows of X with a shaded +/- SEM band
m = mean(X);
e = std(X) / sqrt(size(X,1));
patch([x fliplr(x)], [m+e fliplr(m-e)], C, 'EdgeColor','none', 'FaceAlpha',.3);
plot(x, m, 'Color',C, 'LineWidth',1);
end

function style_raw_axes()
set(gca, 'XColor','none', 'YColor','none');
axis tight
end

function style_mean_axes(C)
set(gca, 'XTick',0, 'YColor','none', 'XColor',C);
axis tight
end