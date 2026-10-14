function [Output, pmtr] = Bursts_detection_BRP(raw, giv, pmtr)
% BURSTS_DETECTION_BRP detects oscillatory bursts in a single LFP/EEG channel.
%
% Bursts are local maxima in the time-frequency power (Morlet wavelets) that
% exceed a per-frequency percentile threshold. Burst onset/offset is defined
% by where power drops below a (possibly different) percentile. Bursts that
% are too short, too wide in frequency, or that overlap in time with a
% stronger burst at a nearby frequency are removed.
%
% Reference: Schuma, Landau, Diester, Anderson, & Karvat: Burst-related potentials as 
%            a temporal anchor for cognition, Comm. Biol. 2026
%
% INPUTS
%   raw   - FieldTrip raw data structure with ONE channel and ONE trial
%           (fields: trial, time, label, fsample, sampleinfo)
%   giv   - structure; giv.fs = sampling frequency (Hz). If empty or missing
%           the field, raw.fsample is used.
%   pmtr  - (optional) structure with detection parameters (defaults below):
%     .foi          frequencies of interest (Hz)                  [15:30]
%     .ptl_high     percentile threshold for burst peaks           [90]
%     .ptl_dur      percentile threshold for burst duration        [ptl_high]
%     .SD_dur       SDs added to ptl_dur threshold                 [0]
%     .min_dur      minimal burst duration (cycles)                [1]
%     .wave_width   wavelet width (cycles)                         [7]
%     .muscle       sample indices of muscle artefacts to remove   []
%     .art_thresh   amplitude threshold for artefacts (data units, e.g. uV) [250]
%     .lastSample   ignore bursts peaking after this sample        [end of data]
%     .min_f_gap    min. frequency gap (Hz) between overlapping bursts, per foi
%                                                   [ceil(foi/4), minimum 4]
%     .max_span     max. allowed frequency span of a burst (Hz)    [foi(end)-foi(1)]
%     .output_power if 1, return the TFR (median-normalized) and the
%                   band-pass filtered data in Output                [0]
%
% OUTPUT
%   Output.burstDists  - bursts x 18 matrix; columns:
%      1 = peak frequency (Hz)            2 = peak time (sample)
%      3 = begin (sample)                 4 = end (sample)
%      5 = nearest trough (sample)        6 = duration (cycles)
%      7 = duration (ms)                  8 = power (relative to percentile)
%      9 = power (absolute)              10 = number of troughs
%     11 = number of peaks               12 = begin, high percentile (sample)
%     13 = end, high percentile (sample) 14 = frequency index (internal)
%     15 = lowest span frequency (Hz)    16 = highest span frequency (Hz)
%     17 = lowest peak frequency (Hz)    18 = highest peak frequency (Hz)
%   Output.TroughPeak  - per burst: {trough samples, peak samples}
%   Output.PeakSpan    - per burst: frequencies of the power maxima over time
%   Output.arts        - logical vector marking artefact samples
%   Output.col_names   - names of the burstDists columns
%   Output.pmtr        - parameters used
%   Output.pwr, Output.filtd - only if pmtr.output_power == 1
%
% Requires: FieldTrip; Image Processing Toolbox (imregionalmax);
%           Statistics and Machine Learning Toolbox (prctile).
% Helper functions: gross_artifact_removal, get_burst_durations,
%                   getFreqSpan, dispRMVprev, showCurrnetTime

t_func = tic;

%% Input checks and parameters
assert(numel(raw.label) == 1 && numel(raw.trial) == 1, ...
    'Input must be a single-channel, single-trial FieldTrip structure');
if nargin < 2 || isempty(giv) || ~isfield(giv, 'fs'); fs = raw.fsample; else; fs = giv.fs; end
if nargin < 3; pmtr = []; end
if ~isfield(pmtr,'foi');          pmtr.foi          = 15:30; end
if ~isfield(pmtr,'ptl_high');     pmtr.ptl_high     = 90; end
if ~isfield(pmtr,'ptl_dur');      pmtr.ptl_dur      = pmtr.ptl_high; end
if ~isfield(pmtr,'SD_dur');       pmtr.SD_dur       = 0; end   % Little et al. 2019 used 1.75
if ~isfield(pmtr,'min_dur');      pmtr.min_dur      = 1; end
if ~isfield(pmtr,'wave_width');   pmtr.wave_width   = 7; end
if ~isfield(pmtr,'muscle');       pmtr.muscle       = []; end
if ~isfield(pmtr,'art_thresh');   pmtr.art_thresh   = 250; end
if ~isfield(pmtr,'lastSample');   pmtr.lastSample   = raw.sampleinfo(end); end
if ~isfield(pmtr,'min_f_gap');    pmtr.min_f_gap    = max(ceil(pmtr.foi./4), 4); end
if ~isfield(pmtr,'max_span');     pmtr.max_span     = pmtr.foi(end) - pmtr.foi(1); end
if ~isfield(pmtr,'output_power'); pmtr.output_power = 0; end

foi = pmtr.foi;
nF  = length(foi);
assert(numel(pmtr.min_f_gap) == nF, 'pmtr.min_f_gap must have one value per foi');
%% Time-frequency representation
% One extra frequency is added below and above foi, so that local maxima at
% the edges of foi can be detected; these rows are removed afterwards.
f1 = max([1, foi(1) - diff(foi([1,2]))]);
f2 = foi(end) + diff(foi([nF-1,nF]));

disp(['Started calculating power at ' showCurrnetTime])
cfg        = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foi    = [f1 foi f2];
cfg.toi    = 'all';
cfg.width  = pmtr.wave_width;
cfg.pad    = 'nextpow2';
try
    TFR = ft_freqanalysis(cfg, raw);
catch
    warning('Padding to the next power of 2 failed; computing without padding');
    cfg = rmfield(cfg, 'pad');
    TFR = ft_freqanalysis(cfg, raw);
end
pwr_raw = squeeze(TFR.powspctrm);
clear TFR

%% Artefact removal
cfg             = [];
cfg.in_channel  = 1;
cfg.fs          = fs;
cfg.muscle      = pmtr.muscle;
cfg.art_thresh  = pmtr.art_thresh;
cfg.art_removal = round(fs*0.5);
pwr_clean = gross_artifact_removal(cfg, raw, pwr_raw);
arts      = isnan(gross_artifact_removal(cfg, raw, raw.trial{1}));
clear pwr_raw

%% Burst peaks: local maxima above the percentile threshold
disp(['Started calculating burst peaks at ' showCurrnetTime])
t_BW = tic;
pwr_pctl_high = prctile(pwr_clean, pmtr.ptl_high, 2);
pwr_pctl_dur  = prctile(pwr_clean, pmtr.ptl_dur, 2) + pmtr.SD_dur*std(pwr_clean, 0, 2, 'omitnan');

pwr_masked = pwr_clean ./ pwr_pctl_high;
pwr_masked(pwr_masked < 1 | isnan(pwr_masked)) = 0;

pwr_for_dur = pwr_clean ./ pwr_pctl_dur;
pwr_for_dur(pwr_for_dur < 1 | isnan(pwr_for_dur)) = 0;

BW = imregionalmax(pwr_masked);
BW([1,end],:) = [];            % remove the padding frequencies
[i, j] = find(BW);
i(j > pmtr.lastSample) = [];   % remove bursts after the last sample of interest
j(j > pmtr.lastSample) = [];

pwr_clean([1,end],:)     = [];
pwr_masked([1,end],:)    = [];
pwr_for_dur([1,end],:)   = [];
pwr_pctl_high([1,end],:) = [];
disp(['Calculating burst peaks took ' num2str(toc(t_BW)) ' seconds']);

%% Band-pass filtered data, to detect troughs and peaks
disp(['Started filtering at ' showCurrnetTime])
t_filter = tic;
dat   = raw.trial{1};
filtd = zeros(nF, length(dat));
for f = 1:nF
    try
        filtd(f,:) = ft_preproc_bandpassfilter(dat, raw.fsample, foi(f) + [-2 2]);
    catch
        filtd(f,:) = bandpass(dat, foi(f) + [-2 2], raw.fsample); % requires Signal Processing Toolbox
    end
end
disp(['Filtering took ' num2str(toc(t_filter)) ' seconds']);

%% Characterize each burst
disp(['Started burstDists structure at ' showCurrnetTime])
nB         = length(j);
burstDists = zeros(nB, 18);
PeakSpan   = cell(nB, 1);
outTick    = cell(nB, 2);       % {troughs, peaks}
too_short  = false(nB, 1);
n          = size(pwr_masked, 2);
disp(['Total bursts: ' num2str(nB)]);
t_trials = tic;
prev = fprintf(' ');

for x = 1:nB
    if ~mod(x, 10000); prev = dispRMVprev(['Finished ' num2str(x) ' bursts in ' num2str(toc(t_trials)) ' sec'], prev); end
    burstDists(x,1)  = foi(i(x));  % peak frequency
    burstDists(x,14) = i(x);       % frequency index
    burstDists(x,2)  = j(x);       % peak time

    % begin/end according to the peak-detection percentile
    [t_before_high, t_after_high, ~, ~] = get_burst_durations(i, j, x, pmtr, i, filtd, n, pwr_masked);
    % begin/end according to the duration percentile
    [t_before, t_after, pks, trgh] = get_burst_durations(i, j, x, pmtr, i, filtd, n, pwr_for_dur);

    burstDists(x,3)  = t_before;
    burstDists(x,4)  = t_after;
    burstDists(x,12) = t_before_high;
    burstDists(x,13) = t_after_high;
    burstDists(x,6)  = (t_after - t_before) ./ ceil(fs ./ burstDists(x,1));
    too_short(x)     = burstDists(x,6) < pmtr.min_dur;

    if ~too_short(x)
        burstDists(x,7)  = (t_after - t_before) / fs * 1000;
        burstDists(x,8)  = pwr_masked(i(x), j(x));
        burstDists(x,9)  = pwr_clean(i(x), j(x));
        outTick{x,1}     = trgh;
        outTick{x,2}     = pks;
        burstDists(x,10) = length(trgh);
        burstDists(x,11) = length(pks);

        % nearest trough to the power maximum (search +-1.5 cycles)
        wished_inds = ceil(fs / burstDists(x,1) * 1.5);
        inds = max([1, j(x) - wished_inds]) : min([length(dat), j(x) + wished_inds]);
        sig  = filtd(i(x), :);
        locs = [];
        val  = [];
        for si = 2:(length(inds)-1)   % faster than findpeaks for short segments
            if sig(inds(si)) < sig(inds(si)-1) && sig(inds(si)) < sig(inds(si)+1)
                locs = [locs, inds(si)];          %#ok<AGROW>
                val  = [val, abs(sig(inds(si)))]; %#ok<AGROW>
            end
        end
        dist = abs(locs - burstDists(x,2));
        nearest_trough = locs(dist == min(dist));
        if length(nearest_trough) > 1     % equidistant troughs: take the deeper one
            [~, tmpi] = max(val(dist == min(dist)));
            nearest_trough = nearest_trough(tmpi);
        end
        if isempty(nearest_trough)        % (very rare) no trough: use the power maximum
            nearest_trough = burstDists(x,2);
        end
        burstDists(x,5) = nearest_trough;
    end

    % frequency span of the burst
    col = pwr_masked(:, burstDists(x,2));
    [~, ind_above, ind_below] = getFreqSpan(col, i(x));
    burstDists(x,15) = foi(i(x) - ind_below);
    burstDists(x,16) = foi(i(x) + ind_above);
    f_ind = (i(x) - ind_below) : (i(x) + ind_above);
    if numel(f_ind) > 1
        curr = pwr_masked(f_ind, burstDists(x,3):burstDists(x,4));
        [M, I] = max(curr);
        I(M == 0) = [];
        burstDists(x,17) = foi(f_ind(min(I)));
        burstDists(x,18) = foi(f_ind(max(I)));
        PeakSpan{x} = foi(f_ind(I));
    else   % single-frequency burst
        burstDists(x,[17,18]) = foi(i(x));
        PeakSpan{x} = foi(i(x));
    end
end

%% Remove bursts that are too short (time) or too wide (frequency)
too_wide = (burstDists(:,16) - burstDists(:,15)) > pmtr.max_span;
bst = burstDists;
bst(too_short | too_wide, :)   = [];
outTick(too_short | too_wide, :) = [];
PeakSpan(too_short | too_wide)   = [];

%% Prune overlapping bursts at nearby frequencies (keep the higher-energy one)
pwr = pwr_clean ./ pwr_pctl_high;
disp('.');
prev = fprintf(' ');
t_prune = tic;
for bi = 1:size(bst, 1)
    if ~mod(bi, 2500); prev = dispRMVprev(['Finished pruning ' num2str(bi) ' bursts in ' num2str(toc(t_prune)) ' sec'], prev); end
    if bst(bi,1) > 0   % skip bursts that were already rejected
        cond1 = bst(:,3) >= bst(bi,3) & bst(:,3) <= bst(bi,4); % begins within the current burst
        cond2 = bst(:,4) <= bst(bi,4) & bst(:,4) >= bst(bi,3); % ends within the current burst
        cond3 = bst(:,3) <= bst(bi,3) & bst(:,4) >= bst(bi,4); % contains the current burst
        cooc = cond1 | cond2 | cond3;
        cooc(bi) = 0;
        cooc = find(cooc);
        for ci = 1:length(cooc)
            cj = cooc(ci);
            if abs(bst(bi,1) - bst(cj,1)) <= pmtr.min_f_gap(bst(bi,14))            
                % energy estimated as relative power x duration; if equal, keep bi
                e1 = pwr(bst(bi,14), bst(bi,2)) * (bst(bi,4) - bst(bi,3));
                e2 = pwr(bst(cj,14), bst(cj,2)) * (bst(cj,4) - bst(cj,3));
                if e2 > e1
                    bst(cj,3) = min([bst(cj,3), bst(bi,3)]);
                    bst(cj,4) = max([bst(cj,4), bst(bi,4)]);
                    bst(bi,:) = 0;   % reject bi
                    break            % no need to compare a rejected burst further
                else
                    bst(bi,3) = min([bst(cj,3), bst(bi,3)]);
                    bst(bi,4) = max([bst(cj,4), bst(bi,4)]);
                    bst(cj,:) = 0;   % reject cj
                end
            end
        end
    end
end
keep = bst(:,1) > 0;
outTick  = outTick(keep, :);
PeakSpan = PeakSpan(keep);
bst      = bst(keep, :);
disp('.');
disp(['Total bursts rejected: ' num2str(size(burstDists,1) - size(bst,1))]);

%% Output
Output.pmtr       = pmtr;
Output.burstDists = bst;
Output.TroughPeak = outTick;
Output.PeakSpan   = PeakSpan;
Output.arts       = arts;
Output.col_names  = {'maxima freq (Hz)', 'maxima timing (sample)', ...
    'begin (sample)', 'end (sample)', 'nearest trough', 'duration (period)', ...
    'duration (ms)', 'power (rel pctl)', 'power (absolute)', 'number of troughs', 'number of peaks', ...
    'begin high percentile (sample)', 'end high percentile (sample)', ...
    'ignore', 'low Span freq (Hz)', 'high Span freq (Hz)', 'low Peak freq (Hz)', 'high Peak freq (Hz)'}';
if pmtr.output_power
    Output.pwr   = pwr_clean ./ median(pwr_clean, 2, 'omitnan');
    Output.filtd = filtd;
end
disp('.');
disp(['The call to ' mfilename ' took ' num2str(toc(t_func)) ' seconds']);
end