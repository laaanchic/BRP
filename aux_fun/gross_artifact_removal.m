function out_data = gross_artifact_removal(cfg, in_reference, in_data)
% Detects "artifacts" in in_reference as points in which the absolute value 
% exceeds cfg.art_thresh, and removes cfg.art_removal points before and 
% after them (puts NaNs) from in_data and outputs to ouT_data.
% Use as out_data = gross_artifact_removal(cfg, in_reference, in_data)
% Parameters:
% -----------
% in_reference      = N_chan x N_samples double, the reference LFP matrics.
%
% in_data           = N_chan/freq x N_samples doubel, the matrix to remeove
%                     artifacts from.
%
% cfg.art_thresh    = double, the artifact threshold in uV. Default = 200;
%
% cfg.art_removal   = integer, number of sanmples before and after an
%                     artifact to be removed. Default = 500.
% 
% cfg.in_channel    = integer. The channel from in_reference to be used for 
%                     artifact detection. Default = 1.
%
% cfg.monitor_artifacts = bool, if true each artifact will be plotted for
%                         inspection. Default = 0.
%
% cfg.fs            = double, the sampling frequency in Hz (relevant only
%                     for plotting). Default = 976.5625 Hz.
%
% Last updated: 26/05/2019, by Golan Karvat

% Defaults:
if ~isfield(cfg, 'art_thresh'); cfg.art_thresh = 200; end
if ~isfield(cfg, 'art_removal'); cfg.art_removal = 500; end
if ~isfield(cfg, 'in_channel'); cfg.in_channel = 1; end
if ~isfield(cfg, 'monitor_artifacts'); cfg.monitor_artifacts = 0; end
if ~isfield(cfg, 'fs'); cfg.fs = 500; end %976.5625

try ref = in_reference.trial{1};
catch ref = in_reference; end
tic
art_thresh = cfg.art_thresh;
art_removal = cfg.art_removal;
arts = find(abs(ref(cfg.in_channel,:))>art_thresh);
try
    for ai = 1:size(cfg.muscle,1) % add muscle artifacts
        arts = [arts, cfg.muscle(ai,1):cfg.muscle(ai,2)];
    end
catch disp('No muscles artifact structure provided'); end
arts = sort(unique(arts));
fs = cfg.fs;
out_data = in_data;
prev_nan = sum(isnan(out_data(1,:)));

if cfg.monitor_artifacts
    bla = ref(cfg.in_channel,:);
    a = find(diff(arts)>5)+1;
    figure(12); clf
    for i = 1:length(a)
        clf; hold on
        inds = arts(a(i))+[-art_removal:art_removal];
        plot(inds/fs,bla(inds));
        inds2 = find(abs(bla(inds))>art_thresh);
        %         inds2 = find(isnan(bla(inds)));
        scatter (inds(inds2)/fs,bla(inds(inds2)),'r','fill');
        %         scatter (inds(inds2)/fs,zeros(size(inds2)),'r','fill');
        title(num2str(i)); set(gca,'xlim',[inds(1) inds(end)]/fs);
        disp(['Figure ' num2str(i) ' out of ' num2str(length(a)) ' , press any key to continue']);
        pause;
    end
end

for a = arts
    if a <= art_removal
        out_data(:, [1 : a + art_removal]) = nan;
    elseif a+art_removal >= size(out_data,ndims(out_data)) %always looking at the last dimension
        out_data(:, [a - art_removal : end]) = nan;
    else
        out_data(:,a + [-art_removal : art_removal]) = nan;
    end
end
post_nan = sum(isnan(out_data(1,:)));
% disp (['Detecting ' num2str(length(arts)) ' atrifact samples and removing ' num2str(post_nan-prev_nan) ' samples took ' num2str(toc) ' seconds']);