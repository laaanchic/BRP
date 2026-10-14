function output = ft_GK_redefinefreq_v04(cfg, freq)
% Redfining trials of one long freq structure.
% Assuming all trials have the same duration.

if ~isfield(cfg, 'fs'); fs = 500; else; fs = cfg.fs; end
if ~isfield(cfg, 'baselinewindow'); do_BL = false; else; do_BL = true; end
if isfield(cfg,'resamp') % the original freq structure was made with "resampling"
    orig_fs = cfg.resamp.orig_fs;
    prev_time = cfg.resamp.prev_time;
    curr_ind = cfg.resamp.curr_ind;  
else
    orig_fs = fs;
    prev_time = freq.time;
    curr_ind = 1:length(freq.time);
end
if ~isfield(cfg, 'trl') 
    output=freq;
    warning('cfg.trl was not defined. The output is identical to the input');
else
    if ~isfield(freq,'powspctrm'), freq.powspctrm=freq.fourierspctrm; nameflag=1;
    else nameflag=0; end
    output = freq;
    trl = cfg.trl;
    ind1 = nearest(curr_ind,trl(1,1));
    ind2 = nearest(curr_ind,trl(1,2));
    offf = trl(1,3)/orig_fs;
    t = ((ind1:ind2)-ind1)/fs+offf;
    
    nTr = size(trl,1);
    nCh = length(freq.label);
    nFr = length(freq.freq);    
    nt = ind2-ind1;
    siz = [nTr, nCh, nFr, nt+1];    
    new_pow = zeros(siz); 
    if ndims(freq.powspctrm)< 4
        tmp = zeros(1,nCh,nFr,length(freq.powspctrm));
        tmp(1,:,:,:) = freq.powspctrm;
        freq.powspctrm = tmp;
        clear tmp
    end
    for tr = 1:nTr
        ind1 = nearest(curr_ind,trl(tr,1));
        ind2 = ind1+nt; %nearest(curr_ind,trl(tr,2));        
        new_pow(tr,:,:,:) = freq.powspctrm(1,:,:,ind1:ind2);
    end
    
    if do_BL
        blind = t>=cfg.baselinewindow(1) & t<=cfg.baselinewindow(end);
        bl = nanmean(new_pow(:,:,:,blind),4);
        if isempty (cfg.blmethod), cfg.blmethod='absolute'; end
            switch cfg.blmethod
                case 'absolute'
                    new_pow = new_pow - bl;
                case 'relative'
                    new_pow = new_pow ./ bl;
                case 'dB'
                    new_pow = 10*log10(new_pow ./ bl);
            end
    end
    output.powspctrm = new_pow;
    output.time = t;
    if size(trl,2)>3; output.trialinfo = trl(:,4:end); end
    output.cumtapcnt = ones(nTr, nFr);
    if nameflag, output.fourierspctrm = output.powspctrm; output = rmfield(output,'powspctrm'); end
end