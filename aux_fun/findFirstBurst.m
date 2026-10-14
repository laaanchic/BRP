function [trlOut, comm] = findFirstBurst(trl,bst,fs,lockto)
trlOut = []; % sample of the furst burst in the trial
comm = []; % commulative timpoints of bursts
for ti = 1:size(trl,1)
    beg = trl(ti,1)-trl(ti,3); % rletive to 0
    flsh = beg + round(fs*trl(ti,9));
    b1 = bst(bst(:,lockto)>=beg & bst(:,lockto)<=flsh, lockto);
    if ~isempty(b1) % only first burst in the period
        trlOut = [trlOut; b1(1)];
        comm = [comm; b1(1)-beg];
    end 
end
end