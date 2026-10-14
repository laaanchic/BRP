function output = findBurstMinima(data, bst, fs)

output = bst;
for bi = 1:size(bst,1)
    sw = ceil(fs/bst(bi,1));
    inds = bst(bi,2) + (-sw:sw);
    [~,I] = min(data.trial{1}(inds));
    output(bi,15) = bst(bi,2)-sw+I-1;    
end

