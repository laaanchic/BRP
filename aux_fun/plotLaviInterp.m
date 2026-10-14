function [abba2, f2, lavi2] = plotLaviInterp (f_lavi, lavi)
[f2, lavi2] = interLog(f_lavi, lavi);
% f2 = f_lavi; lavi2=lavi;
abba2 = ABBA(lavi2, f2); abba2=abba2{1};
yline(median(lavi2),':k');
c2 = [108,192,12; 83,174,244]/255;
for bi = 1:size(abba2,1)
    if abba2(bi,9)==1 
        c = c2(1,:); 
    else
        c = c2(2,:); 
    end
    inds = abba2(bi,1):abba2(bi,2);
    plot(f2(inds), lavi2(inds), 'linewidth', 2, 'color', c);
end
plot(f2, lavi2, 'k', 'linewidth', .5);
set(gca,'xscale','log','xtick',[2:2:10,15,20,30],'xlim',[3.1,31],'tickdir','out');

function [x2, y2] = interLog(x,y)
logx = log10(x);
logx_new = linspace(min(logx), max(logx), 10*numel(x));
x2 = 10.^logx_new;
y2 = interp1(logx, y, logx_new, 'pchip');   % or 'spline' 'pchip'
end
end