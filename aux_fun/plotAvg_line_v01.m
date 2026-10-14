function [m, inds] = plotAvg_line_v01(data,C,z,ylm,do_bl,factor,examp,previnds)
% C = colormap
% z = color
do_bl = 1;
factor = 1;
examp = 0;
previnds=[];
% if nargin<3, z=1; end
% if isempty(z), z=1; end
% if size(C,1)<2, z=1; end
% if nargin<5, do_bl=1; end
% if nargin<6, factor=1; end
% if nargin<7, examp=0; end
% if nargin<8, previnds=[]; end
clr = C(z,:);

x = data.time{1};
y = squeeze(cat(3, data.trial{:}))*factor;

% baseline correction to the first 10%
if do_bl
    blind = round(0.2*length(x));
    bl = mean(y(1:blind,:));
    y=y-bl;
end

% % smooth a little
% sm_fact = round(0.02*length(x));
% % smoother = gausswin(sm_fact);
% smoother = ones(1,sm_fact);
% smoother = smoother/sum(smoother);
% % y = conv2(smoother,1,y,'same');
m = mean(y,2)';
e = std(y,[],2)';
% v = var(y,[],2)';
% e = sqrt(v / size(y,2));
% e = sqrt(v);
% E = [m+e; m-e];
% [~,~,E] = ttest(y');
% patch([x flip(x)], [m-e, flip(m+e)], clr, 'edgecolor','none','facealpha',0.1)
% patch([x flip(x)], [E(1,:), flip(E(2,:))], clr, 'edgecolor','none','facealpha',0.3)
if examp % plot individual examples of bursts
    if isempty (previnds)
        [~,I] = sort(corr(y,m'),'descend');
        n = length(I);
        inds = round(n.*(0.1:0.1:0.25));
        inds = I(inds);
    else inds = previnds;
    end
    for i = 1:length(inds)
        cc = clr + 0.15*i;
        cc(cc>1)=1;
        plot(x, y(:,inds(i)), 'color', cc, 'linewidth', 0.25)
%         plot(x, y(:,inds(i)), 'color', [0 0 0]+0.4+0.1*i, 'linewidth', 0.5)
    end
end
plot(x,m,'color',clr,'linewidth',0.5)
ytck = get(gca,'ytick');
set(gca,'xlim',x([1 end]));%, 'xtick',[],'ytick', ytck(2:end));
if nargin>3 && ~isempty(ylm)
    ylim(ylm)
end
grid on
