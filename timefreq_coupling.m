function cfc = timefreq_coupling(f,cfg)
% f is a cell array of spm meeg files
%
% computes time-freq and cross freq coupling 
% [and MMN if that's of use] of spm meeg files.
%
% optionals cfg struct with optional
% cfg.sensors  = reduction method: pca, mean or sep
% cfg.resamp   = up/down sample by fixed val (+/- n)
% cfg.freqs    = frequency vector, eg. 1:.25:100
% cfg.baseline = do baselining [0/1]
% cfg.basetime = start / stop in secs, eg. [-.1 1]
% cfg.compute  = {'abs','pow','phase','complex','plv'}
% cfg.MMN      = [0/1]
%
%
%
%


try cfg; catch; cfg = []; end

[MMN,t,F,T,Y] = DoMMN(f,cfg);
y = squeeze(cat(3,Y{:}));
y = inner(y);
f = F{1}{1};
t ;

% set time
time = t;
T    = [100 300]/1000;
T    = [findthenearest(T(1),time):findthenearest(T(2),time)];

%
Q  = @squeeze;
m  = @mean;
nf = size(y,4);

yY = y;
yY = permute(yY,[2 1 3 4]);
yY = Q(m(yY(:,:,T,:),3));

CC = 10;

for i = 1:CC
    n = 0;
    for f1 = 1:nf
        for f2 = 1:nf
            n = n + 1;
            fprintf('running corr %d of %d for condition %d of %d\n',n,(nf*nf),i,size(yY,2));
            [r{i}(f1,f2),p{i}(f1,f2)] = corr(Q(yY(:,i,f1)),Q(yY(:,i,f2)));
        end
    end
end

DoTFPlot(t,f,y)
DoCFCPlot(r,f);
DoCFCPlot_p(p,f);


cfc.r = r;
cfc.p = p;
cfc.f = f;
cfc.t = t;
cfc.T = T;
cfc.y = y;

return

end

function DoTFPlot(t,f,y)
% plot time frequency averaged over group

figure,
nc = size(y,2);
ha = tight_subplot(2,nc/2,[.01 .03],[.1 .01],[.01 .01]);

for i = 1:nc
    axes(ha(i)),
    %subplot(2,nc/2,i),
    M = max(y(:));
    plotbert(t,f,squeeze(mean(y(:,i,:,:),1))',M,['Condition ',num2str(i)]);
    set(gcf,'inverthardcopy','off');

end
end

function DoCFCPlot(r,f)
% plot correlation coefficients 

n = length(r);
figure,

k = round(length(f)*.015);
h = @(x)HighResMeanFilt(x,1,k);
fprintf('smoothing with kern of %d\n',k);

for i = 1:n
    subplot(2,round(n/2),i),
    imagesc(h(r{i}));
    set(gca,'YDir','normal');
    title(['Frequency coupling: trial ',num2str(i)],'fontsize',18);

    NumTicks = 8;
    box off ; alpha(.5);
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks),'XTickLabel',round(linspace(f(1),f(end),NumTicks)));
    set(gca,'YTick',linspace(L(1),L(2),NumTicks),'YTickLabel',round(linspace(f(1),f(end),NumTicks)));
    set(gca,'fontsize',18);
    
end
whitebg(1,'k'); alpha(.5);
set(gcf,'inverthardcopy','off');

end

function DoCFCPlot_p(p,f)
% plot p-values thresholded to <.1, otherwise dark blue

n = length(p);
figure,

k = round(length(f)*.015);
h = @(x)HighResMeanFilt(x,1,k);
fprintf('smoothing with kern of %d\n',k);

for i = 1:n
    subplot(2,round(n/2),i),
    A = (p{i});
    A(A>.1) = inf;
    imagesc(h(A));
    caxis([0 .1])

    colormap(flipud(jet))
    
    set(gca,'YDir','normal');
    title(['Frequency coupling: trial ',num2str(i)],'fontsize',18);

    NumTicks = 8;
    box off ; alpha(.5);
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks),'XTickLabel',round(linspace(f(1),f(end),NumTicks)));
    set(gca,'YTick',linspace(L(1),L(2),NumTicks),'YTickLabel',round(linspace(f(1),f(end),NumTicks)));
    set(gca,'fontsize',18);
    
end
whitebg(1,'k'); alpha(.5);
set(gcf,'inverthardcopy','off');

end