function cfc = timefreq_coupling(cfg)
% f is a cell array of spm meeg files
%
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

try cfg;         catch; return;        end

try fvec = cfg.freqvec; catch; fvec = [cfg.F{1}(1) cfg.F{1}(end)]; end
try tvec = cfg.timevec; catch; tvec = [cfg.T{1}(1) cfg.T{1}(end)]; end
try cvec = cfg.conditn; catch; cvec = 1:size(cfg.Y,2);             end

Y = cfg.Y;
F = cfg.F;
T = cfg.T;

%y = squeeze(cat(4,Y{:}));
y = inner(Y);
f = F{1};
t = T{1};

% returns
cfc.f = f;
cfc.t = t;
cfc.T = T;
cfc.y = y;

% subspace time
time = t;
T    = tvec;%/1000;
T    = [findthenearest(T(1),time):findthenearest(T(2),time)];
fprintf('using time window %d to %d\n',tvec(1),tvec(2));

% subspace frequency
freq = f;
F    = fvec;
F    = [findthenearest(F(1),freq):findthenearest(F(2),freq)];
fprintf('using frequency window %d to %d\n',fvec(1),fvec(2));

% roi/foi
y    = y(:,:,T,F);

%
Q  = @squeeze;
m  = @mean;
nf = size(y,4);

yY = y;
yY = permute(yY,[2 1 3 4]);
yY = Q(m(yY(:,:,:,:),3));

% Conditions [subspc]
CC = 1:size(yY,1);
CC = CC(cvec);


for i = 1:length(CC)
    n = 0;
    for f1 = 1:nf
        for f2 = 1:nf
            n = n + 1;
            fprintf('running corr %d of %d for condition %d of %d\n',n,(nf*nf),i,length(CC));
            [r{i}(f1,f2),p{i}(f1,f2)] = corr(Q(yY(CC(i),:,f1))',Q(yY(CC(i),:,f2))');
        end
    end
end

DoTFPlot(t(T),f(F),y)
DoCFCPlot(r,f(F));
DoCFCPlot_p(p,f(F));

% cfc returns
cfc.cfc.r = r;
cfc.cfc.p = p;
cfc.cfc.f = f(F);
cfc.cfc.t = t(T);
cfc.cfc.y = y;

return

end

function DoTFPlot(t,f,y)
% plot time frequency averaged over group

figure,
nc = size(y,2);
%ha = tight_subplot(2,nc/2,[.01 .03],[.1 .01],[.01 .01]);

for i = 1:nc
    %axes(ha(i)),
    subplot(2,nc/2,i),
    M = mean(y,1);
    M = max(M(:));
    plotbert(t,f,squeeze(mean(y(:,i,:,:),1))',M,['Condition ',num2str(i)]);
    set(gcf,'inverthardcopy','off');
    whitebg(1,'k'); %alpha(.5);
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
    whitebg(1,'k'); alpha(.5);
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
    whitebg(1,'k'); alpha(.5);
end
whitebg(1,'k'); alpha(.5);
set(gcf,'inverthardcopy','off');

end
