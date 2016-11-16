function [F,T,Y,FT] = SPM_2TF(D,cfg)
% TF for SPM MEEG Data which is either robust averaged or not
% Runs on both sensor (MEEG or LFP) and source data
% 
% optionals cfg struct with optional
% cfg.sensors  = reduction method: pca, mean, sep or ica
% cfg.resamp   = up/down sample by fixed val (+/- n)
% cfg.freqs    = frequency vector, eg. 1:.25:100
% cfg.baseline = do baselining [0/1]
% cfg.basetime = start / stop in secs, eg. [-.1 1]
% cfg.type     = 'evoked' or 'induced'
% cfg.compute  = {'abs','pow','phase','complex','plv'}
%
% AS


% Run over multiple subjects if D is cell
if iscell(D) && length(D) > 1
    for i = 1:length(D)
        clear dfile; 
        dfile = spm_eeg_load(D{i});
        [F(i,:),T(i,:),Y(i,:),FT(i,:)] = SPM_2TF(dfile,cfg);
    end
    return;
end


% inputs
try cfg; catch cfg = []; end

if isfield(cfg,'sensors'); meth = cfg.sensors; else meth = 'pca';    end
if isfield(cfg,'resamp');  rs   = cfg.resamp ; else rs   = 1;        end
if isfield(cfg,'freqs');   Freqs= cfg.freqs;   else Freqs= 1;        end

if isfield(cfg,'baseline');bl   = cfg.baseline;else bl   = 0;        end
if isfield(cfg,'basetime');bs   = cfg.basetime;else bs   = [-.1 .1]; end
if isfield(cfg,'type'    );type = cfg.type;    else type = 'evoked'; end
if isfield(cfg,'space'   );spce = cfg.space;   else spce = 'sensor'; end

global n ncc
n   = []; % n spatial modes if PCA
ncc = []; % n spatial modes if ICA


% foi
if length(Freqs) == 1;
    if ~isempty(D.frequencies)
         Freqs = D.frequencies;
    else Freqs = 1:.2:100;
    end
end

fprintf('Calculating for frequencies %d to %d in %d steps\n',Freqs(1),Freqs(end),(Freqs(2)-Freqs(1)));

% data source
switch spce
    case 'sensor';  SS = D;
    case 'source';  SS = GetSS(D);
end

% Data dimensions
S    = size(SS);
ns   = S(1);
nc   = S(3);
nsam = S(2);

% check whether [robust] averaged
if nc > 25; 
    nc = length(D.condlist);
    for i = 1:(nc)
        ic{i}  = D.indtrial(D.condlist{i});
        domeeg = 1;
    end
else
    ic = num2cell(1:nc);
    domeeg = 0;
end


% FT
FT          = D.fttimelock;
dt          = 1/D.fsample;
FS          = 1/dt;        % sampling frequenct
t           = D.time;      % peristim time
cfg.fsample = FS*rs;       % adjusted sample rate
t           = linspace(t(1), t(end),rs*length(t));
smooth      = 4;

% retain evoked or induced matrix
switch type
    case lower('evoked') ; TF = 'evagram'; fprintf('Calculating evoked\n');
    case lower('induced'); TF = 'agram';   fprintf('Calculating induced\n');
end


for i = 1:nc
    
    % remove non neurophys channels
    if i == 1 && domeeg;
        D = DoMEEG(D);
    end
    
    % select trial data
    Yc = SS(:,:,ic{i});   
    cfg.sampletimes  = t;
    cfg.kill_evagram = 0;
    
    % baseline
    if bl;  cfg.baseline     = 'subtract';
            cfg.basetimes    = bs;
    else    cfg.baseline     = 'none';
    end
        
    switch meth
        case 'pca';
            YY   = HighResMeanFilt(full(Yc),rs,smooth);
            TF_Y = eig_engine(YY,cfg,Freqs);

            Y{i} = double(TF_Y.(TF));
            F{1} = TF_Y.freqs;
            T{1} = TF_Y.sampletimes;
            
        case 'ica';
            YY   = HighResMeanFilt(full(Yc),rs,smooth);
            TF_Y = ica_engine(YY,cfg,Freqs);
            
            Y{i} = double(TF_Y.(TF));
            F{1} = TF_Y.freqs;
            T{1} = TF_Y.sampletimes;
            
        case 'mean';
            YY   = HighResMeanFilt(full(Yc),rs,smooth);
            YY   = mean(YY);
            TF_Y = bert_singlechannel([YY;YY],cfg,Freqs);
            
            Y{i} = double(TF_Y.(TF));
            F{1} = TF_Y.freqs;
            T{1} = TF_Y.sampletimes;
        
        case {'sep','separate'};
            for j = 1:ns
                clc;
                fprintf('Calculating TF for condition %d of %d, sensor %d of %d\n',i,nc,j,ns);
                
                % resmaple as required
                YY = HighResMeanFilt(full(Yc(j,:)),rs,smooth);
                YY = full(YY);

                % input trials*samples
                TF_Y = bert_singlechannel([YY;YY],cfg,Freqs);
                
                %Y = double(TF_Y.agram); % real
                Y{i,j} = double(TF_Y.(TF));
                
                F{1,1} = TF_Y.freqs;
                T{1,1} = TF_Y.sampletimes;
                
            end
    end
    
end
end



function Y = eig_engine(YY,cfg,Freqs)
% simple call to my svd based pca, then SM's bert_singlechan
%
%
% opts on bert_singlechan
% cfg.compute - 
%         'abs' -absolute value of the Hilbert trace (default) 
%         'pow' - Compute the power of the Hilbert trace
%         'phase' - COmpute the phase of the Hilbert trace     
%         'complex' keep the complex elements
%         'plv' - compute the PLV across


if  ndims(YY) == 2;
    % pricipal eigenmode over YY'
    fprintf('Already averaged; taking principal eigenmmode to timefreq\n');
    YY = PEig(YY');
    YY = YY-YY(1);
    Y  = bert_singlechannel([YY';YY'],cfg,Freqs);
    
    
elseif ndims(YY)  > 2;
    % SVD based dim reduction over 3D (pre-trial averaged) data
    fprintf('NOT averaged; taking principal eigenmmode over 3D to timefreq\n');
    YY  = checkcheck(YY);
    YY  = squeeze(YY);
    YY  = YY'; 
    Y   = bert_singlechannel(YY,cfg,Freqs);
 
end

end


function Y = ica_engine(YY,cfg,Freqs)
% calls fastica for channel space reduction
%
%
% opts on bert_singlechan
% cfg.compute - 
%         'abs' -absolute value of the Hilbert trace (default) 
%         'pow' - Compute the power of the Hilbert trace
%         'phase' - COmpute the phase of the Hilbert trace     
%         'complex' keep the complex elements
%         'plv' - compute the PLV across


if  ndims(YY) == 2;
    % pricipal eigenmode over YY'
    %fprintf('Already averaged; doing ica\n');
    %YY = PEig(YY');
    %Y  = bert_singlechannel([YY';YY'],cfg,Freqs);
    
    
elseif ndims(YY)  > 2;
    % SVD based dim reduction over 3D (pre-trial averaged) data
    fprintf('Data not averaged; doing ica over sources\n');
    
    YY  = ica_check(YY);
    YY  = squeeze(YY);
    YY  = YY'; 
    Y   = bert_singlechannel(YY,cfg,Freqs);
 
end

end

function [NY] = ica_check(Y)
global ncc

if isempty(ncc)
    fprintf('finding number of ica for dataset\n');
    QQ = squeeze(mean(Y,3));
    [Cc,A,W] = fastica(QQ);
    ncc = size(Cc,1);
end

fprintf('decomposing channel space to %d components\n',size(ncc,1));
clear A W QQ

for i = 1:size(Y,3)
    QQ      = squeeze(Y(:,:,i));
    [C,A,W] = fastica(QQ);
    C       = C(1:ncc,:);
    
    NY(:,:,i) = C;
end

NY = squeeze(mean(Y,1));
% ( C'*pinv(W)' )'

end

function D = DoMEEG(ID)

fprintf('removing non neurophys channels\n');

% Find MEG & EEG channels
MEG = strfind(ID.chanlabels,'MEG');
MEG = find(~cellfun(@isempty,MEG));

EEG = strfind(ID.chanlabels,'EEG');
EEG = find(~cellfun(@isempty,EEG));

% Neurophys channels only
D = ID(sort(unique([EEG MEG])),:,:);

end

function [NY] = checkcheck(Y)
% reduce channel space to minimum eigenmodes needed to explain 90% variance
% per trial 
global n

lim = 10000;
jsz = ceil(size(Y,2) / lim);

if size(Y,2)*size(Y,3) < lim
    NY = Y;
else
    fprintf('reducing channel space\n');
    if isempty(n)
        QQ     = squeeze(mean(Y,3));
        [eQ,n] = PEig90(QQ');
    end
    
    for i = 1:size(Y,3)
        eQ        = Y(:,:,i);
        eQ        = eQ';
        eQ        = PEig(eQ,n);
        NY(:,:,i) = eQ';
    end
end   


end

function SS = GetSS(D)
% Convert from D(sens ,samps,trials)
%         to   D(verts,samps,trials)

inv     = D.inv{end};
js      = 1:length(inv.inverse.M);
scale   = inv.inverse.scale;
J       = inv.inverse.J;  
U       = inv.inverse.U; % spatial projector 
T       = inv.inverse.T; %
TT      = T*T';
M       = inv.inverse.M(js, :);
Ic      = inv.inverse.Ic;
It      = inv.inverse.It;
Np      = length(It);
trial   = D.condlist;
clabell = {};

for i = 1:numel(trial)
    c       = D.indtrial(trial{i}, 'GOOD');
    clabell = [clabell D.conditions(c)];
    
    Nt    = length(c);
    fprintf('reconstructing trial %d\n',i);
    
    for j = 1:Nt
        for k = 1:length(U)
            Y       = D(Ic{k},It,c(j));
            UY{k,1} = U{k}*Y*scale(k);
        end
        Y = spm_cat(UY);
        
        if j > 1
            MY{i} = MY{i} + M*Y;
        else
            MY{i} = M*Y;
        end
    end
end

MY = squeeze(inner(MY));
SS = permute(MY,[3 2 1]);

end

