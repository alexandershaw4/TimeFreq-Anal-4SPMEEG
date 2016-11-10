function [F,T,Y,FT] = SPM_2TF(D,cfg)
% TF for SPM MEEG Data either robust averaged or not
% 
% optionals cfg struct with optional
% cfg.sensors  = reduction method: pca, mean or sep
% cfg.resamp   = up/down sample by fixed val (+/- n)
% cfg.freqs    = frequency vector, eg. 1:.25:100
% cfg.baseline = do baselining [0/1]
% cfg.basetime = start / stop in secs, eg. [-.1 1]
% cfg.compute  = {'abs','pow','phase','complex','plv'}
%
% AS
try cfg; catch cfg = []; end

if isfield(cfg,'sensors'); meth = cfg.sensors; else meth = 'pca';    end
if isfield(cfg,'resamp');  rs   = cfg.resamp ; else rs   = 1;        end
if isfield(cfg,'freqs');   Freqs= cfg.freqs;   else Freqs= 1;        end

if isfield(cfg,'baseline');bl   = cfg.baseline;else bl   = 0;        end
if isfield(cfg,'basetime');bs   = cfg.basetime;else bs   = [-.1 .1]; end


% foi
if length(Freqs) == 1;
    if ~isempty(D.frequencies)
         Freqs = D.frequencies;
    else Freqs = 1:.2:100;
    end
end

% Data dimensions
S    = size(D);
ns   = S(1);
nc   = S(3);
nsam = S(2);

% check whether robust averaged
if nc > 25; 
    nc = length(D.condlist);
    for i = 1:(nc)
        ic{i} = D.indtrial(D.condlist{i});
    end
else
    ic = num2cell(1:10);
end


% FT
FT          = D.fttimelock;
dt          = 1/D.fsample;
FS          = 1/dt;        % sampling frequenct
t           = D.time;      % peristim time
cfg.fsample = FS*rs;       % adjusted sample rate
t           = linspace(t(1), t(end),rs*length(t));
smooth      = 4;

for i = 1:nc
    Yc = D(:,:,ic{i});
    cfg.sampletimes  = t;
    cfg.kill_evagram = 0;
    
    if bl;  cfg.baseline     = 'subtract';
            cfg.basetimes    = bs;
    else    cfg.baseline     = 'none';
    end
        
    switch meth
        case 'pca';
            YY   = HighResMeanFilt(full(Yc),rs,smooth);
            TF_Y = eig_engine(YY,cfg,Freqs);

            Y{i} = double(TF_Y.evagram);
            F{i} = TF_Y.freqs;
            T{i} = TF_Y.sampletimes;
            
        case 'mean';
            YY   = HighResMeanFilt(full(Yc),rs,smooth);
            YY   = mean(YY);
            TF_Y = bert_singlechannel([YY;YY],cfg,Freqs);
            
            Y{i} = double(TF_Y.evagram);
            F{i} = TF_Y.freqs;
            T{i} = TF_Y.sampletimes;
        
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
                Y{i,j} = double(TF_Y.evagram);
                
                F{i,j} = TF_Y.freqs;
                T{i,j} = TF_Y.sampletimes;
                
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
    Y  = bert_singlechannel([YY';YY'],cfg,Freqs);
    
    
elseif ndims(YY)  > 2;
    % SVD based dim reduction over 3D (pre-trial averaged) data
    fprintf('NOT averaged; taking principal eigenmmode over 3D to timefreq\n');
    X  = PEig(VecRetainDim(permute(YY,[2 1 3])),1:6);
    YY = mean(X);
    Y  = bert_singlechannel([YY;YY],cfg,Freqs);
 
end

end

