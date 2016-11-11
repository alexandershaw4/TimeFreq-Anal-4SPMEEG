function [MMN,t,F,T,Y] = DoMMN(f,cfg)
% Called by timefreq_coupling
% Calls SPM_2TF
%
% Calculates MMN for roving paradigm
%
%

try cfg;             catch; cfg= []; end
if ~isfield(cfg,'MMN'); cfg.MMN = 0; end

for i = 1:length(f)
    D = spm_eeg_load(f{i});
    [Data,time,F{i},T{i},Y{i}] = MakeMMN(D,cfg);
    
    if cfg.MMN;
        MMN{i,:} = Data;
        T{i}     = time;
    end
end

if cfg.MMN;
    MMN = squeeze(cat(3,MMN{:}));
    t   = T{i};
    DoPlots(MMN,t);
else
    MMN = [];
    t   = T{i}{1};
end

end

function [MMN,time,F,T,Y] = MakeMMN(D,cfg)

try cfg.freqs;    catch; cfg.freqs = 1:.25:60; end
try cfg.sensors;  catch; cfg.sensors  = 'pca'; end
try cfg.baseline; catch; cfg.baseline = 1;     end

[F,T,Y] = SPM_2TF(D,cfg);  clc;

MMN{1}  = PEig(Y{1}'-Y{2}');
MMN{2}  = PEig(Y{2}'-Y{3}');
MMN{3}  = PEig(Y{3}'-Y{4}');
MMN{4}  = PEig(Y{4}'-Y{5}');
MMN{5}  = PEig(Y{5}'-Y{6}');
MMN{6}  = PEig(Y{6}'-Y{7}');
MMN{7}  = PEig(Y{7}'-Y{8}');
MMN{8}  = PEig(Y{8}'-Y{9}');
MMN{9}  = PEig(Y{9}'-Y{10}');

time = T{1};

end

function DoPlots(M,t)

%M  = cat(2,M{:});
M = innercell(M);

M = squeeze(mean(M,1));
M = M';

figure
mesh(M);
set(gca,'fontsize',18);
xlabel('Mismatch','fontsize',18);
ylabel('Samples','fontsize',18);
zlabel('Magnitude','fontsize',18);

figure
contourf(M);
set(gca,'fontsize',18);
xlabel('Mismatch','fontsize',18);
ylabel('Samples','fontsize',18);

figure
ribbon(M)
set(gca,'fontsize',18);
xlabel('Mismatch','fontsize',18);
ylabel('Samples','fontsize',18);

figure,
corr_mat_tight(M)


end
