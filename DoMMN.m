function [MMN,t] = DoMMN(f,cfg)


try cfg;             catch; cfg= []; end
if ~isfield(cfg,'MMN'); cfg.MMN = 0; end

for i = 1:length(f)
    Y = cfg.Y(i,:);
    T = cfg.T(i,:);
    [Data,time] = MakeMMN(Y,T);

    MMN{i,:} = Data;
    T{i}     = time;
    
end

MMN = squeeze(cat(3,MMN{:}));
t   = T{i};
DoPlots(MMN,t);

end



function [MMN,time] = MakeMMN(Y,T)

for i = 1:length(Y)-1
    MMN{i}  = PEig(abs(Y{i}')-abs(Y{i+1}'));
    time    = T{1};
end

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
