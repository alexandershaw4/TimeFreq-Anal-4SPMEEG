function plot_tf_roving(cfc)

y = cfc.Y;
t = cfc.T;
f = cfc.F;

if iscell(cfc.Y) && ndims(cfc.Y) == 2
    clear y
    y = inner(cfc.Y);
end

if iscell(t); t = t{1}; end
if iscell(f); f = f{1}; end

y = permute(y,[2 1 3 4]);

nt = size(y,1);
ha = tight_subplot(2,nt/2,[.01 .03],[.1 .01],[.01 .01]);
M  = max(spm_vec(mean(y(1,:,:,:),2)));

for i = 1:size(y,1)
    trial = squeeze( mean (y(i,:,:,:),2) );
    
    axes(ha(i)), plotbert(t,f,trial',M,['Condition ' num2str(i)]);
end
