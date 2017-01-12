function bert_tf(D)

D              = spm_eeg_load(D);
cfg.fsample    = D.fsample;

nc             = size(D,1);
nt             = length(D.condlist);

fprintf('there are %d channels\n',nc);
fprintf('there are %d conditions\n',nt);

for i = 1:nc
    for j = 1:nt
        M = squeeze(D(i,:,D.indtrial(D.condlist{j})))';
        fprintf('there are %d trials for condition %d\n',size(M,1),j);

        cfg.compute    = 'pow';
        cfg.central    = 'median';
        cfg.start_time = 0;
        cfg.endtime    = .35;
        cfg.venumber   = 1;

        timefreq{i,j} = bert_singlechannel(M,cfg,5:0.2:60,[-2 0]);
    end
end

save('TimeFreq','timefreq');