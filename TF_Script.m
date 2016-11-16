


cfg.freqs = 1:.25:60;
cfg.MMN   = 0;
cfg.type  = 'evoked';

f = dir('mTrials*.mat');

cfc = timefreq_coupling({f.name},cfg);