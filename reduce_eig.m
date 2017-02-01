function y = reduce_eig(x,d)

[~, n] = PEig90(x);
ncyc = 0;

while n > d
    ncyc = ncyc + 1;
    x = HighResMeanFilt(x,1,2);
    [~, n] = PEig90(x);
    y = x;
end

fprintf('Finished after %d cycles\n',ncyc);