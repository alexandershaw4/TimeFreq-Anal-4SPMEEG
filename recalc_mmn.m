function [y,t,maxx,meann] = recalc_mmn(cfc,T,B)

try t = cfc.t; catch t = cfc.T; end
try y = cfc.y; catch y = cfc.Y; end
try f = cfc.f; catch f = cfc.F; end
Q = @squeeze;

if iscell(y); y = inner(y); end
if iscell(f); f = f{1}; end
if iscell(t); t = t{1}; end

y = abs(y);

T = [findthenearest(T(1),t):findthenearest(T(2),t)];
B = [findthenearest(B(1),t):findthenearest(B(2),t)]; 

[NC,NSub,Nsamp,Nf] = size(y);

G = @(a,b)(ones([size(a,1),1])*mean(b,1));
f = @(b,a)(100 * (b -  a) ./ a);

for i = 1:(NC-1)
    for j = 1:NSub
        DEV  = Q(y(i,j,T,:));
        OFF1 = Q(y(i,j,B,:));
        
        STD  = Q(y(i+1,j,T,:));
        OFF2 = Q(y(i+1,j,B,:));
        
        % whole trial subraction
        deviant = DEV - G(DEV,OFF1);
        standrd = STD - G(STD,OFF2);
        
        MMN{i,j} = PEig(deviant-standrd);
        
        % max diff
        maxx(i,j)  = max(PEig(DEV'));
        meann(i,j) = mean(PEig(DEV'));
    end
end

y = MMN;
t = T;
        
