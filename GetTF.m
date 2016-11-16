function ER = GetTF(cfc,T,F)
% CFC if a struc returned by SPM_2TF or do_tf_script
% T = time of interest  e.g. [.1 .3]
% F = freqs of interest e.g. [1 30]
% Gets the samples from subspace time and frequency data [e.g. the ERP]
% 
% 

y  = cfc.y;
t  = cfc.t;
f  = cfc.f;
s  = f(2)-f(1);

dt = [findthenearest(T(1),t):findthenearest(T(2),t)];
df = [findthenearest(F(1),f):findthenearest(F(2),f)];
Q  = @squeeze;

[nc,Nsub,nsamp,nf] = size(y);

for i = 1:nc
    for j = 1:Nsub
        temp = Q(y(i,j,dt',df));
        
        temp = Q(mean(temp,1));
        
        ER(i,j,:) = temp;
        F         = f(df);
        %T         = t(dt);
        
    end
end

ER.ER = ER;
ER.t  = T;
ER.f  = F;
ER.Q  = Q;