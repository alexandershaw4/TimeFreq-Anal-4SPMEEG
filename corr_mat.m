function r = corr_mat(x)
% Input full double matrix, x, to plot correlation matrix and return
% corresponding r and p values in structure. Diagonals are histograms with
% std normal fit line (help histfit). 
% Non diagonals are scatters (p<.05 = red)
%
% AS2016

[~,xi] = size(x);
npl    = xi*(xi-1);
v      = 1:xi;
n      = 0;

for i  = 1:xi
    X  = x(:,i);            %# dependent
    yi = v(~ismember(v,i)); %# ind of indep

    for j = 1:xi
        if i == j; Y = X;
        else       Y = x(:,j);
        end;       n = n + 1;
        
        [R,p] = corr(X,Y);
        
        subplot(xi,xi,n),
        if i == j; histfit(X); [R,p]=deal(0);
        else       scatter(X,Y); lsline;
                   if p < .05; scatter(X,Y,'r'); lsline; end
                   title(sprintf('r = %d,\n p = %d',R,p));
        end
        
        % save correl vals
        r.R(i,j) = R;
        r.p(i,j) = p;
        
    end
end
            
        
        
        
        
        
        
        
        

