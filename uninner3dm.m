function y = uninner3dm(x,orig)
% Reverse effect of inner3dm
%
% test it, a la:
%
% X = [{randn(4,3)} {randn(4,3)}];
% X(:,:,2) = X;
% X(2,:,:) = X; % X is a 3D cell comprising doubles
%
% mX = inner3dm(X); % mX is a double, 5D matrix
% 
% test = permute(mX,[4 5 1 2 3]);
% testcell = spm_unvec(test,X); % testcell is the same as X
%
%

dX = permute(x,[4 5 6 1 2 3]);
y  = spm_unvec(dX,orig);