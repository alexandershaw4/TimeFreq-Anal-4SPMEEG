function y = HighResMeanFilt(mIn,varargin)
% upscales matrices and filters using a mean-window
% 'mIn' is the matrix
% varargin{1} is the upscale ratio n:1 (e.g. enter '2' = 2:1, double)
% varargin{2} is the SmoothWindow on the filter
% note: use negative values for downsampling. Also works on 3D data
% AS 2015

if ndims(mIn) > 2; % handle high dimentions
    mOut = size(mIn);
    mIn  = VecRetainDim(mIn,1);
end
    

x     = mIn;
try n = varargin{1}; catch n = 2; end
try m = varargin{2}; catch m = 1; end

% resample 1st dim
if n > 0; nn = 1; else nn = (0 - n) + 1; n = 1; end
dX = resample(x,n,nn,2);

% upsample 2nd dim
for i = 1:size(dX,1)
    ddX(i,:) = resample(dX(i,:),n,nn,2);
end

% smooth (mean of window, func below)
y = smoothmat(ddX,m,m);

    try mOut; % reinstate orginal dimentions 
        y = spm_unvec(y, ones([mOut(1)*(n*nn) mOut(2)*(n*nn) mOut(3)]));
    end
end

function matrixOut = smoothmat(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr; 
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

A = isnan(matrixIn);
matrixIn(A) = 0;

nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;

end