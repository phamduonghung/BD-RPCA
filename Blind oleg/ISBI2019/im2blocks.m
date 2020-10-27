function [z,K] = im2blocks(g,N,dN)

% IM2BLOCKS Partition to blocks
%   IM2BLOCKS(g,N,dN), with 'g' being a 2-D image/matrix, partitions the
%   former into N(1)-by-N(2) blocks, with a spatial overlap of dN(1) and
%   dN(2) pixels in the column and row directions, respectively. The re-
%   sulting blocks are returned as a 1-D cell array, with its structural
%   arrangement corresponding to a lexicographic (column-major) ordering
%   of the blocks.
%
%   If 'N' is a scalar, it is assumed that N(1)=N(2)=N. In a similar man-
%   ner, if 'dN' is a scalar, it is assumed that dN(1)=dN(2)=dN. When the
%   'dN' argument is ommitted, the function assumes a 50% overlap between
%   adjaccent blocks.
%
%   Instead of a 2-D array, 'g' can be defined to be a matrix/image size,
%   in which case the function returns the starting and ending indices of
%   the blocks in the form of a four-column matrix 'z'. In this case, the
%   k-th block of some "implied" matrix 'g' (of the specified size) would 
%   be defined as 'g(z(k,1):z(k,2),z(k,3):z(k,4))'.
%
%                       [z,K] = im2blocks(g,N,dN)
%   Input:
%                   g - either a 2-D array or an image/matrix size        
%                   N - block dimension
%                  dN - block overlap (default: dN = floor(N/2)) 
%   Output:
%                   z - either a cell array of the resulting blocks or a 
%                       matrix of block indices
%                   K - number of blocks in the column and row directions
%
%   See also GETBLOCKS, PARTIDOM
%
%   Written by Oleg Michailovich, August, 2017

cls={'numeric'};
ats={'integer','positive','size',[1 2]};
fname='im2blocks';

if isscalar(N)
    N=N([1 1]);
end
validateattributes(N,cls,ats,fname,'N')

if nargin>2
    if isscalar(dN)
        dN=dN([1 1]);
    end
    validateattributes(dN,cls,ats,fname,'dN')
    
    if any((N-dN)<1)
        error('Step size cannot exceed block size!')
    end
else
    dN=floor(N/2);
end

if ismatrix(g)
    if numel(g)==2
        validateattributes(g,cls,ats,fname,'g')
        n=g;
        c=true;
    else
        n=size(g);
        c=false;
    end
else
    error('The first input must be either a matrix or matrix dimensions!')
end

if any((n-N)<1)
    error('Block size cannot exceed image size!')
end

K=floor((n-N)./dN+1);
n0=floor((n-((K-1).*dN+N))/2);
I=[n0(1)+1 n0(1)+N(1) n0(2)+1 n0(2)+N(2)];
M=prod(K);
I=I(ones(M,1),:);

for k=1:M
    i1=mod(k-1,K(1));
    i2=(k-i1-1)/K(1);
    
    I(k,1:2)=I(k,1:2)+i1*dN(1);
    I(k,3:4)=I(k,3:4)+i2*dN(2);
end

if c
    z=I;
else
    z=cell(1,M);
    for k=1:M
        z{k}=g(I(k,1):I(k,2),I(k,3):I(k,4));
    end
end

end