function C = mtimes3d(A,B)
% 3d matrix multiplication
% efficiently computes:
% for k=1:size(A,3)
%       C(:,:,k) = A(:,:,k) * B(:,:,k)
% end

if size(A,2)~=size(B,1) || (all([size(A,3) size(B,3)]>1) && size(B,3)~=size(A,3))
    error('matices size mismatch');
end

if all( [size(A,1) size(A,2) size(B,2)] < 10 )        
% implement by multiplying A with each column vector of B
    if size(B,2) > size(A,1)
        % swap A,B and calculate C = (B'*A')'
        tmp = A;
        A = permute(conj(B),[2 1 3]);
        B = permute(conj(tmp),[2 1 3]);
        swapped = 1;
    else
        swapped = 0;
    end
    
    nPage = max(size(A,3),size(B,3));
    nRow = size(A,1);
    nCol =  size(B,2);
    
    C = zeros(nRow,nCol,nPage);
    for colInd = 1:nCol
        W = permute(B(:,colInd,:),[2 1 3]);
        C(:,colInd,:) = sum( bsxfun(@times,A,W) , 2);
    end
    
    if swapped
        C = permute(conj(C),[2 1 3]);
    end


else 
    % implement by multiplying each page of A with each page of B
    nPage = max(size(A,3),size(B,3));
    nRow = size(A,1);
    nCol =  size(B,2);
    C = zeros(nRow,nCol,nPage);
    if size(B,3)==1, B = repmat(B,[1 1 nPage]);end
    if size(A,3)==1, A = repmat(A,[1 1 nPage]);end
    for pageInd = 1:nPage
        C(:,:,pageInd) = A(:,:,pageInd)*B(:,:,pageInd);
    end
end