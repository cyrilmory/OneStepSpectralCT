function [] = PlayIteratesTransposed(iterates,mat_inf, mat_sup)
[ObjectSize, nmat, niterates] = size(iterates);
ObjectSize = sqrt(ObjectSize);

% Normalize each material to uint8
mat_rescaled = uint8(zeros(size(iterates)));
for mat=1:nmat
    mat_rescaled(:,mat,:) = uint8((iterates(:,mat,:) - mat_inf(mat)) * 255 / (mat_sup(mat) - mat_inf(mat)));
end

% Play sequences side by side
tmp=permute(reshape(mat_rescaled, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);
% tmp=flip(tmp,1);
% tmp=flip(tmp,2);

implay(reshape(tmp, [ObjectSize, ObjectSize*nmat, niterates]));

% % Write it
% writeVideoResult(double(reshape(mat_rescaled, [ObjectSize, ObjectSize * nmat, niterates])), 'WeidingerNoRegul256_2000iter.avi', 0,255);
end