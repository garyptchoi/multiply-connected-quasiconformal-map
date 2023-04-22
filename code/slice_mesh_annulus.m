function [v_sliced,f_sliced] = slice_mesh_annulus(v,f,slice_path,bd)
% Slice a mesh with annulus topology to obtain a simply-connected open mesh
%
% Input:
% v: nv x 2 or 3 vertex coordinates
% f: nf x 3 triangulation (in anti-clockwise orientation)
% slice_path: np x 1 vertex indices of the slice path
% bd: the inner boundary
%
% Output:
% v_sliced: (nv+np) x 2 or 3 vertex coordinates of the sliced mesh
% f_sliced: nf x 3 triangulation of the sliced mesh
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, 
%     "Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory."
%     Journal of Scientific Computing, 87(3), 70, 2021.
%
% Copyright (c) 2021, Gary P. T. Choi
% https://math.mit.edu/~ptchoi

%%

nv_ori = length(v);
nf_ori = length(f);

% add one artificial vertex at the inner center
centroid = sum(v(bd,:),1)/length(bd);
v = [v; centroid];
f = [f; bd([2:end,1]), bd, (nv_ori+1)*ones(length(bd),1)];
slice_path = [nv_ori+1, slice_path];

nv = length(v);
nf = length(f);
np = length(slice_path);
v_sliced = [v; v(slice_path(2:np),:)];
e = [f(:,1),f(:,2); f(:,2),f(:,3); f(:,3),f(:,1)];
f_sliced = f;

%%
%         p1
% ----------------------
%         || 
%  update || keep
%         || 
% ----------------------
%         pn

for i = 2:np
    pid = slice_path(i);
    dir2 = find((e(:,2) == slice_path(i-1) & e(:,1) == slice_path(i)));
    fid = mod(dir2-1,nf)+1;
    id = find(f(fid,:) == pid);
    f_sliced(fid,id) = nv+i-1;
    
    %% also update the other faces on that side
    flag = 1;
    while flag
        %% continue
        switch id
            case 1   
                fid = find(e(:,1) == f(fid,1) & e(:,2) == f(fid,3));
            case 2
                fid = find(e(:,1) == f(fid,2) & e(:,2) == f(fid,1));
            case 3
                fid = find(e(:,1) == f(fid,3) & e(:,2) == f(fid,2));
        end
        fid = mod(fid-1,nf)+1;
        id = find(f(fid,:) == pid);
        f_sliced(fid,id) = nv+i-1;
        
        flag = (sum(ismember(f(fid,:),slice_path)) == 1);
    end
end

% remove the artificially added center and the faces
v_sliced(nv_ori+1,:) = [];
f_sliced = f_sliced(1:nf_ori,:);

% correct the indices
f_sliced(f_sliced >= nv_ori+2) = f_sliced(f_sliced >= nv_ori+2) - 1;
