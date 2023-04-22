function map = annulus_map(v,f,bd_outer,bd_inner)
% Compute a map of a triangle mesh with annulus topology onto an annulus
%
% Input:
% v: nv x 2 or 3 vertex coordinates of the input surface
% f: nf x 3 triangulation of the input surface
%
% Output:
% map: nv x 2 vertex coordinates of the annulus map
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, 
%     "Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory"
%     Journal of Scientific Computing, 87(3), 70, 2021.
%
% Copyright (c) 2021, Gary P. T. Choi
% https://math.mit.edu/~ptchoi

%% slice the mesh 
e = [f(:,[1,2]);f(:,[2,3]);f(:,[3 1])];
nv = length(v);
ne = length(e);

% arbitrarily choose a vertex at the inner boundary
inner_bd_vid = bd_inner(1);

% find the closest vertex at the outer boundary
[~,outer_bd_vid_temp] = min(sum((v(bd_outer,:) - repmat(v(inner_bd_vid,:),length(bd_outer),1)).^2,2));
outer_bd_vid = bd_outer(outer_bd_vid_temp);
% find the shortest path between the two vertices
G = sparse(e(:,1),e(:,2),ones(ne,1),nv,nv);
G = (G + G');
[~, slice_path, ~] = graphshortestpath(G,inner_bd_vid,outer_bd_vid);

% ensure no other vertices in slice_path are on the boundaries
[~,id1,~] = intersect(slice_path,bd_inner);
[~,id2,~] = intersect(slice_path,bd_outer);
slice_path = slice_path(max(id1):min(id2));

% slice mesh along the shortest path 
[v_sliced,f_sliced] = slice_mesh_annulus(v,f,slice_path,bd_inner);

%% Map to a rectangle and then an annulus
corner = [slice_path(1); slice_path(length(slice_path)); nv+length(slice_path); nv+1];

rect_sliced = rectangular_conformal_map_periodic(v_sliced,f_sliced,corner);
rect = rect_sliced(1:nv,:);

% Apply an exponential map to map the rectangle to an annulus
map_z = exp(2*pi*(complex(rect(:,1),rect(:,2))-max(rect(:,1))));
map = [real(map_z),imag(map_z)];
