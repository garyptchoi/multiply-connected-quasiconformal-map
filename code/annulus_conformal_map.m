function map = annulus_conformal_map(v,f)
% Compute a conformal map of a triangle mesh with annulus topology onto 
% a 2D annulus domain using the annulus conformal map (ACM) method in [1].
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

bd = meshboundaries(f); % get all boundaries
nb = length(bd); % number of boundaries
if nb ~= 2
    error('The input mesh is not topologically equivalent to annulus!');
end

% determine which boundary is the outer one
bd1 = bd{1};
bd2 = bd{2};
if size(v,2) == 2
    % use the Pick's formula
    area1 = abs(sum(v(bd1,1).*v(bd1([2:length(bd1),1]),2)) - ...
        sum(v(bd1([2:length(bd1),1]),1).*v(bd1,2)));
    area2 = abs(sum(v(bd2,1).*v(bd2([2:length(bd2),1]),2)) - ...
        sum(v(bd2([2:length(bd2),1]),1).*v(bd2,2)));
    if area1 > area2
        outer_bd_id = 1;
        inner_bd_id = 2;
    else
        outer_bd_id = 2;
        inner_bd_id = 1;
    end        
else
    % pick the longer boundary as the outer boundary
    bd1_length = sum(sqrt(sum((v(bd1,:) - v(bd1([2:end,1]),:)).^2,2)));
    bd2_length = sum(sqrt(sum((v(bd2,:) - v(bd2([2:end,1]),:)).^2,2)));
    if bd1_length > bd2_length
        outer_bd_id = 1;
        inner_bd_id = 2;
    else
        outer_bd_id = 2;
        inner_bd_id = 1;
    end
end

%% slice the mesh 
e = [f(:,[1,2]);f(:,[2,3]);f(:,[3 1])];
nv = length(v);
ne = length(e);

% arbitrarily choose a vertex at the inner boundary
inner_bd_vid = bd{inner_bd_id}(1);

% find the closest vertex at the outer boundary
[~,outer_bd_vid_temp] = min(sum((v(bd{outer_bd_id},:) - repmat(v(inner_bd_vid,:),length(bd{outer_bd_id}),1)).^2,2));
outer_bd_vid = bd{outer_bd_id}(outer_bd_vid_temp);
% find the shortest path between the two vertices
G = sparse(e(:,1),e(:,2),ones(ne,1),nv,nv);
G = (G + G');
[~, slice_path, ~] = graphshortestpath(G,inner_bd_vid,outer_bd_vid);

% ensure no other vertices in slice_path are on the boundaries
[~,id1,~] = intersect(slice_path,bd{inner_bd_id});
[~,id2,~] = intersect(slice_path,bd{outer_bd_id});
slice_path = slice_path(max(id1):min(id2));

% slice mesh along the shortest path 
[v_sliced,f_sliced] = slice_mesh_annulus(v,f,slice_path,bd{inner_bd_id});

%% Map to a rectangle and then an annulus
corner = [slice_path(1); slice_path(length(slice_path)); nv+length(slice_path); nv+1];

rect_sliced = rectangular_conformal_map_periodic(v_sliced,f_sliced,corner);
rect = rect_sliced(1:nv,:);

% Apply an exponential map to map the rectangle to an annulus
map_z = exp(2*pi*(complex(rect(:,1),rect(:,2))-max(rect(:,1))));
map = [real(map_z),imag(map_z)];

%% quasi-conformal composition for improving the conformality
% compute the Beltrami coefficient
mu = beltrami_coefficient(map, f, v); 
fixed = [bd1;bd2];

% compose the map with another quasi-conformal map to cancel the distortion
A = generalized_laplacian(map,f,mu); 
A(fixed,:) = 0;
A(fixed,fixed) = diag(ones(length(fixed),1));

bx = zeros(length(v),1); 
by = bx;
bx(fixed) = map(fixed,1);
by(fixed) = map(fixed,2);

% solve the generalized Laplace equation
map_x = A\bx;
map_y = A\by;
map = [map_x,map_y];
