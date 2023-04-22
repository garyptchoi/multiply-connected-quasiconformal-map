function map = poly_annulus_conformal_map(v,f)
% Compute a conformal map of a multiply-connected open mesh onto a 2D 
% circle domain using the poly-annulus conformal map (PACM) method in [1].
%
% Input:
% v: nv x 2 or 3 vertex coordinates of a multiply-connected surface
% f: nf x 3 triangulations of a multiply-connected surface
% 
% Output:
% map: nv x 2 vertex coordinates of the conformal parameterization
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
if nb < 2
    error('The input mesh is not a multiply-connected surface.');
end

bd_list = 1:nb;

% find outer boundary
bd_outer_id = 1;
bd_outer = bd{bd_outer_id};
bd_outer_length = sum(sqrt(sum((v(bd_outer,:) - v(bd_outer([2:end,1]),:)).^2,2)));
for i = 2:nb
    bd_next = bd{i};
    bd_next_length = sum(sqrt(sum((v(bd_next,:) - v(bd_next([2:end,1]),:)).^2,2)));
    if bd_next_length > bd_outer_length
        bd_outer_id = i;
        bd_outer_length = bd_next_length;
    end
end

%% Initial annulus map without the quasi-conformal composition step
map = v;
temp = setdiff(bd_list, bd_outer_id);
bd_scanned = [bd_outer_id,temp(1)];
v_filled = map;
f_filled = f;
count = 1;
for i = setdiff(bd_list,[bd_outer_id,temp(1)])
    % fill all holes except one
    centroid = sum(map(bd{i},:),1)/length(bd{i});
    if ~inpolygon(centroid(:,1),centroid(:,2),map(bd{i},1),map(bd{i},2))
        [cx,cy,~] = find_inner_circle(map(bd{i},1),map(bd{i},2));
        centroid(:,1:2) = [cx,cy];
    end
    v_filled = [v_filled; centroid];
    f_filled = [f_filled; bd{i}([2:end,1]), bd{i}, (length(map) + count)*ones(length(bd{i}),1)];
    count = count + 1;
end
%  compute annulus map
map_filled = annulus_map(v_filled,f_filled,bd{bd_outer_id}, bd{temp(1)});

%% Handle the holes one by one
while ~isempty(setdiff(bd_list, bd_scanned))
    v_new = map_filled(1:length(map),:);

    temp = setdiff(bd_list, bd_scanned);
    if ~isempty(temp)
        % pick a new hole to exclude from filling
        bd_scanned = [bd_scanned,temp(1)];
        v_new_filled = v_new;
        f_new_filled = f;
        count = 1;
        for i = setdiff(bd_list,[bd_outer_id,temp(1)])
            % fill all holes except one
            centroid = sum(v_new(bd{i},:),1)/length(bd{i});
            if ~inpolygon(centroid(:,1),centroid(:,2),v_new(bd{i},1),v_new(bd{i},2))
                [cx,cy,~] = find_inner_circle(v_new(bd{i},1),v_new(bd{i},2));
                centroid = [cx,cy];
            end
            v_new_filled = [v_new_filled; centroid];
            f_new_filled = [f_new_filled; bd{i}([2:end,1]), bd{i}, (length(v_new) + count)*ones(length(bd{i}),1)];
            count = count + 1;

        end

        map_filled = annulus_map(v_new_filled,f_new_filled,bd{bd_outer_id},bd{temp(1)});

    end
end

%% Mobius transformation on the disk to reduce area distortion
% (Choi et al., SIAM J. Imaging Sci. 2020)
map = mobius_area_correction_disk(v,f,map_filled(1:length(v),:));

%% Enforce exact circularity
map_z = complex(map(:,1),map(:,2));
for i = setdiff(1:nb,bd_outer_id)
    [cx,cy,~] = find_inner_circle(real(map_z(bd{i})),imag(map_z(bd{i})));
    c = complex(cx,cy);
    r = min(abs(map_z(bd{i}) - c)); 
    map_z(bd{i}) = (map_z(bd{i}) - c)./abs(map_z(bd{i}) - c)*r + c;
end
map = [real(map_z), imag(map_z)];

%% quasi-conformal composition for improving the conformality
% compute the Beltrami coefficient
mu = beltrami_coefficient(map, f, v); 
fixed = [];
for i = 1:nb
    fixed = [fixed; bd{i}];
end

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

%% optimal rotation to align with the input mesh
corners = zeros(4,1);
[~,corners(1)] = min(abs(v(:,1) - min(v(:,1))) + abs(v(:,2) - min(v(:,2))));
[~,corners(2)] = min(abs(v(:,1) - max(v(:,1))) + abs(v(:,2) - min(v(:,2))));
[~,corners(3)] = min(abs(v(:,1) - max(v(:,1))) + abs(v(:,2) - max(v(:,2))));
[~,corners(4)] = min(abs(v(:,1) - min(v(:,1))) + abs(v(:,2) - max(v(:,2))));
[U, ~, ~] = Kabsch([map(corners,1:2),zeros(4,1)]', [v(corners,1:2),zeros(4,1)]');
map = (U(1:2,1:2)*map')';

end
