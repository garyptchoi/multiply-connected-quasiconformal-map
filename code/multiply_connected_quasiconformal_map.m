function [map, map_planar] = multiply_connected_quasiconformal_map(v1,f1,bdy1_all,blm1_all,ilm1,v2,f2,bdy2_all,blm2_all,ilm2)

% Compute a landmark-matching quasi-conformal map with arclength 
% parameterized boundary condition from a multiply-connected surface 
% to another using the methods in [1] and [2] with modifications.
%
% Input:
% v1: nv1 x 2 or 3 vertex coordinates of shape 1
% f1: nf1 x 3 triangulations of shape 1
% bdy1_all: a 1 x k cell array with each cell containing vertex indices of a boundary curve of shape 1
% blm1_all: a 1 x k cell array with each cell containing vertex indices of boundary landmarks on the k-th boundary curve of shape 1
% ilm1: indices of interior landmarks of shape 1
% v2: nv1 x 2 or 3 vertex coordinates of shape 2
% f2: nf1 x 3 triangulations of shape 2
% bdy2_all: a 1 x k cell array with each cell containing vertex indices of a boundary curve of shape 2
% blm2_all: a 1 x k cell array with each cell containing vertex indices of boundary landmarks on the k-th boundary curve of shape 2
% ilm2: indices of interior landmarks of shape 2
%
% Output:
% map: nv1 x 2 or 3 vertex coordinates of the resulting quasiconformal map
% map_planar: nv1 x 2 vertex coordinates of the quasiconformal map on the plane
% 
% If you use this code in your own work, please cite the following papers:
%
% [1] G. P. T. Choi, 
%     "Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory".
%     Journal of Scientific Computing, 87(3), 70, 2021.
% 
% [2] G. P. T. Choi and L. Mahadevan, 
%     "Planar Morphometrics Using Teichmüller Maps".
%     Proceedings of the Royal Society A, 474(2217), 20170905, 2018. 
%
% Copyright (c) 2023, Gary P. T. Choi

if size(v1,2) == 3
    % compute a conformal parameterization using [1]
    map1 = poly_annulus_conformal_map(v1,f1);
else
    map1 = v1;
end

if size(v2,2) == 3
    % compute a conformal parameterization using [1]
    map2 = poly_annulus_conformal_map(v2,f2);
else
    map2 = v2;
end

if length(ilm1) > size(ilm1,1)
    ilm1 = ilm1';
end
if length(ilm2) > size(ilm2,1)
    ilm2 = ilm2';
end
    
bdy1_id_all = [];
bdytarget_all = [];

for k = 1:length(bdy1_all)
    % handle the boundary one by one based on [2]
    
    bdy1 = bdy1_all{k};
    bdy2 = bdy2_all{k};
    blm1 = blm1_all{k};
    blm2 = blm2_all{k};
    
    bdyk_id = [];
    bdyk_target = [];
    
    if length(bdy1) > size(bdy1,1)
        bdy1 = bdy1';
    end
    if length(bdy2) > size(bdy2,1)
        bdy2 = bdy2';
    end
    if length(blm1) > size(blm1,1)
        blm1 = blm1';
    end
    if length(blm2) > size(blm2,1)
        blm2 = blm2';
    end
    
    if ~isempty(setdiff(blm1,bdy1))
        error('Some boundary landmarks do not lie on the correct boundary.');
    end
    
    if ~isempty(setdiff(blm2,bdy2))
        error('Some boundary landmarks do not lie on the correct boundary.');
    end
    
    for i = 1:length(blm1)
        % find the segment between two boundary landmarks
        if i < length(blm1)
            startid1 = blm1(i);
            endid1 = blm1(i+1);
            startid2 = blm2(i);
            endid2 = blm2(i+1);
        else
            startid1 = blm1(i);
            endid1 = blm1(1);
            startid2 = blm2(i);
            endid2 = blm2(1);
        end
        
        id11 = find(bdy1 == startid1);
        id12 = find(bdy1 == endid1);
        id21 = find(bdy2 == startid2);
        id22 = find(bdy2 == endid2);
        
        % check if the order is correct
        if id11 < id12
            segment1 = (id11:id12)';
            if ~isempty(intersect(bdy1(segment1(2:end-1)),blm1))
                % contain other landmarks, i.e. incorrect
                segment1 = [id12:length(bdy1), 1:id11]';
            end
        else
            segment1 = (id12:id11)';
            if ~isempty(intersect(bdy1(segment1(2:end-1)),blm1))
                % contain other landmarks, i.e. incorrect
                segment1 = [id11:length(bdy1), 1:id12]';
            end
        end
        
        if id21 < id22
            segment2 = (id21:id22)';
            if ~isempty(intersect(bdy2(segment2(2:end-1)),blm2))
                % contain other landmarks, i.e. incorrect
                segment2 = [id22:length(bdy2), 1:id21]';
            end
        else
            segment2 = (id22:id21)';
            if ~isempty(intersect(bdy2(segment2(2:end-1)),blm2))
                % contain other landmarks, i.e. incorrect
                segment2 = [id21:length(bdy2), 1:id22]';
            end
        end
        
        % find the accumulated length
        segmentlength1 = sqrt(sum((map1(bdy1(segment1(2:end)),:)-map1(bdy1(segment1(1:end-1)),:)).^2,2));
        
        % interpolate the boundary points on the target segment
        trange = [0;cumsum(segmentlength1)]/(sum(segmentlength1));
        trange(end) = 1; % for numerical stability
        [pt,~,~] = interparc(trange,map2(bdy2(segment2),1),map2(bdy2(segment2),2),'linear');
        
        bdyk_id = [bdyk_id; bdy1(segment1)];
        bdyk_target = [bdyk_target; pt];   
    end
    bdy1_id_all = [bdy1_id_all; bdyk_id];
    bdytarget_all = [bdytarget_all; bdyk_target];
end

%% Landmark matching quasiconformal map with the above boundary correspondence as in [1]
point = [bdy1_id_all;ilm1];
target = [bdytarget_all; map2(ilm2,:)];

mu = zeros(length(f1),1);
converge = 0;
iteration_count = 0;
while converge == 0
    iteration_count = iteration_count + 1;
    % compose the map with another quasi-conformal map to cancel the distortion
    A = generalized_laplacian(map1,f1,mu); 
    A(point,:) = 0;
    A(point,point) = diag(ones(length(point),1));

    bx = zeros(length(map1),1); 
    by = bx;
    bx(point) = target(:,1);
    by(point) = target(:,2);

    % solve the generalized Laplace equation
    map_x = A\bx;
    map_y = A\by;
    map_planar = [map_x,map_y];

    update_mu = beltrami_coefficient(map1,f1,map_planar);
    count = length(find(abs(update_mu) >= 1));
    display(['Overlap # = ',int2str(count)]);
    if count == 0 
        converge = 1;
    else
        update_mu(abs(update_mu)>1) = 0.9*update_mu(abs(update_mu)>1)./abs(update_mu(abs(update_mu)>1));
        update_mu = mean(abs(update_mu)).*(cos(angle(real(update_mu)+sqrt(-1)*imag(update_mu)))+...
            sqrt(-1)*sin(angle(real(update_mu)+sqrt(-1)*imag(update_mu))));
        mu = update_mu;
    end
    
    if iteration_count > 500
        warning('The algorithm fails to converge after 500 iterations. Terminated.');
        converge = 1;
    end
end

map = map_planar;

if size(v2,2) == 3
    % get the final surface mapping result using interpolation
    F1 = scatteredInterpolant(map2,v2(:,1),'linear');
    F2 = scatteredInterpolant(map2,v2(:,2),'linear');
    F3 = scatteredInterpolant(map2,v2(:,3),'linear');

    map = zeros(size(v1,1),3);
    map(:,1) = F1(map_planar(:,1),map_planar(:,2));
    map(:,2) = F2(map_planar(:,1),map_planar(:,2));
    map(:,3) = F3(map_planar(:,1),map_planar(:,2));
end