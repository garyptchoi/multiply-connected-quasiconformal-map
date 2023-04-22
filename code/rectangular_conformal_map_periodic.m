function map = rectangular_conformal_map_periodic(v,f,corner)
% Compute the rectangular conformal mapping using the fast method in [1] 
% with the top and bottom boundaries enforced to be periodic in x as 
% described in [2].
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% corner: 4 x 1 vertex indices for the four corners of the rectangle, with anti-clockwise orientation
% 4 - 3
% |   |
% 1 - 2
% 
% Output:
% map: nv x 2 vertex coordinates of the rectangular conformal parameterization
% 
% If you use this code in your own work, please cite the following papers:
% [1] T. W. Meng, G. P.-T. Choi and L. M. Lui, 
%     "TEMPO: Feature-Endowed Teichmüller Extremal Mappings of Point Clouds."
%     SIAM Journal on Imaging Sciences, 9(4), pp. 1922-1962, 2016.
%
% [2] G. P. T. Choi, 
%     "Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory."
%     Journal of Scientific Computing, 87(3), 70, 2021.
%
% Copyright (c) 2021, Gary P. T. Choi
% https://math.mit.edu/~ptchoi

%% Initial setup

nv = length(v);

if size(v,1) < size(v,2)
    v = v';
end
if size(f,1) < size(f,2)
    f = f';
end
if size(v,2) == 2
    v = [v, zeros(nv,1)];
end

bd = meshboundaries(f);
bdy_index = bd{1};

% rearrange the boundary indices
id = find(bdy_index == corner(1));
bdy_index = bdy_index([id:end,1:id-1]);

id1 = 1;
id2 = find(bdy_index == corner(2));
id3 = find(bdy_index == corner(3));
id4 = find(bdy_index == corner(4));

if id2 > id3
    % the boundary orientation is wrong, meaning the input f has wrong orientation
    % we correct the orientation and start over
    warning('The input triangulations are with clockwise orientation!');
    f = fliplr(f);
    bd = meshboundaries(f);
    bdy_index = bd{1};
    id = find(bdy_index == corner(1));
    bdy_index = bdy_index([id:end,1:id-1]);
    id1 = 1;
    id2 = find(bdy_index == corner(2));
    id3 = find(bdy_index == corner(3));
    id4 = find(bdy_index == corner(4));
end

%% Step 1: Mapping the input mesh onto the unit disk

bdy_length = sqrt((v(bdy_index,1) - v(bdy_index([2:end,1]),1)).^2 + ...
            (v(bdy_index,2) - v(bdy_index([2:end,1]),2)).^2 + ...
            (v(bdy_index,3) - v(bdy_index([2:end,1]),3)).^2);
partial_edge_sum = zeros(length(bdy_length),1);

% arc-length parameterization boundary constraint
for i = 2:length(bdy_length)
    for j = 1:i-1
        partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length);
bdy = exp(theta*1i);

% disk harmonic map
M = cotangent_laplacian(v,f);
[mrow,mcol,mval] = find(M(bdy_index,:));
M = M - sparse(bdy_index(mrow),mcol,mval,nv, nv) + ...
        sparse(bdy_index,bdy_index,ones(length(bdy_index),1),nv, nv);
c = zeros(nv,1); 
c(bdy_index) = bdy;
z = M\c;
disk = [real(z),imag(z)]; 

if sum(sum(isnan(disk))) ~= 0
    % use tutte embedding instead
    disk = tutte_map(v,f,bdy_index,bdy); 
end

%% Step 2: Mapping the unit disk to the unit square

% compute the generalized Laplacian
mu = beltrami_coefficient(disk,f,v);
Ax = generalized_laplacian(disk,f,mu); 
Ay = Ax;

% set the boundary constraints
bx = zeros(nv,1); 
by = bx;

bottom = bdy_index(id1:id2);
right = bdy_index(id2:id3);
top = bdy_index(id3:id4);
left = bdy_index([id4:end,id1]);

% top and bottom must have the same x-coordinates
Ax(bottom,:) = 0;
Ax = Ax + sparse([bottom;bottom],[bottom;flipud(top)],[ones(length(bottom),1);-ones(length(top),1)],nv,nv);

% top and bottom y-coordinates must differ by 1
Ay([top;bottom],:) = 0;
Ay([top;bottom],[top;bottom]) = diag(ones(length([top;bottom]),1));
by(top) = 1;

% left must have x = 0, right must have x = 1
Ax([left;right],:) = 0;
Ax([left;right],[left;right]) = diag(ones(length([left;right]),1));
bx(right) = 1;

% solve the generalized Laplace equation
square_x = Ax\bx;
square_y = Ay\by;

%% Step 3: Optimize the width of the square to achieve a conformal map

w_opt = fminbnd(@(w) sum(abs(beltrami_coefficient([w*square_x,square_y],f,v)).^2),0,10);
map = [w_opt*square_x, square_y];

end