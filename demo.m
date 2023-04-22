% Multiply-connected quasiconformal map
%
% Compute a landmark-matching quasi-conformal map with arclength 
% parameterized boundary condition from a multiply-connected surface 
% to another using the methods in [1] and [2] with modifications.
%
% Main function:
% map = multiply_connected_quasiconformal_map(...
%     v1,f1,bdy1_all,blm1_all,ilm1,v2,f2,bdy2_all,blm2_all,ilm2)
%
% If you use this code in your own work, please cite the following papers:
%
% [1] G. P. T. Choi, 
%     "Efficient Conformal Parameterization of Multiply-Connected Surfaces
%      Using Quasi-Conformal Theory".
%     Journal of Scientific Computing, 87(3), 70, 2021.
% 
% [2] G. P. T. Choi and L. Mahadevan, 
%     "Planar Morphometrics Using Teichmüller Maps".
%     Proceedings of the Royal Society A, 474(2217), 20170905, 2018. 
%
% Copyright (c) 2023, Gary P. T. Choi

addpath(genpath('code'));
addpath('data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D example

% load two planar shapes with the same topology
load('example_2D.mat');

figure;
patch('Faces',f1,'Vertices',v1,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on; 
plot(v1(cell2mat(blm1_all),1),v1(cell2mat(blm1_all),2),'ro','MarkerFaceColor','r');
plot(v1(ilm1,1),v1(ilm1,2),'bo','MarkerFaceColor','b');
title('Planar shape 1');

figure;
patch('Faces',f2,'Vertices',v2,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot(v2(cell2mat(blm2_all),1),v2(cell2mat(blm2_all),2),'mo','MarkerFaceColor','m');
plot(v2(ilm2,1),v2(ilm2,2),'yo','MarkerFaceColor','y');
title('Planar shape 2');

% compute a landmark-matching quasiconformal map
map = multiply_connected_quasiconformal_map(...
    v1,f1,bdy1_all,blm1_all,ilm1,v2,f2,bdy2_all,blm2_all,ilm2);

figure;
patch('Faces',f1,'Vertices',map,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot(map(cell2mat(blm1_all),1),map(cell2mat(blm1_all),2),'ro','MarkerFaceColor','r');
plot(map(ilm1,1),map(ilm1,2),'bo','MarkerFaceColor','b');
title('Planar mapping result');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D example

% load two 3D surfaces with the same topology
load('example_3D.mat');

figure;
patch('Faces',f1,'Vertices',v1,'FaceColor','none','EdgeColor','k');
axis equal tight off; 
hold on; 
plot3(v1(cell2mat(blm1_all),1),v1(cell2mat(blm1_all),2),v1(cell2mat(blm1_all),3),'ro','MarkerFaceColor','r');
plot3(v1(ilm1,1),v1(ilm1,2),v1(ilm1,3),'bo','MarkerFaceColor','b');
title('Surface 1');

figure;
patch('Faces',f2,'Vertices',v2,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot3(v2(cell2mat(blm2_all),1),v2(cell2mat(blm2_all),2),v2(cell2mat(blm2_all),3),'mo','MarkerFaceColor','m');
plot3(v2(ilm2,1),v2(ilm2,2),v2(ilm2,3),'yo','MarkerFaceColor','y');
title('Surface 2');

%% compute a conformal parameterization using [1]
map1 = poly_annulus_conformal_map(v1,f1);
map2 = poly_annulus_conformal_map(v2,f2);

figure;
patch('Faces',f1,'Vertices',map1,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot(map1(cell2mat(blm1_all),1),map1(cell2mat(blm1_all),2),'ro','MarkerFaceColor','r');
plot(map1(ilm1,1),map1(ilm1,2),'bo','MarkerFaceColor','b');
title('Conformal parameterization of Surface 1');

figure;
patch('Faces',f2,'Vertices',map2,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot(map2(cell2mat(blm2_all),1),map2(cell2mat(blm2_all),2),'mo','MarkerFaceColor','m');
plot(map2(ilm2,1),map2(ilm2,2),'yo','MarkerFaceColor','y');
title('Conformal parameterization of Surface 2');

%% compute a landmark-matching quasiconformal map
[map,map_planar] = multiply_connected_quasiconformal_map(...
    v1,f1,bdy1_all,blm1_all,ilm1,v2,f2,bdy2_all,blm2_all,ilm2);

figure;
patch('Faces',f1,'Vertices',map_planar,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot(map_planar(cell2mat(blm1_all),1),map_planar(cell2mat(blm1_all),2),'ro','MarkerFaceColor','r');
plot(map_planar(ilm1,1),map_planar(ilm1,2),'bo','MarkerFaceColor','b');
title('Planar mapping result');

figure;
patch('Faces',f1,'Vertices',map,'FaceColor','none','EdgeColor','k');
axis equal tight off
hold on;
plot3(map(cell2mat(blm1_all),1),map(cell2mat(blm1_all),2),map(cell2mat(blm1_all),3),'ro','MarkerFaceColor','r');
plot3(map(ilm1,1),map(ilm1,2),map(ilm1,3),'bo','MarkerFaceColor','b');
title('Surface mapping result');
