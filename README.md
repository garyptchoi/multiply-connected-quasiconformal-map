# Multiply-Connected Quasiconformal Map

<img src = "https://github.com/garyptchoi/multiply-connected-quasiconformal-map/blob/main/cover.png" height="360" />

This code computes a landmark-matching quasi-conformal map with arclength parameterized boundary condition from a multiply-connected surface to another using the methods in [1] and [2] with modifications.

Any comments and suggestions are welcome. 

If you use this code in your own work, please cite the following papers:

[1] G. P. T. Choi, 
    "[Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory.](https://doi.org/10.1007/s10915-021-01479-y)"
    Journal of Scientific Computing, 87(3), 70, 2021.

[2] G. P. T. Choi and L. Mahadevan, 
    "[Planar Morphometrics Using Teichmuller Maps.](https://doi.org/10.1098/rspa.2017.0905)"
    Proceedings of the Royal Society A, 474(2217), 20170905, 2018. 

Copyright (c) 2023, Gary Pui-Tung Choi

===============================================================

Usage:
* `[map,map_planar] = multiply_connected_quasiconformal_map(v1,f1,bdy1_all,blm1_all,ilm1,v2,f2,bdy2_all,blm2_all,ilm2)`

Input:
* `v1`: nv1 x 2 or 3 vertex coordinates of surface 1
* `f1`: nf1 x 3 triangulations of surface 1
* `bdy1_all`: a 1 x k cell array with each cell containing vertex indices of a boundary curve of surface 1
* `blm1_all`: a 1 x k cell array with each cell containing 1 x b_k vertex indices of boundary landmarks on the k-th boundary curve of surface 1
* `ilm1`: 1 x i vertex indices of interior landmarks of surface 1
* `v2`: nv2 x 2 or 3 vertex coordinates of surface 2
* `f2`: nf2 x 3 triangulations of surface 2
* `bdy2_all`: a 1 x k cell array with each cell containing vertex indices of a boundary curve of surface 2
* `blm2_all`: a 1 x k cell array with each cell containing 1 x b_k vertex indices of boundary landmarks on the k-th boundary curve of surface 2
* `ilm2`: 1 x i vertex indices of interior landmarks of surface 2

Output:
* `map`: nv1 x 2 or 3 vertex coordinates of the resulting quasiconformal map
* `map_planar`: nv1 x 2 vertex coordinates of the quasiconformal map on the plane

================================================================

Remarks:
* The `meshboundaries` function can be used for finding the boundaries of the two input surfaces. However, please make sure that the boundaries and the landmarks are in correct order.

================================================================

Dependencies:
* [Poly-Annulus Conformal Map](https://github.com/garyptchoi/poly-annulus-conformal-map)
* [find_inner_circle.m](https://www.mathworks.com/matlabcentral/fileexchange/32543-maximum-inscribed-circle-using-voronoi-diagram)
* [Kabsch.m](https://www.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm)
* [interparc.m](https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc)
