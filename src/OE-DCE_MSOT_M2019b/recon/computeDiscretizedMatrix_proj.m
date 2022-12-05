function A_mat_p = computeDiscretizedMatrix_proj( x_pt, y_pt, R_pt, theta, image_width, n, lt )

nn = n*n ; % number of columns of the matrix
n_angles = size( x_pt, 1 ) ; % number of points of the curve
Dxy = image_width / ( n - 1 ) ; % sampling distance in x and y

x_pt_unrot = x_pt * cos( theta ) - y_pt * sin( theta ) ; % horizontal position of the points of the curve in the original grid (not rotated)
y_pt_unrot = x_pt * sin( theta ) + y_pt * cos( theta ) ; % vertical position of the points of the curve in the original grid (not rotated)

x_aux_1 = [ x_pt_unrot; zeros( 1, lt ) ] ;
y_aux_1 = [ y_pt_unrot; zeros( 1, lt ) ] ;
x_aux_2 = [ zeros( 1, lt ) ; x_pt_unrot ] ;
y_aux_2 = [ zeros( 1, lt ) ; y_pt_unrot ] ;
dist_aux = sqrt( ( x_aux_1 - x_aux_2 ).^2 + ( y_aux_1 - y_aux_2 ).^2 ) ;
clear x_aux_1 x_aux_2 y_aux_1 y_aux_2
d_pt = dist_aux( 2 : n_angles, : ) ; % length of the segments of the curve
clear dist_aux

vec_int = (1/2) * ( [ d_pt; zeros( 1, lt ) ] + [ zeros( 1, lt ); d_pt ] ) ./ R_pt ; % vector for calculating the integral
clear d_pt

x_pt_pos_aux = ( x_pt_unrot + ( image_width/2 ) ) / Dxy + 1 ; % horizontal position of the points of the curve in normalized coordinates
clear x_pt_unrot
y_pt_pos_aux = ( y_pt_unrot + ( image_width/2 ) ) /Dxy + 1 ; % vertical position of the points of the curve in normalized coordinates
clear y_pt_unrot

x_pt_pos_bef = floor( x_pt_pos_aux ) ; % horizontal position of the point of the grid at the left of the point (normalized coordinates)
x_pt_pos_aft = floor( x_pt_pos_aux+1 ) ; % horizontal position of the point of the grid at the right of the point (normalized coordinates)
x_pt_dif_bef = x_pt_pos_aux - x_pt_pos_bef ;
clear x_pt_pos_aux
y_pt_pos_bef = floor( y_pt_pos_aux ) ; % vertical position of the point of the grid below of the point (normalized coordinates)
y_pt_pos_aft = floor( y_pt_pos_aux+1 ) ; % vertical position of the point of the grid above of the point (normalized coordinates)
y_pt_dif_bef = y_pt_pos_aux - y_pt_pos_bef ;
clear y_pt_pos_aux

pos_sq_1x = x_pt_pos_bef; pos_sq_1y = y_pt_pos_bef; % position of the first point of the square
pos_sq_2x = x_pt_pos_aft; pos_sq_2y = y_pt_pos_bef; % position of the second point of the square
pos_sq_3x = x_pt_pos_bef; pos_sq_3y = y_pt_pos_aft; % position of the third point of the square
pos_sq_4x = x_pt_pos_aft; pos_sq_4y = y_pt_pos_aft; % position of the fourth point of the square
in_pos_sq_1 = (pos_sq_1x>0)&(pos_sq_1x<=n)&(pos_sq_1y>0)&(pos_sq_1y<=n); % boolean determining if the first point of the square is inside the grid
in_pos_sq_2 = (pos_sq_2x>0)&(pos_sq_2x<=n)&(pos_sq_2y>0)&(pos_sq_2y<=n); % boolean determining if the first point of the square is inside the grid
in_pos_sq_3 = (pos_sq_3x>0)&(pos_sq_3x<=n)&(pos_sq_3y>0)&(pos_sq_3y<=n); % boolean determining if the first point of the square is inside the grid
in_pos_sq_4 = (pos_sq_4x>0)&(pos_sq_4x<=n)&(pos_sq_4y>0)&(pos_sq_4y<=n); % boolean determining if the first point of the square is inside the grid
Pos_sq_1_t = n*(pos_sq_1x-1)+pos_sq_1y; % one dimensional position of the first points of the squares in the grid
Pos_sq_2_t = n*(pos_sq_2x-1)+pos_sq_2y; % one dimensional position of the first points of the squares in the grid
Pos_sq_3_t = n*(pos_sq_3x-1)+pos_sq_3y; % one dimensional position of the first points of the squares in the grid
Pos_sq_4_t = n*(pos_sq_4x-1)+pos_sq_4y; % one dimensional position of the first points of the squares in the grid
clear pos_sq_1x pos_sq_2x pos_sq_3x pos_sq_4x
Pos_sq_1_t_vec = reshape(Pos_sq_1_t,1,n_angles*lt); % Pos_triang_1_t in vector form
Pos_sq_2_t_vec = reshape(Pos_sq_2_t,1,n_angles*lt); % Pos_triang_1_t in vector form
Pos_sq_3_t_vec = reshape(Pos_sq_3_t,1,n_angles*lt); % Pos_triang_1_t in vector form
Pos_sq_4_t_vec = reshape(Pos_sq_4_t,1,n_angles*lt); % Pos_triang_1_t in vector form
clear Pos_sq_1_t Pos_sq_2_t Pos_sq_3_t Pos_sq_4_t

weight_sq_1 = (1-x_pt_dif_bef).*(1-y_pt_dif_bef).*vec_int; % weight of the first point of the triangle
weight_sq_2 = (x_pt_dif_bef).*(1-y_pt_dif_bef).*vec_int; % weight of the second point of the triangle
weight_sq_3 = (1-x_pt_dif_bef).*(y_pt_dif_bef).*vec_int; % weight of the third point of the triangle
weight_sq_4 = (x_pt_dif_bef).*(y_pt_dif_bef).*vec_int; % weight of the fourth point of the triangle
weight_sq_1_t_vec = reshape(weight_sq_1,1,n_angles*lt); % weight_sq_1 in vector form
weight_sq_2_t_vec = reshape(weight_sq_2,1,n_angles*lt); % weight_sq_1 in vector form
weight_sq_3_t_vec = reshape(weight_sq_3,1,n_angles*lt); % weight_sq_1 in vector form
weight_sq_4_t_vec = reshape(weight_sq_4,1,n_angles*lt); % weight_sq_1 in vector form
clear weight_sq_1 weight_sq_2 weight_sq_3 weight_sq_4

Row_Matrix = ( linspace( 1, lt, lt )' * ones( 1, n_angles ) )' ; % rows of the sparse matrix
Row_Matrix_vec = reshape( Row_Matrix, 1, n_angles * lt ) ; % rows of the sparse matrix in vector form
clear Row_Matrix

A_mat_p = sparse( [ Pos_sq_1_t_vec( in_pos_sq_1 ) Pos_sq_2_t_vec( in_pos_sq_2 ) Pos_sq_3_t_vec( in_pos_sq_3 ) Pos_sq_4_t_vec( in_pos_sq_4 ) ],...
    [ Row_Matrix_vec( in_pos_sq_1 ) Row_Matrix_vec( in_pos_sq_2 ) Row_Matrix_vec( in_pos_sq_3 ) Row_Matrix_vec( in_pos_sq_4 ) ],...
    [ weight_sq_1_t_vec( in_pos_sq_1 ) weight_sq_2_t_vec( in_pos_sq_2 ) weight_sq_3_t_vec( in_pos_sq_3 ) weight_sq_4_t_vec( in_pos_sq_4 ) ], nn, lt ) ;



























