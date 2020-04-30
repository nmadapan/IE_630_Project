%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the data for the IMRT Problem. 
% 
% Assumptions:
%   1. There is one tumor, three organs at risk and normal tissue. 
%   2. The region is sub-divided into cuboidal voxels. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%% CHANGE THESE VARIABLES
DIRECTION = 'mxy'; % Collimator direction
WITH_AREA = true; % If true, it will compute Joules/sec, otherwise, Joules/m^2/sec

%% Initialization: Collimator
CY = 1; % CZ = CY % Half side of the the collimator. 
CN = 2*CY / 0.5; % No. of points on each axes (Y and Z) of collimator. 

%% Initialization: Ellipse
EX = 10; % Half x-axis of the Ellipse
EY = 5; % Half y-axis of the Ellipse
EZ = 1; % Half z-axis of the Ellipse
% dd = 0.4;
% dxe = dd; dye = dd; dze = dd;
% Nx = ceil(2*EX/dxe); if(rem(Nx, 2) == 0); Nx = Nx + 1; end
% Ny = ceil(2*EY/dye); if(rem(Ny, 2) == 0); Ny = Ny + 1; end
% Nz = ceil(2*EZ/dze); if(rem(Nz, 2) == 0); Nz = Nz + 1; end
N = 20;
Nx = N + 1;
Ny = N + 1;
Nz = N + 1;
dxe = 2*EX/Nx; % Range along x-axis is 20
dye = 2*EY/Ny; % Range along y-axis is 10
dze = 2*EZ/Nz; % Range along z-axis is 4
% TODO: Need to play with this threshold. 
% d_thresh = norm([dxe, dye, dze]);
% d_thresh = min([dxe, dye, dze]);
% d_thresh = sqrt(2)/2.81; % Based on the collimator resolution. 

%% Initialization: Properties
F = 10;
ALPH_T = 0.7;
ALPH_N = 0.1;
ALPH_O = 0.15;
COUNT = 0;

assert(any(strcmp({'x', 'y', 'xy', 'mxy'}, DIRECTION)), ...
    'Error! DIRECTION should be "x", "y" or "xy"');

if(strcmp(DIRECTION, 'x'))
    UNIT_DISTANCE = dxe;
end
if(strcmp(DIRECTION, 'y'))
    UNIT_DISTANCE = dye;
end
if(strcmp(DIRECTION, 'xy') || strcmp(DIRECTION, 'mxy'))
    UNIT_DISTANCE = norm([dxe, dye]);    
end
UNIT_AREA = UNIT_DISTANCE * dze;

%% Initialization: Create symbols and ellipse equations
tic;
syms x y z E T O1 O2 O3;
E = x^2/(EX^2) + y^2/(EY^2) - 1; % Elliptical region which is the body. 
T = x^2 + y^2 + z^2 - 1; % Tumor
O1 = (x-3)^2 + 4*y^2 + z^2/4 - 1; % Organ at risk 
O2 = 4*(x-3)^2 + 4*(y-3)^2 + 4*z^2 - 1; % Organ at risk
O3 = 4*x^2 + 9*(y-2)^2 + 16*z^2 - 1; % Organ at risk
O4 = (x+3)^2 + 4*y^2 + z^2/4 - 1;
O5 = 4*(x+3)^2 + 4*(y-3)^2 + 4*z^2 - 1;
O6 = 4*x^2 + 9*(y+2)^2 + 16*z^2 - 1;
O7 = 4*(x+3)^2 + 4*(y+3)^2 + 4*z^2 - 1;
O8 = 4*(x-3)^2 + 4*(y+3)^2 + 4*z^2 - 1;

fprintf('Creating sybmbols: %.02f seconds\n', toc);

%% Plot tumor, critical organs, normal tissue
% Note that the plotting is done on x-axis. 
tic
figure;
hold on;
grid on;
ezplot(E, [-1*EX, EX, -1 * EY, EY])
ezplot(subs(T, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O1, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O2, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O3, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O4, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O5, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O6, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O7, z, 0), [-1*EX, EX, -1 * EY, EY])
ezplot(subs(O8, z, 0), [-1*EX, EX, -1 * EY, EY])
xlim([-1*(EX+0.5), EX+0.5])
ylim([-1*(EY+0.5), EY+0.5])
fprintf('Plotting: %.02f seconds\n', toc);

%% Create the BODY MESH
% Make the precision in the ellipse is smaller than that of the collimator.
% So one cell in the collimator will hit only one voxel. 
tic;
xe_pts = linspace(-1*EX, EX, Nx);
ye_pts = linspace(-1*EY, EY, Ny);
ze_pts = linspace(-1*EZ, EZ, Nz);
[Xe, Ye, Ze] = meshgrid(xe_pts, ye_pts, ze_pts);
Xe = Xe(:); Ye = Ye(:); Ze = Ze(:);
P = [Xe, Ye, Ze];
% Plot the mesh
scatter(Xe, Ye, 'g')
fprintf('Creating ellipse mesh: %.02f seconds\n', toc);

%% Determine the necessary flags
% Find the points that fall inside the ellipse
tic;
in_eps_flags = double(subs(E, {x, y}, {P(:,1), P(:,2)})) <= 0;
% Plot the points that fall inside the ellipse
scatter(P(in_eps_flags, 1), P(in_eps_flags, 2), 'm')

% Find the points that are in OAR
in_oar_flags = false;
o1 = double(subs(O1, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o1 <= 0;
o2 = double(subs(O2, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o2 <= 0;
o3 = double(subs(O3, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o3 <= 0;
o4 = double(subs(O4, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o4 <= 0;
o5 = double(subs(O5, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o5 <= 0;
o6 = double(subs(O6, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o6 <= 0;
o7 = double(subs(O7, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o7 <= 0;
o8 = double(subs(O8, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_oar_flags = in_oar_flags | o8 <= 0;

clear o1 o2 o3
% Plot the points that fall inside the OAR
scatter(P(in_oar_flags, 1), P(in_oar_flags, 2), 'k')

% Find the points inside the tumor
t1 = double(subs(T, {x, y, z}, {P(:,1), P(:,2), P(:, 3)}));
in_tumor_flags = t1 <= 0;
clear t1
% Plot the points that fall inside the OAR
scatter(P(in_tumor_flags, 1), P(in_tumor_flags, 2), 'r')

% Find the points inside the normal tissue
in_normal_flags = (in_eps_flags & ~in_tumor_flags) & ~in_oar_flags;
temp = P(:, 3) == 0; % Find the points on z = 0 only. And make a 2D plot. 
Q = P(temp, :);
scatter(Q(in_normal_flags(temp), 1), Q(in_normal_flags(temp), 2), 'b')
clear temp Q
fprintf('Computing necessary flags: %.02f seconds\n', toc);

%% Discretize the collimator.
tic;
% It is in the Y-Z plane. If not, we will rotate the points back. 
zc_pts = linspace(-1*CY, CY, CN+1); 
yc_pts = linspace(-1*CY, CY, CN+1); 
xc_pts = EX * ones(size(yc_pts));

% Collimator plane equation: ax + by + cz + d = 0: x = 10
if(strcmp(DIRECTION, 'x'))
    col_plane = [1, 0, 0, -1 * EX]; 
end
if(strcmp(DIRECTION, 'xy'))
    col_plane = [1, 1, 0, -1 * EX]; 
end
if(strcmp(DIRECTION, 'mxy'))
    col_plane = [1, -1, 0, -1 * EX];
end
if(strcmp(DIRECTION, 'y'))
    col_plane = [0, 1, 0, -1 * EX]; 
end
fprintf('Discretize collimator: %.02f seconds\n', toc);

%% Loop over and find the matrix
tic;
% TODO: Need to tune this parameter. 
d_thresh = min([dxe, dye, dze]);
% Project all voxels in the body on to the collimator. 
[Q, ~] = point_to_plane(col_plane, P); 

all([Q, ones(size(Q,1),1)] * col_plane' == 0)

% If y: rotate the x,y in Q by -90 degrees. 
% If xy: rotate the x,y in Q by -45 degrees. 
if(strcmp(DIRECTION, 'x'))
    Q_prime = Q;
end
if(strcmp(DIRECTION, 'xy')) % -45 % x cos() - y sin(), x sin() + y cos()
    % [sqrt(1/2)*(x+y), sqrt(1/2)*(-x+y), z]
    tx = sqrt(1/2) * (Q(:, 1) + Q(:, 2));
    ty = sqrt(1/2) * (-1*Q(:, 1) + Q(:, 2));
    Q_prime = [tx, ty, Q(:, 3)];
end
if(strcmp(DIRECTION, 'mxy')) % -135 % x cos() - y sin(), x sin() + y cos()
    theta = -135;
    tx = cosd(theta)*Q(:, 1) - sind(theta)*Q(:, 2);
    ty = sind(theta)*Q(:, 1) + cosd(theta)*Q(:, 2);
    Q_prime = [tx, ty, Q(:, 3)];
end
if(strcmp(DIRECTION, 'y')) % -90 % x cos() - y sin(), x sin() + y cos()
    % [y, -x, z]
    Q_prime = [Q(:, 2), -1*Q(:, 1), Q(:, 3)];
end

% Find if the projected points fall inside the collimator.
[in_col_flags, indices] = where_in_collimator(Q_prime, yc_pts, zc_pts);

% Mark the voxels that are already visited. 
visited_flags = false * ones(size(P, 1), 1);
visited_flags(~in_eps_flags) = true; % Ignore the voxels that are outside the ellipse. 

% Find out the intensity at each voxel. 
intensity_list = zeros(size(P, 1), 1);
% Intensity of voxels outside the ellipse is marked -1. 
intensity_list(~in_eps_flags) = -1; 

% Loop over all points
for idx = 1 : size(P, 1)
   if(visited_flags(idx)); continue; end; % Points that are already visited are ignored. 
   if(~in_col_flags(idx)) % If projected point doesnt fall on collimator, intensity is zero. 
       visited_flags(idx) = true;
       intensity_list(idx) = 0; % This variable is already defaulted to zero.
       continue; 
   end
   
   v1 = P(idx, :); % Point at consideration
   v2 = Q(idx, :); % Projected point on the collimator
   if(isequal(v1, v2))
       % Intensity is exactly equal to initial intensity. 
       visited_flags(idx) = true;
       intensity_list(idx) = F;
       continue
   end
   
   % Find the voxels that are close to the 3D line joined by v1 and v2.
   dists = point_to_line_distance(P, v1, v2);
   dists(~in_eps_flags) = inf;
   flags = dists <= d_thresh;
   
   if(strcmp(DIRECTION, 'xy') || strcmp(DIRECTION, 'x'))
       % flags = flags & P(:, 1) >= P(idx,1); % Remove the voxels that are to the left of v1.
       % Sort the voxels w.r.t x-axis. 
       [new_P, I] = sortrows(P(flags, :), 1);
   end
   if(strcmp(DIRECTION, 'y') || strcmp(DIRECTION, 'mxy'))
       % flags = flags & P(:, 2) >= P(idx,2); % Remove the voxels that are to the left of v1.
       % Sort the voxels w.r.t y-axis. 
       [new_P, I] = sortrows(P(flags, :), 2); 
   end

   % Unsort to sort, Use I directly. Sort to unsort, use J instead. 
   temp = sortrows([I, (1:numel(I))'], 1);
   J = temp(:, 2);
   
   % Sanity check. v1 should be the first point in the new set of points. 
   % We are not removing the points that are to the left or below so, this
   % assert makes no sense. 
   % assert(isequal(new_P(1, :), v1), 'Error: Initial point mismatch')
   
   % Find new set of tumor flags. 
   new_tumor_flags = in_tumor_flags(flags);
   new_tumor_flags = new_tumor_flags(I);
   
   % Find new set of normal flags. 
   new_normal_flags = in_normal_flags(flags);
   new_normal_flags = new_normal_flags(I);
   
   % Find new set of OAR flags. 
   new_oar_flags = in_oar_flags(flags);
   new_oar_flags = new_oar_flags(I);
   
   % Combine these three flags. 1 -> normal, 2 -> tumor and 3 -> OAR.
   combined_info = zeros(size(new_tumor_flags));
   combined_info(new_normal_flags == 1) = 1;
   combined_info(new_tumor_flags == 1) = 2;
   combined_info(new_oar_flags == 1) = 3;
   
   % Run a loop backward and upadte the intensities of the respective
   % voxels. 
   temp_intensity_list = zeros(size(combined_info, 1), 1);
   for vid = size(combined_info, 1):-1:1
       COUNT = COUNT + 1;
       if(vid == size(combined_info, 1))
           temp_intensity_list(vid) = F;
           continue;
       end
       if(combined_info(vid) == 1); alph = ALPH_N; end
       if(combined_info(vid) == 2); alph = ALPH_T; end
       if(combined_info(vid) == 3); alph = ALPH_O; end
       temp_intensity_list(vid) = temp_intensity_list(vid + 1) * exp(-1 * alph * UNIT_DISTANCE);
   end
   
   % Unsort the temp_intensity_list w.r.t to orignal array P. Use J. 
   intensity_list(flags) = temp_intensity_list(J); % Sorted to unsorted.
   
   visited_flags(flags) = true;
end
fprintf('Percentage visited: %.02f %% \n', 100*COUNT/size(P, 1))
fprintf('Rest: %.02f seconds\n', toc);

if(WITH_AREA)
    intensity_list = UNIT_AREA * intensity_list;
end

dir_name = ['data_', num2str(N)];
try
   mkdir(dir_name)
catch exp 
end

mat_fname = fullfile(dir_name, ['data_', DIRECTION]);
% assignin('base', [DIRECTION, '_indices'], indices);
save(mat_fname, 'indices', 'intensity_list', 'in_normal_flags', ...
    'in_eps_flags', 'in_tumor_flags', 'in_oar_flags', 'yc_pts', 'zc_pts', ...
    'xc_pts', 'xe_pts', 'ye_pts', 'ze_pts', 'in_col_flags', 'P', 'dxe', 'dye', 'dze')

flags = P(:, 3) == 0;
figure;
imagesc(reshape(intensity_list(flags), numel(xe_pts), numel(ye_pts)))
colorbar;

real_vol = (pi * EX * EY * 2 * EZ);
est_vol = sum(in_eps_flags) * dxe * dye * dze;
diff_vol = est_vol - real_vol;

fprintf('Volume of ellipse: %f\n', real_vol);
fprintf('Est. Volume of ellipse: %f\n', est_vol);
fprintf('Diff: %f\n', diff_vol);

flags = (flags & in_eps_flags);
real_area = (pi * EX * EY);
est_area = sum(flags) * dxe * dye;
diff_area = est_area - real_area;

fprintf('Area of ellipse: %f\n', real_area);
fprintf('Est. Area of ellipse: %f\n', est_area);
fprintf('Diff: %f\n', diff_area);

dV = dxe * dye * dze; 
fprintf('### Tumor sanity check ###\n')
fprintf('Real volume: %.02f\n', pi * 4 / 3)
fprintf('Est. volume: %.02f\n', sum(in_tumor_flags) * dV)

function [in_flags, indices] = where_in_collimator(V, yr, zr)
    %%%%%%%%%%%%%%%%%%%%
    % Description:
    %   Determine the voxels in V, that are effected by the collimator. 
    %   This function checks if the voxels fall with in the range of
    %   collimator y and z axes. 
    % Input:
    %   V - matrix of size N x 3. Each row is an [x, y, z] vector. 
    %   yr - grid of points on collimator along y-direction. It MUST be a
    %   sorted array. 
    %   zr - grid of points on collimator along z-direction. It MUST be a
    %   sorted array. 
    %%%%%%%%%%%%%%%%%%%%
    assert(issorted(yr), 'Error! yr should be a sorted array. ')
    assert(issorted(zr), 'Error! zr should be a sorted array. ')
    assert(size(V, 2) == 3, 'Error! V should have three columns. ')
    
    in_flags = in_collimator(V, yr, zr);
    
    % Initialize
    y = V(:, 2); z = V(:, 3); % Column vectors
    num_y = numel(yr); num_z = numel(zr); % Scalars
    
    % Find the y indices of all points on the collimator. 
    mat_yr = repmat(yr, numel(y), 1); % Repmat yr to a matrix
    [~, Ky] = sort([mat_yr, y], 2); % Find where y fits into yr
    [row, col] = find(Ky == num_y + 1); % The index of num_y+1 is the desired index of y. 
    temp = sortrows([row, col], 1);
    final_y_ids = temp(:, 2);
    % Substract 1 when values in y are exactly equal to values in yr. 
    final_y_ids = final_y_ids - sum(y == yr, 2); 
    
    % Find the z indices of all points on the collimator. 
    mat_zr = repmat(zr, numel(z), 1); % Repmat zr to a matrix
    [~, Kz] = sort([mat_zr, z], 2); % Find where z fits into zr
    [row, col] = find(Kz == num_z + 1); % The index of num_z+1 is the desired index of z. 
    % find function randomly shuffles the row ids, so they should be
    %   sorted first. 
    temp = sortrows([row, col], 1);
    final_z_ids = temp(:, 2);
    % Substract 1 when values in z are exactly equal to values in zr. 
    final_z_ids = final_z_ids - sum(z == zr, 2); % sum(z == zr, 2) is vector of ones/zeros.
    
    % Convert (i, j) into an scalar index computed in a column first format. 
    indices = num_y * (final_z_ids - 1) + final_y_ids;
    
    % Finally, set the indices of the points outside the collimator region
    % to zero. Such points might have valid indices due to the way the
    % indices are computed. 
    indices(~in_flags) = 0;
    
    % Sanity check on the indices. 
    assert(min(indices) == 0 && max(indices) == numel(yr)*numel(zr), 'Error! Invalid indices.')
end

function flag = in_collimator(V, yr, zr)
    %%%%%%%%%%%%%%%%%%%%
    % Description:
    %   Determine the voxels in V, that are effected by the collimator. 
    %   This function checks if the voxels fall with in the range of
    %   collimator y and z axes. 
    % Input:
    %   V - matrix of size N x 3. Each row is an [x, y, z] vector. 
    %   yr - grid of points on collimator along y-direction. 
    %   zr - grid of points on collimator along z-direction.
    %%%%%%%%%%%%%%%%%%%%
    y = V(:, 2); z = V(:, 3);
    min_yr = min(yr); max_yr = max(yr);
    min_zr = min(zr); max_zr = max(zr);
    % Initialize
    flag = true * ones(size(y));
    % Check if y falls between y min and y max of collimator
    flag = flag & (y <= max_yr) & (y >= min_yr);
    % Check if z falls between z min and z max of collimator
    flag = flag & (z <= max_zr) & (z >= min_zr);
end