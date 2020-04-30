close all;
clear;
clc;

global N_coll intensity_list indices

%% DATA X
load('data_x.mat')
N_coll = numel(yc_pts) * numel(zc_pts);
dV = dxe * dye * dze;

% Plot the heatmap
z = P(:, 3);
flags = z == 0;
X = reshape(intensity_list(flags), numel(xe_pts), numel(ye_pts));
figure;
imagesc(X)
colorbar;

% D matrix for normal tissue
[D_N1, ~] = return_dmat(in_normal_flags, in_col_flags);
D_N1 = [D_N1, zeros(size(D_N1, 1), 2 * N_coll)];
% D matrix for tumor
[D_T1, T1_ids] = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
[D_O1, O1_ids] = return_dmat(in_oar_flags, in_col_flags);
D_O1 = [D_O1, zeros(size(D_O1, 1), 2 * N_coll)];

%% DATA X-Y
load('data_xy.mat')
N_coll = numel(yc_pts) * numel(zc_pts);

% Plot the heatmap
z = P(:, 3);
flags = z == 0;
X = reshape(intensity_list(flags), numel(xe_pts), numel(ye_pts));
figure;
imagesc(X)
colorbar;

% D matrix for normal tissue
[D_N2, ~] = return_dmat(in_normal_flags, in_col_flags);
D_N2 = [zeros(size(D_N2, 1), N_coll), D_N2, zeros(size(D_N2, 1), N_coll)];
% D matrix for tumor
[D_T2, T2_ids] = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
[D_O2, ~] = return_dmat(in_oar_flags, in_col_flags);
D_O2 = [zeros(size(D_O2, 1), N_coll), D_O2, zeros(size(D_O2, 1), N_coll)];

%% DATA Y

load('data_y.mat')
N_coll = numel(yc_pts) * numel(zc_pts);

% Plot the heatmap
z = P(:, 3);
flags = z == 0;
X = reshape(intensity_list(flags), numel(xe_pts), numel(ye_pts));
figure;
imagesc(X)
colorbar;

% D matrix for normal tissue
[D_N3, ~] = return_dmat(in_normal_flags, in_col_flags);
D_N3 = [zeros(size(D_N3, 1), 2 * N_coll), D_N3];
% D matrix for tumor
[D_T3, T3_ids] = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
[D_O3, ~] = return_dmat(in_oar_flags, in_col_flags);
D_O3 = [zeros(size(D_O3, 1), 2 * N_coll), D_O3];

%% DATA MXY (Minus X Plus Y)

load('data_mxy.mat')
N_coll = numel(yc_pts) * numel(zc_pts);

% Plot the heatmap
z = P(:, 3);
flags = z == 0;
X = reshape(intensity_list(flags), numel(xe_pts), numel(ye_pts));
figure;
imagesc(X)
colorbar;

% D matrix for normal tissue
[D_N4, ~] = return_dmat(in_normal_flags, in_col_flags);
D_N4 = [zeros(size(D_N4, 1), 2 * N_coll), D_N4];
% D matrix for tumor
[D_T4, T4_ids] = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
[D_O4, ~] = return_dmat(in_oar_flags, in_col_flags);
D_O4 = [zeros(size(D_O4, 1), 2 * N_coll), D_O4];

%%
D_T = [D_T1, D_T2, D_T3, D_T4];
 
save('mlop_data_voxel_wise', 'D_N1', 'D_N2', 'D_N3', 'D_N4', 'D_O1', 'D_O2', ...
     'D_O3', 'D_O4', 'D_T', 'dV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D, nz_ids] = return_dmat(flag1, flag2)
    global N_coll intensity_list indices
    flags = flag1 & flag2;
    K = sum(flags);
    D = zeros(K, N_coll);
    D((indices(flags)-1)*K + (1:K)') = intensity_list(flags);
    nz_ids = find(flags);
end
