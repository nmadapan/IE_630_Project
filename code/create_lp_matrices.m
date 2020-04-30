close all;
clear;
clc;

global N_coll intensity_list indices

%% DATA X
load('data_x.mat')
N_coll = numel(yc_pts) * numel(zc_pts);

% Plot the heatmap
z = P(:, 3);
flags = z == 0;
X = reshape(intensity_list(flags), numel(xe_pts), numel(ye_pts));
figure;
imagesc(X)
colorbar;

% D matrix for normal tissue
D_N1 = return_dmat(in_normal_flags, in_col_flags);
D_N1 = [D_N1, zeros(size(D_N1, 1), 2 * N_coll)];
Dt_N1 = sum(D_N1, 1);
% D matrix for tumor
D_T1 = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
D_O1 = return_dmat(in_oar_flags, in_col_flags);
D_O1 = [D_O1, zeros(size(D_O1, 1), 2 * N_coll)];
Dt_O1 = sum(D_O1, 1);

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
D_N2 = return_dmat(in_normal_flags, in_col_flags);
D_N2 = [zeros(size(D_N2, 1), N_coll), D_N2, zeros(size(D_N2, 1), N_coll)];
Dt_N2 = sum(D_N2, 1);
% D matrix for tumor
D_T2 = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
D_O2 = return_dmat(in_oar_flags, in_col_flags);
D_O2 = [zeros(size(D_O2, 1), N_coll), D_O2, zeros(size(D_O2, 1), N_coll)];
Dt_O2 = sum(D_O2, 1);

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
D_N3 = return_dmat(in_normal_flags, in_col_flags);
D_N3 = [zeros(size(D_N3, 1), 2 * N_coll), D_N3];
Dt_N3 = sum(D_N3, 1);
% D matrix for tumor
D_T3 = return_dmat(in_tumor_flags, in_col_flags);
% D matrix for organ at risk
D_O3 = return_dmat(in_oar_flags, in_col_flags);
D_O3 = [zeros(size(D_O3, 1), 2 * N_coll), D_O3];
Dt_O3 = sum(D_O3, 1);

D_T = [D_T1, zeros(size(D_T1, 1), 2 * N_coll); ...
       zeros(size(D_T2, 1), N_coll), D_T2, zeros(size(D_T2, 1), N_coll); ...
       zeros(size(D_T3, 1), 2*N_coll), D_T3];
Dt_T = sum(D_T, 1);

C_N = Dt_N1 + Dt_N2 + Dt_N3;
C_O = Dt_O1 + Dt_O2 + Dt_O3;
   
save('mlop_data', 'D_N1', 'D_N2', 'D_N3', 'D_O1', 'D_O2', ...
     'D_O3', 'D_T', 'Dt_N1', 'Dt_N2', 'Dt_N3', 'Dt_O1', ...
     'Dt_O2', 'Dt_O3', 'Dt_T', 'C_N', 'C_O');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = return_dmat(flag1, flag2)
    global N_coll intensity_list indices
    flags = flag1 & flag2;
    K = sum(flags);
    D = zeros(K, N_coll);
    D((indices(flags)-1)*K + (1:K)') = intensity_list(flags);
end