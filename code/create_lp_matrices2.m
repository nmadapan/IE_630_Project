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
K = 5;
D_N1 = return_dmat(in_normal_flags, in_col_flags, K);
D_N1 = [D_N1, zeros(size(D_N1, 1), 2 * K)];
Dt_N1 = sum(D_N1, 1);
Vt_N1 = sum(D_N1~=0) * dV;
% D matrix for tumor
D_T1 = return_dmat(in_tumor_flags, in_col_flags, K);
V_T1 = sum(D_T1~=0) * dV;
% D matrix for organ at risk
D_O1 = return_dmat(in_oar_flags, in_col_flags, K);
D_O1 = [D_O1, zeros(size(D_O1, 1), 2 * K)];
Dt_O1 = sum(D_O1, 1);
Vt_O1 = sum(D_O1~=0) * dV;

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
D_N2 = return_dmat(in_normal_flags, in_col_flags, K);
D_N2 = [zeros(size(D_N2, 1), K), D_N2, zeros(size(D_N2, 1), K)];
Dt_N2 = sum(D_N2, 1);
Vt_N2 = sum(D_N2~=0) * dV;
% D matrix for tumor
D_T2 = return_dmat(in_tumor_flags, in_col_flags, K);
V_T2 = sum(D_T2~=0) * dV;
% D matrix for organ at risk
D_O2 = return_dmat(in_oar_flags, in_col_flags, K);
D_O2 = [zeros(size(D_O2, 1), K), D_O2, zeros(size(D_O2, 1), K)];
Dt_O2 = sum(D_O2, 1);
Vt_O2 = sum(D_O2~=0) * dV;

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
D_N3 = return_dmat(in_normal_flags, in_col_flags, K);
D_N3 = [zeros(size(D_N3, 1), 2 * K), D_N3];
Dt_N3 = sum(D_N3, 1);
Vt_N3 = sum(D_N3~=0) * dV;
% D matrix for tumor
D_T3 = return_dmat(in_tumor_flags, in_col_flags, K);
V_T3 = sum(D_T3~=0) * dV;
% D matrix for organ at risk
D_O3 = return_dmat(in_oar_flags, in_col_flags, K);
D_O3 = [zeros(size(D_O3, 1), 2 * K), D_O3];
Dt_O3 = sum(D_O3, 1);
Vt_O3 = sum(D_O3~=0) * dV;

D_T = [D_T1, zeros(size(D_T1, 1), 2 * K); ...
       zeros(size(D_T2, 1), K), D_T2, zeros(size(D_T2, 1), K); ...
       zeros(size(D_T3, 1), 2*K), D_T3];
V_T = [V_T1, zeros(size(V_T1, 1), 2 * K); ...
       zeros(size(V_T2, 1), K), V_T2, zeros(size(V_T2, 1), K); ...
       zeros(size(V_T3, 1), 2*K), V_T3];
Dt_T = sum(D_T, 1);
Vt_T = sum(V_T, 1);

C_N = Dt_N1 + Dt_N2 + Dt_N3;
C_O = Dt_O1 + Dt_O2 + Dt_O3;
C_T = Dt_T;

V_N = Vt_N1 + Vt_N2 + Vt_N3;
V_O = Vt_O1 + Vt_O2 + Vt_O3;
V_T = Vt_T;

% fprintf('### Tumor sanity check ###\n')
% fprintf('Real volume: %.02f\n', pi * 4 / 3)
% fprintf('Est. volume: %.02f\n', sum(in_tumor_flags) * dV)

save('mlop_data2', 'D_N1', 'D_N2', 'D_N3', 'D_O1', 'D_O2', ...
     'D_O3', 'D_T', 'Dt_N1', 'Dt_N2', 'Dt_N3', 'Dt_O1', ...
     'Dt_O2', 'Dt_O3', 'Dt_T', 'C_T', 'C_N', 'C_O', 'V_T', 'V_N', 'V_O');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = return_dmat(flag1, flag2, L)
    global N_coll intensity_list indices
    flags = flag1 & flag2;
    K = sum(flags);
    D = zeros(K, N_coll);
    D((indices(flags)-1)*K + (1:K)') = intensity_list(flags);
    D = custom_reshape(D, L);
end

function D = custom_reshape(D, K)
    D = reshape(D, size(D, 1), K, size(D, 2)/K);
    D = sum(D, 3);
end