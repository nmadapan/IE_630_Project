clc;
close all;
clear;

%% Import data
load('mlop_data.mat')

load('data_x.mat', 'intensity_list', 'P', 'indices')
intensity_list_x = intensity_list;
P_x = P;
indices_x = indices;

load('data_y.mat', 'intensity_list', 'P', 'indices')
intensity_list_y = intensity_list;
P_y = P;
indices_y = indices;

load('data_xy.mat', 'intensity_list', 'P', 'indices', 'xe_pts', 'ye_pts', 'ze_pts')
intensity_list_xy = intensity_list;
P_xy = P;
indices_xy = indices;

% Plot the heatmap
z = P_xy(:, 3);
flags = z == 0;
X = reshape(intensity_list_xy(flags), numel(xe_pts), numel(ye_pts));
figure;
imagesc(X)
colorbar;

% intensity_list_mxy = reshape(intensity_list, [numel(xe_pts), numel(ye_pts), numel(ze_pts)]);
% intensity_list_mxy = intensity_list_mxy(end:-1:1, :, :);
% intensity_list_mxy = intensity_list_mxy(:);
% P_mxy = P;
% indices_mxy = reshape(indices, [numel(xe_pts), numel(ye_pts), numel(ze_pts)]);
% indices_mxy = indices_mxy(end:-1:1, :, :);
% indices_mxy = indices_mxy(:);
% 
% % Plot the heatmap
% z = P_mxy(:, 3);
% flags = z == 0;
% X = reshape(intensity_list_mxy(flags), numel(xe_pts), numel(ye_pts));
% figure;
% imagesc(X)
% colorbar;

%%
bo1 = 0.1*20.1*0.001*(4/3)* pi * (1*0.5*2)
bo2 = 0.1*20.1*0.001*(4/3) * pi * (0.5*0.5*0.5)
bo3 = 0.1*20.1*0.001*(4/3) * pi * (0.5*(1/3)*0.25)
bn1 = 0.2*17.2*0.001*(10*2*2)
bn3 = 0.2*17.2*0.001*(5*2*2)
bn2 = 0.2*17.2*0.001*(6*2*2)
bt = 1.05*27*0.001*(4/3) * pi * (1)*(1/2)

lb = zeros(size(C_N)); % Lower bound

b = 1e-7 * [bo1;bo2;bo3;bn1;bn2;bn3;-1*bt] * 100;
C= [C_O; C_N; -1 * Dt_T];
A = [Dt_O1; Dt_O2; Dt_O3; Dt_N1; Dt_N2; Dt_N3; -1 * Dt_T];

fitnessfcn = @(x)[x*C'];
options = optimoptions('gamultiobj','UseVectorized',true);
x = gamultiobj(fitnessfcn,numel(C_N),A,b,[],[],lb,[],options);

%% Plotting
% figure;
% for idx = 1 : size(x, 1)
%    y = x(idx, :) ;
%    y1 = y(1:25); y2 = y(26:50); y3 = y(51:75);
%    subplot(1, 3, 1)
%    imagesc(reshape(y1, 5, 5));
%    colorbar;
%    subplot(1, 3, 2)
%    imagesc(reshape(y2, 5, 5));
%    colorbar;
%    subplot(1, 3, 3)
%    imagesc(reshape(y3, 5, 5));
%    colorbar;
%     pause(0.5)
% end

%% Plotting
figure;
clims = [0, 1e-3];
for idx = 1 : size(x, 1)
   y = x(idx, :) ;
   y1 = y(1:25); y2 = y(26:50); y3 = y(51:75);
   
   t_fl = indices_x == 0;
   indices_x(t_fl) = 1;
   z1 = intensity_list_x .* y1(indices_x);
   z1(t_fl) = 0;
   flags = P_x(:, 3) == 0;
   z = z1(flags);

   t_fl = indices_y == 0;
   indices_y(t_fl) = 1;
   z2 = intensity_list_y .* y2(indices_y);
   z2(t_fl) = 0;
   flags = P_y(:, 3) == 0;
   z = z + z2(flags);

   t_fl = indices_xy == 0;
   indices_xy(t_fl) = 1;
   z3 = intensity_list_xy .* y3(indices_xy);
   z3(t_fl) = 0;
   flags = P_xy(:, 3) == 0;
   z = z + z3(flags);
   
   imagesc(reshape(z, numel(xe_pts), numel(ye_pts)), clims)
   colorbar;   
   
   pause(0.5)
end