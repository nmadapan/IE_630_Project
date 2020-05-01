close all;
clear;
clc;

% load('x_80.mat')
% x = [x(:, 1:75), zeros(size(x, 1), 25)];

load('x_106.mat')
x = x(:, 1:100);

%% Import stuff
load('..\data_x.mat', 'intensity_list', 'P', 'indices', 'in_tumor_flags', 'in_oar_flags', 'in_normal_flags')
intensity_list_x = intensity_list;
P_x = P;
indices_x = indices;
in_tumor_flags_x = in_tumor_flags;
in_normal_flags_x = in_normal_flags;
in_oar_flags_x = in_oar_flags;

load('..\data_y.mat', 'intensity_list', 'P', 'indices', 'in_tumor_flags', 'in_oar_flags', 'in_normal_flags')
intensity_list_y = intensity_list;
P_y = P;
indices_y = indices;
in_tumor_flags_y = in_tumor_flags;
in_normal_flags_y = in_normal_flags;
in_oar_flags_y = in_oar_flags;

load('..\data_xy.mat', 'intensity_list', 'P', 'indices', 'xe_pts', 'ye_pts', 'ze_pts', 'in_tumor_flags', 'in_oar_flags', 'in_normal_flags')
intensity_list_xy = intensity_list;
P_xy = P;
indices_xy = indices;
in_tumor_flags_xy = in_tumor_flags;
in_normal_flags_xy = in_normal_flags;
in_oar_flags_xy = in_oar_flags;

load('..\data_mxy.mat', 'intensity_list', 'P', 'indices', 'xe_pts', 'ye_pts', 'ze_pts', 'in_tumor_flags', 'in_oar_flags', 'in_normal_flags')
intensity_list_mxy = intensity_list;
P_mxy = P;
indices_mxy = indices;
in_tumor_flags_mxy = in_tumor_flags;
in_normal_flags_mxy = in_normal_flags;
in_oar_flags_mxy = in_oar_flags;

%% Plotting
figure;
clims = [0, 3];
x = x(:, 1:100);
values = [66, 60, 61, 45];
for idx = values
   idx
   y = x(idx, :) ;
   y1 = y(1:25); y2 = y(26:50); y3 = y(51:75); y4 = y(76:100);
   
   z = zeros(size(find(P_x(:, 3)==0)));
   z = reshape(z, numel(xe_pts), numel(ye_pts));
   
   % Other variables
   tumor_intensity = 0.0;
   normal_intensity = 0.0;
   oar_intensity = 0.0;
   
   % x-direction
   t_fl = indices_x == 0;
   indices_x(t_fl) = 1;
   z1 = intensity_list_x .* y1(indices_x);
   z1(t_fl) = 0;
   flags = P_x(:, 3) == 0;
   image_x = reshape(z1(flags), numel(xe_pts), numel(ye_pts));
   z = z + image_x + fliplr(image_x);
   tumor_intensity = tumor_intensity + sum(z1(in_tumor_flags_x));
   normal_intensity = normal_intensity + sum(z1(in_normal_flags_x));
   oar_intensity = oar_intensity + sum(z1(in_oar_flags_x));
   
   % y direction
   t_fl = indices_y == 0;
   indices_y(t_fl) = 1;
   z2 = intensity_list_y .* y2(indices_y);
   z2(t_fl) = 0;
   flags = P_y(:, 3) == 0;
   image_y = reshape(z2(flags), numel(xe_pts), numel(ye_pts));
   z = z + image_y + flipud(image_y);
   tumor_intensity = tumor_intensity + sum(z2(in_tumor_flags_y));
   normal_intensity = normal_intensity + sum(z2(in_normal_flags_y));
   oar_intensity = oar_intensity + sum(z2(in_oar_flags_y));

   % xy direction
   t_fl = indices_xy == 0;
   indices_xy(t_fl) = 1;
   z3 = intensity_list_xy .* y3(indices_xy);
   z3(t_fl) = 0;
   flags = P_xy(:, 3) == 0;
   image_xy = reshape(z3(flags), numel(xe_pts), numel(ye_pts));
   z = z + image_xy + flipud(image_xy);
   tumor_intensity = tumor_intensity + sum(z3(in_tumor_flags_xy));
   normal_intensity = normal_intensity + sum(z3(in_normal_flags_xy));
   oar_intensity = oar_intensity + sum(z3(in_oar_flags_xy));

   % mxy direction
   t_fl = indices_mxy == 0;
   indices_mxy(t_fl) = 1;
   z4 = intensity_list_mxy .* y4(indices_mxy);
   z4(t_fl) = 0;
   flags = P_mxy(:, 3) == 0;
   image_mxy = reshape(z4(flags), numel(xe_pts), numel(ye_pts));
   z = z + image_mxy + flipud(image_mxy);
   tumor_intensity = tumor_intensity + sum(z4(in_tumor_flags_mxy));
   normal_intensity = normal_intensity + sum(z4(in_normal_flags_mxy));
   oar_intensity = oar_intensity + sum(z4(in_oar_flags_mxy));
   
   imagesc(z);
   colorbar;   
   fprintf('Tumor: %.02f\n', tumor_intensity)
   fprintf('Normal: %.02f\n', normal_intensity)
   fprintf('OAR: %.02f\n\n', oar_intensity)
   
   wfname = fullfile('results', ['img_', num2str(idx), '.png']);
   saveas(gcf, wfname)
   
%    pause(0.5)
end