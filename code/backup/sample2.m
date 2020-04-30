
%% Import stuff
load('..\data_x.mat', 'intensity_list', 'P', 'indices')
intensity_list_x = intensity_list;
P_x = P;
indices_x = indices;

load('..\data_y.mat', 'intensity_list', 'P', 'indices')
intensity_list_y = intensity_list;
P_y = P;
indices_y = indices;

load('..\data_xy.mat', 'intensity_list', 'P', 'indices', 'xe_pts', 'ye_pts', 'ze_pts')
intensity_list_xy = intensity_list;
P_xy = P;
indices_xy = indices;

load('..\data_mxy.mat', 'intensity_list', 'P', 'indices', 'xe_pts', 'ye_pts', 'ze_pts')
intensity_list_mxy = intensity_list;
P_mxy = P;
indices_mxy = indices;

%% Plotting
figure;
clims = [0, 3];
% x = x(:, 1:100);
x = rand(1000, 100);
for idx = 1 : size(x, 1)
   idx
   y = x(idx, :) ;
   y1 = y(1:25); y2 = y(26:50); y3 = y(51:75); y4 = y(76:100);
   
   z = zeros(size(find(P_x(:, 3)==0)));
   
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

   t_fl = indices_mxy == 0;
   indices_mxy(t_fl) = 1;
   z3 = intensity_list_mxy .* y3(indices_mxy);
   z3(t_fl) = 0;
   flags = P_mxy(:, 3) == 0;
   z = z + z3(flags);
   
   imagesc(reshape(z, numel(xe_pts), numel(ye_pts)))
   colorbar;   
   
   pause(0.5)
end