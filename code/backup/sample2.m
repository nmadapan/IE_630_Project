%% Plotting
figure;
clims = [0, 3];
for idx = 1 : size(x, 1)
   y = x(idx, :) ;
   y1 = y(1:25); y2 = y(26:50); y3 = y(51:75);
   
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
   
   imagesc(reshape(z, numel(xe_pts), numel(ye_pts)))
   colorbar;   
   
   pause(0.5)
end