close all;
clear;
clc;

a = magic(3)
i = [1, 1, 2,];
j = [1, 1, 2];

a((j-1)*3 + i) = [3, 2, 4];
a
% a([i',j'])
% [i',j']

% rng(2)
% y = round(2 * rand(4, 1) - 1, 1);
% y(3) = 1;
% y
% yr = [-0.9, -0.5, 0, 0.5, 1];
% Ny = numel(yr);
% n_yr = repmat(yr, numel(y), 1);
% [s_n_yr, K] = sort([n_yr, y], 2);
% K == Ny+1
% [row, col] = find(K == Ny + 1);
% A = sortrows([row, col], 1)
% final_y_ids = A(:, 2)
% final_y_ids = final_y_ids - sum(y == yr, 2)