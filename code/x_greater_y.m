close all;
clear;
clc;

v1 = [-1, 5, -3];
v2 = [-2, 2, -1];
v3 = [0, 2, 1];

x_greater_than_y(v3, v1, true)
x_greater_than_y(v3, v2, true)

function res = x_greater_than_y(x, y, x_greater)
%%%%%%%%%%%%%%%%%%%%%%%
% Verifies if t * x </> y
%
% Inputs:
% x: 1D array of real numbers
% y: 1D array of real numbers
% x_greater: boolean variable. If true, it evaluates x>y, else x<y.
%   Default value is set to true. 
%
% Return:
%   true or false. 
%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin == 2)
        x_greater = true;
    end
    F = @(m, a, b) (a*m'-b);
    N = numel(x);

    % Find all the points to check
    t = sort([y(:) ./ x(:); 0]); % 0 is added to the list
    nt = 0.5 * ( t(1:end-1) + t(2:end) );
    values = unique([t; nt; max(t)+1]); % 1 + max value is also added. 
    values = values(values > 0);

    res = F(values(:), x(:), y(:));
    if(~x_greater)
        res(res > 0) = inf;
        res(res < 0) = 1;
        res = any(and(sum(res) > 0, sum(res) <= N));
    else
        res(res > 0) = 1;
        res(res < 0) = inf;
        res = any(and(sum(res) > 0, sum(res) <= N));    
    end
end
