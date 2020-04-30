function [X, d]=point_to_plane(plane, point)
    a = plane(1);
    b = plane(2);
    c = plane(3);
    d = plane(4);
    
    x = point(:, 1);
    y = point(:, 2);
    z = point(:, 3);

    % given an plane equation ax + by + cz + d = 0, project points xyz onto the plane
    % return the coordinates of the new projected points
    % written by Neo Jing Ci, 11/7/18
    x  = x(:)'; y = y(:)'; z = z(:)';
    
    d = -1*d;
    A=[1 0 0 -a; 0 1 0 -b; 0 0 1 -c; a b c 0];
    
    B=[x; y; z; d * ones(1, numel(x))];
    X=A\B;
    X = X(1:3, :);

    d = sqrt(sum((B(1:3, :) - X).^2, 1));
    
    X = X'; d = d';
end