function x_dot_k = SINDy_lib(X,Xi_k)
theta_x_k = [ 1 X(1) X(2) ...
    X(1)^2 X(1)*X(2) X(2)^2 ...
    X(1)^3 X(1)^2*X(2) X(1)*X(2)^3 X(2)^3 ...
    sin(X(1)) sin(X(2)) cos(X(1)) cos(X(2)) ...
    sin(2*X(1)) sin(2*X(2)) cos(2*X(1)) sin(2*X(2)) ...
    sin(3*X(1)) sin(3*X(2)) cos(3*X(1)) cos(3*X(2))];
x_dot_k = theta_x_k * Xi_k;

return