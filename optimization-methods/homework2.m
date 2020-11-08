clc; clear;

syms x1 x2;

y = 100*(x2-x1^2)^2 + (1-x1)^2;

diff(y, x1, 1)
diff(y, x2, 1)
%diff(y, 2)