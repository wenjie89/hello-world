close all; clear; clc;

%% Problem 2

% Problem 2 (a)
x1 = 1; x2 = 0; x3 = 0; x4 = 0; x5 = 0;
x = [x1, x2, x3, x4, x5]';


A = hilb(5);

leftA = sum( abs(A*x) );

rightA = max( sum( abs(A) ) );

% Problem 2 (b)

y1 = 0; y2 = 0; y3 = 0; y4 = 1; y5 = 0;
y = [y1, y2, y3, y4, y5]';


B = inv(A);

leftB = sum( abs(B*y) );

rightB = max( sum( abs(B) ) );

% Problem 2 (c)

z1 = 1; z2 = 1; z3 = 1; z4 = 1; z5 = 1;
z = [z1, z2, z3, z4, z5]';

leftC = max( abs(A*z) );
rightC = max( sum( abs(A),2 ) );

% Problem 2 (d)

w1 = -1; w2 = 1; w3 = -1; w4 = 1; w5 = -1;
w = [w1, w2, w3, w4, w5]';

leftD = max( abs(B*w) )
rightD = max( sum( abs(B),2 ) )


%% 











