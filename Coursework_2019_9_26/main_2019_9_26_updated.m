%% Problem 1
close all; clear; clc;

A = [1.02 1.01; 1.01 1];
b = [2.03 2.01]';
Norm_b = norm(b,Inf);
CondA = cond(A,Inf);

% if b == A*x
%     disp('ok')
% else
%     disp('strange!')
% end


Xhat1 = linspace(2,20,200); 
Xhat2 = linspace(2,20,200);

x = [1 1]';

[xhat1, xhat2] = meshgrid(Xhat1, Xhat2);
xhat = cell(size(xhat1,1),size(xhat1, 2)); % necessary to initialized it!

for i = 1:size(xhat1,1)
    for j = 1:size(xhat1, 2)
      
        xhat{i,j} = [xhat1(i,j) xhat2(i,j)]';
        Norm_rhat(i,j) = norm(b-A*xhat{i,j},Inf);
        Norm_xmxh(i,j) = norm(x-xhat{i,j},Inf);
        goal(i,j) = Norm_xmxh(i,j) - CondA*Norm_rhat(i,j)/Norm_b;

    end
end


figure(1)
clf
imagesc(Xhat1,Xhat2, goal )
colorbar

[goalIdx1, goalIdx2] = find(abs(goal)==max(max(abs(goal))));

goal(goalIdx1, goalIdx2)

xhat{goalIdx1, goalIdx2}



%% Problem 1 (c)

clear; close all; clc;

A = [1.02 1.01; 1.01 1];
b = [2.03 2.01]';
x = [1 1]';
fmax = 1e-4;

f = -fmax:1e-6:fmax;
NormStore1 = zeros(1,length(f));
NormStore2 = zeros(1,length(f));
NormStore3 = zeros(1,length(f));
NormStore4 = zeros(1,length(f));

for i = 1:length(f)
    bhat1 =  b + [-fmax f(i)]';
    y1 = A\bhat1;
    NormStore1(i) = norm(y1-x, Inf);
end

figure(2)
clf
plot(f,NormStore1)
xlim([-1e-4 1e-4])
xlabel('perturbation f_2 when f_1 = -10^{-4}')
ylabel('Infinity norm of y-x')
set(gca, 'FontSize', 15)



for i = 1:length(f)
    bhat2 =  b + [fmax f(i)]';
    y2 = A\bhat2;
    NormStore2(i) = norm(y2-x, Inf);
end

figure(3)
clf
plot(f,NormStore2)
xlim([-1e-4 1e-4])
xlabel('perturbation f_2 when f_1 = 10^{-4}')
ylabel('Infinity norm of y-x')
set(gca, 'FontSize', 15)


for i = 1:length(f)
    bhat3 =  b + [f(i) -fmax]';
    y3 = A\bhat3;
    NormStore3(i) = norm(y3-x, Inf);
end

figure(4)
clf
plot(f,NormStore1)
xlim([-1e-4 1e-4])
xlabel('perturbation f_1 when f_2 = -10^{-4}')
ylabel('Infinity norm of y-x')
set(gca, 'FontSize', 15)



for i = 1:length(f)
    bhat4 =  b + [f(i) fmax]';
    y4 = A\bhat4;
    NormStore4(i) = norm(y4-x, Inf);
end

figure(5)
clf
plot(f,NormStore4)
xlim([-1e-4 1e-4])
xlabel('perturbation f_1 when f_2 = 10^{-4}')
ylabel('Infinity norm of y-x')
set(gca, 'FontSize', 15)



return








y2 = A\bhat2











%% Problem 2
close all; clear; clc;

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





















