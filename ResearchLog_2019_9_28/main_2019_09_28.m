%% Problem 1
close all; clear; clc; 
delete *.eps

% define parameters

eps0 = 8.85e-12; % vacuum permittivity

Lx = 3e-2; % Lx = 3cm
Ly = 3e-2; % Ly = 3cm

x0 = Lx/2; 
y0 = Ly/2; % center of the cylinder 

R = Lx/6;  % radius of the cylinder
rou0 = 8.85e-12;  % charge density in the cylinder


nx = 100;  % number of points 
ny = nx;


% discretization domain

x = linspace(0,Lx,nx+1);
y = linspace(0,Ly,ny+1);

delx = x(2)-x(1);
dely = y(2)-y(1);

[X, Y] = meshgrid(x,y);

% Compute the electric field and electric potential

r  = sqrt((X-x0).^2+(Y-y0).^2);


% E1 = ro_0 / (2*eps0) * r;       % r < R
% E2 = ro_0*R^2/(2*eps0) ./ r;    % r > R
% 
% Phi1 = -ro_0/(4*eps0) * r.*2;   % r < R
% Phi2 = -ro_0/(4*eps0) * R^2 * (1-2*log(R./r)); % r > R 

E1 = rou0 / (2*eps0) * sqrt((X-x0).^2+(Y-y0).^2);       % r < R
E2 = rou0*R^2/(2*eps0) ./ sqrt((X-x0).^2+(Y-y0).^2);    % r > R

V1 = -rou0/(4*eps0) * ((X-x0).^2+(Y-y0).^2);   % r < R
V2 = -rou0/(4*eps0) * R^2 * (1-2*log(R./sqrt((X-x0).^2+(Y-y0).^2))); % r > R 

[idxX, idxY] = find( sqrt( (X-x0).^2 + (Y-y0).^2 ) > R );


E = E1;
V = V1;

for i = 1:length(idxX) 
    E(idxX(i), idxY(i)) = E2(idxX(i), idxY(i));
    V(idxX(i), idxY(i)) = V2(idxX(i), idxY(i));
end


Ex = E.* (X-x0)./r;   % Ex = Er * (x-x0)/r
Ey = E.* (Y-y0)./r;   % Ey = Er * (y-y0)/r

xidx = find(y == y0);
yidx = find(x == x0);




figure(1)
clf 
subplot(1,2,1)
plot(X(xidx,:),V(xidx,:))
xlabel('x [m]')
ylabel('Electrical Potential [V]')
subplot(1,2,2)
plot(X(xidx,:),E(xidx,:))
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
print -depsc fig1.eps

figure(2)
clf
subplot(1,2,1)
plot(Y(:,yidx),V(:,yidx))
xlabel('y [m]')
ylabel('Electrical Potential [V]')
subplot(1,2,2)
plot(Y(:,yidx),E(:,yidx))
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
print -depsc fig2.eps

figure(3)
clf
hold on 
axis equal
contour(X,Y,V)
quiver(X,Y,Ex,Ey)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig3.eps

figure(4) 
clf
imagesc(x,y,V) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig4.eps

figure(5)
clf
imagesc(x,y,E) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig5.eps


%% Problem 2

% parameter definition
a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
b = a;
c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
d = c;

% v = ones(size(V,1),size(V,2));
% v(1,:) = 0; v(size(v,1),:) = 0; v(:,1) = 0; v(:,size(v,2)) = 0; % boundary condition

v = zeros(size(V,1),size(V,2));

% distribution of the charge density
s = rou0*ones(size(X,1),size(X,2))*(-2/delx^2-2/dely^2)^-1/(-eps0);
for i = 1 : length(idxX)
    s(idxX(i), idxY(i)) = 0;
end

epsilon = 1e-10;  % accuracy
dnorms=sum(sum(s.^2)); 

% figure(6)
% clf
% imagesc(x,y,s)
% axis equal
% colorbar

wVect = 1:0.05:2;
vCell = cell(1,length(wVect));
itVect = zeros(1,length(wVect));
exCell = cell(1,length(wVect));
eyCell = cell(1,length(wVect));
eCell = cell(1,length(wVect));


for i = 1:length(wVect)
    
    w = wVect(i);

    [v,it] = sorEE430S19(a,b,c,d,v,s,w,epsilon,nx,ny,dnorms);
    vCell{i} = v;    % store the potential
    itVect(i) = it;  % store the iteration number


    % calculating the electric field from the potential
    [ex, ey] = gradient(v);
    ex = -ex/delx;
    ey = -ey/dely;
    e = sqrt(ex.^2+ey.^2);
    exCell{i} = ex;  % store ex
    eyCell{i} = ey;  % store ey
    eCell{i} = e;    % store e

end



figure(6)
clf
subplot(1,2,1)
plot(X(xidx,:),v(xidx,:))
hold on
plot(X(xidx,:),V(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Potential [V]')
legend('SOR', 'Analytical')
subplot(1,2,2)
plot(X(xidx,:),e(xidx,:))
hold on
plot(X(xidx,:),E(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
legend('SOR', 'Analytical')
print -depsc fig6.eps

figure(7)
clf
subplot(1,2,1)
plot(Y(:,yidx),v(:,yidx))
hold on
plot(Y(:,yidx),V(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Potential [V]')
legend('SOR', 'Analytical')
subplot(1,2,2)
plot(Y(:,yidx),e(:,yidx))
hold on 
plot(Y(:,yidx),E(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
legend('SOR', 'Analytical')
print -depsc fig7.eps


figure(8)
clf
hold on 
axis equal
contour(X,Y,v)
quiver(X,Y,ex,ey)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig8.eps

figure(9) 
clf
imagesc(x,y,v) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig9.eps

figure(10)
clf
imagesc(x,y,e) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig10.eps



figure(11)
clf
semilogy(wVect,itVect)
hold on
plot(wVect(itVect==min(itVect)),itVect(itVect==min(itVect)),'r*')
xlabel('w')
ylabel('Number of iteration')
print -depsc fig11.eps

wb_num = wVect(itVect==min(itVect));
wb_ana = 2/(1+sqrt(1-cos(pi/nx)^2));




% s = -rou0/eps0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
% rosor = cos(pi/nx); w = 2/(1+sqrt(1-rosor^2));
% epsilon = 1e-10; 
% dnorms=sum(sum(s.^2));
% 
% [v,it]=sorEE430S19(a,b,c,d,v,s,w,epsilon,nx,ny,dnorms);
% [ex, ey] = gradient(v);
% ex = -ex/delx;
% ey = -ey/dely;
% e = sqrt(ex.^2+ey.^2);


% clear w a b c d v s epsilon dnorms
% 
% w = 1.1:0.1:2.9;
% itVect = zeros(1,length(w));
% for i = 1:length(w)
%     a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
%     b = a;
%     c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
%     d = c;
%     v = zeros(size(V,1),size(V,2));
%     s = -rou0/eps0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
%     epsilon = 1e-10; 
%     dnorms=sum(sum(s.^2));
%     [v,it]=sorEE430S19(a,b,c,d,v,s,w(i),epsilon,nx,ny,dnorms);
%     itVect(i) = it;
% end





%% Problem 3


v = zeros(size(V,1),size(V,2));

% match boundary conditions
v(1,:) = V(1,:) ;
v(size(v,1),:) = V(size(v,1),:); 
v(:,1) = V(:,1);
v(:,size(v,2)) = V(:,size(v,2)); 

w = wb_num;   % set it as the best value calculated numerically 

[v,~] = sorEE430S19(a,b,c,d,v,s,w,epsilon,nx,ny,dnorms);


% calculating the electric field from the potential
[ex, ey] = gradient(v);
ex = -ex/delx;
ey = -ey/dely;
e = sqrt(ex.^2+ey.^2);


figure(12)
clf
subplot(1,2,1)
plot(X(xidx,:),v(xidx,:))
hold on
plot(X(xidx,:),V(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Potential [V]')
legend('SOR', 'Analytical')
subplot(1,2,2)
plot(X(xidx,:),e(xidx,:))
hold on
plot(X(xidx,:),E(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
legend('SOR', 'Analytical')
print -depsc fig12.eps

figure(13)
clf
subplot(1,2,1)
plot(Y(:,yidx),v(:,yidx))
hold on
plot(Y(:,yidx),V(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Potential [V]')
legend('SOR', 'Analytical')
subplot(1,2,2)
plot(Y(:,yidx),e(:,yidx))
hold on 
plot(Y(:,yidx),E(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
legend('SOR', 'Analytical')
print -depsc fig13.eps


figure(14)
clf
hold on 
axis equal
contour(X,Y,v)
quiver(X,Y,ex,ey)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig14.eps

figure(15) 
clf
imagesc(x,y,v) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig15.eps

figure(16)
clf
imagesc(x,y,e) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig16.eps





