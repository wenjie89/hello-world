%% Problem 1
close all; clear; clc;

e0 = 8.85e-12; % [E/m]
Lx = 0.03; %[m]
Ly = Lx; 
R = Lx/6;
x0 = Lx/2; y0 = Ly/2; % center of the cylinder
rou0 = 8.85e-12; % charge density [C/m^3]

nx = 20; % number of intevals
ny = nx;
delx = Lx/nx;
dely = Ly/ny;

x = 0:delx:Lx;
y = 0:dely:Ly;

[X, Y] = meshgrid(x,y);
r = ((X-x0).^2+(Y-y0).^2).^0.5;

E = zeros(size(X,1),size(X,2));
V = zeros(size(X,1),size(X,2));

for i = 1:size(X,1)
    for j = 1;size(X,2)
        if ((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5 < R || ...
                ((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5 == R || ((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5 > R 
            E(i,j) = rou0/(2*e0)*((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5;
        else
            E(i,j) = rou0/(2*e0)*R^2./((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5;    
        end
    end
end


for i = 1:size(X,1)
    for j = 1;size(X,2)
        if ((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5 < R || ...
                ((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5 == R || ((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5 > R 
            V(i,j) = -rou0/(4*e0)*((X(i,j)-x0)^2+(Y(i,j)-y0)^2);
        else
            V(i,j) = -rou0/(4*e0)*R^2*(1-2*log(R./((X(i,j)-x0)^2+(Y(i,j)-y0)^2)^0.5));    
        end
    end
end
        





[Ex, Ey] = gradient(V);
Ex = -Ex/delx;
Ey = -Ey/dely;
xidx = find(Y(:,1)==Ly/2);
yidx = find(X(1,:)==Lx/2);

figure(1)
subplot(1,2,1)
plot(X(xidx,:),V(xidx,:))
xlabel('x [m]')
ylabel('Electrical Potential [V]')
subplot(1,2,2)
plot(X(xidx,:),E(xidx,:))
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
%print -depsc fig1.eps

figure(2)
subplot(1,2,1)
plot(Y(:,yidx),V(:,yidx))
xlabel('y [m]')
ylabel('Electrical Potential [V]')
subplot(1,2,2)
plot(Y(:,yidx),E(:,yidx))
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
%print -depsc fig2.eps

figure(3)
hold on 
axis equal
contour(X,Y,V)
quiver(X,Y,Ex,Ey)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
%print -depsc fig3.eps

figure(4) 
imagesc(x,y,V) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
%print -depsc fig4.eps

figure(5)
imagesc(x,y,E) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
%print -depsc fig5.eps

return
%% Problem 2

a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
b = a;
c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
d = c;
v = ones(size(V,1),size(V,2));
v(1,:) = 0; v(size(v,1),:) = 0; v(:,1) = 0; v(:,size(v,2)) = 0; % boundary condition
s = -rou0/e0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
rosor = cos(pi/nx); w = 2/(1+sqrt(1-rosor^2));
epsilon = 1e-10; 
dnorms=sum(sum(s.^2));

[v,it]=sorEE430S19(a,b,c,d,v,s,w,epsilon,nx,ny,dnorms);
[ex, ey] = gradient(v);
ex = -ex/delx;
ey = -ey/dely;
e = sqrt(ex.^2+ey.^2);

figure(6)
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

figure(88)
hold on 
axis equal
contour(X,Y,v)
quiver(X,Y,ex,ey)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig88.eps

figure(99) 
imagesc(x,y,v) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig99.eps

figure(100)
imagesc(x,y,e) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig100.eps

clear w a b c d v s epsilon dnorms

w = 1.1:0.1:2.9;
IT = zeros(1,length(w));
for i = 1:length(w)
    a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
    b = a;
    c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
    d = c;
    v = zeros(size(V,1),size(V,2));
    s = -rou0/e0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
    epsilon = 1e-10; 
    dnorms=sum(sum(s.^2));
    [v,it]=sorEE430S19(a,b,c,d,v,s,w(i),epsilon,nx,ny,dnorms);
    IT(i) = it;
end

figure(9)
plot(w,IT)
xlabel('w')
ylabel('Number of iteration')
w_opt = w(IT==min(IT)); % optimal w
print -depsc fig9.eps


%% problem 3

a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
b = a;
c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
d = c;
v = zeros(size(V,1),size(V,2));

% match boundary conditions
v(1,:) = V(1,:) ;
v(size(v,1),:) = V(size(v,1),:); 
v(:,1) = V(:,1);
v(:,size(v,2)) = V(:,size(v,2)); 

s = -rou0/e0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
epsilon = 1e-30; 
dnorms=sum(sum(s.^2));
[v_opt,it]=sorEE430S19(a,b,c,d,v,s,w_opt,epsilon,nx,ny,dnorms);
[ex_opt, ey_opt] = gradient(v_opt);
ex_opt = -ex_opt/delx;
ey_opt = -ey_opt/dely;
e_opt = sqrt(ex_opt.^2+ey_opt.^2);

figure(10)
subplot(1,2,1)
plot(X(xidx,:),v_opt(xidx,:))
hold on 
plot(X(xidx,:),V(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Potential [V]')
legend('SOR (optimized)','Analytical')
subplot(1,2,2)
plot(X(xidx,:),e_opt(xidx,:))
hold on
plot(X(xidx,:),E(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
legend('SOR (optimized)','Analytical')
print -depsc fig10.eps

figure(11)
subplot(1,2,1)
plot(Y(:,yidx),v_opt(:,yidx))
hold on
plot(Y(:,yidx),V(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Potential [V]')
legend('SOR (optimized)','Analytical')
subplot(1,2,2)
plot(Y(:,yidx),e_opt(:,yidx))
hold on
plot(Y(:,yidx),E(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
legend('SOR (optimized)','Analytical')
print -depsc fig11.eps

figure(120)
hold on 
axis equal
contour(X,Y,v_opt)
quiver(X,Y,ex_opt,ey_opt)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig120.eps

figure(130) 
imagesc(x,y,V) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig130.eps

figure(140)
imagesc(x,y,E) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig140.eps

%% Problem 4
close all; clear; clc;
delete *.eps

e0 = 8.85e-12; % [E/m]
Lx = 0.03; %[m]
Ly = Lx; 
R = Lx/6;
x0 = Lx/4; y0 = Ly/4;
rou0 = 8.85e-12; %[C/m^3]

nx = 20; % number of intevals
ny = nx;
delx = Lx/nx;
dely = Ly/ny;

x = 0:delx:Lx;
y = 0:dely:Ly;

[X, Y] = meshgrid(x,y);

r = ((X-x0).^2+(Y-y0).^2).^0.5;

V = -rou0/(4*e0)*r.^2;
E = rou0/(2*e0)*r;

[Ex, Ey] = gradient(V);
Ex = -Ex/delx;
Ey = -Ey/dely;

xidx = find(Y(:,1)==Ly/4);
yidx = find(X(1,:)==Lx/4);

figure(12)
subplot(1,2,1)
plot(X(xidx,:),V(xidx,:))
xlabel('x [m]')
ylabel('Electrical Potential [V]')
subplot(1,2,2)
plot(X(xidx,:),E(xidx,:))
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
print -depsc fig12.eps

figure(13)
subplot(1,2,1)
plot(Y(:,yidx),V(:,yidx))
xlabel('y [m]')
ylabel('Electrical Potential [V]')
subplot(1,2,2)
plot(Y(:,yidx),E(:,yidx))
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
print -depsc fig13.eps

figure(14)
hold on 
axis equal
contour(X,Y,V)
quiver(X,Y,Ex,Ey)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig14.eps

figure(16) 
imagesc(x,y,V) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig16.eps

figure(17)
imagesc(x,y,E) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig17.eps

a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
b = a;
c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
d = c;
v = ones(size(V,1),size(V,2));
v(1,:) = 0; v(size(v,1),:) = 0; v(:,1) = 0; v(:,size(v,2)) = 0; % boundary condition
s = -rou0/e0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
rosor = cos(pi/nx); w = 2/(1+sqrt(1-rosor^2));
epsilon = 1e-10; 
dnorms=sum(sum(s.^2));

[v,it]=sorEE430S19(a,b,c,d,v,s,w,epsilon,nx,ny,dnorms);
[ex, ey] = gradient(v);
ex = -ex/delx;
ey = -ey/dely;
e = sqrt(ex.^2+ey.^2);

a = 1/dely^2*(-2/delx^2-2/dely^2)^-1;
b = a;
c = 1/delx^2*(-2/delx^2-2/dely^2)^-1;
d = c;
v = zeros(size(V,1),size(V,2));

% match boundary conditions
v(1,:) = V(1,:) ;
v(size(v,1),:) = V(size(v,1),:); 
v(:,1) = V(:,1);
v(:,size(v,2)) = V(:,size(v,2)); 

s = -rou0/e0*(-2/delx^2-2/dely^2)^(-1)*ones(size(v,1),size(v,2));
epsilon = 1e-30; 
dnorms=sum(sum(s.^2));

[v_opt,it]=sorEE430S19(a,b,c,d,v,s,w,epsilon,nx,ny,dnorms);
[ex_opt, ey_opt] = gradient(v_opt);
ex_opt = -ex_opt/delx;
ey_opt = -ey_opt/dely;
e_opt = sqrt(ex_opt.^2+ey_opt.^2);

figure(18)
subplot(1,2,1)
plot(X(xidx,:),v_opt(xidx,:))
hold on 
plot(X(xidx,:),V(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Potential [V]')
legend('SOR (optimized)','Analytical')
subplot(1,2,2)
plot(X(xidx,:),e_opt(xidx,:))
hold on
plot(X(xidx,:),E(xidx,:),'r--')
xlabel('x [m]')
ylabel('Electrical Field [V/m]')
legend('SOR (optimized)','Analytical')
print -depsc fig18.eps

figure(19)
subplot(1,2,1)
plot(Y(:,yidx),v_opt(:,yidx))
hold on
plot(Y(:,yidx),V(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Potential [V]')
legend('SOR (optimized)','Analytical')
subplot(1,2,2)
plot(Y(:,yidx),e_opt(:,yidx))
hold on
plot(Y(:,yidx),E(:,yidx),'r--')
xlabel('y [m]')
ylabel('Electrical Field [V/m]')
legend('SOR (optimized)','Analytical')
print -depsc fig19.eps

figure(20)
hold on 
axis equal
contour(X,Y,v_opt)
quiver(X,Y,ex_opt,ey_opt)
legend('Potential', 'Field')
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig20.eps

figure(21) 
imagesc(x,y,v_opt) % magnitude of potential
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig21.eps

figure(22)
imagesc(x,y,e_opt) % magnitude of electric field
colorbar
xlabel('x [m]')
ylabel('y [m]')
print -depsc fig22.eps

