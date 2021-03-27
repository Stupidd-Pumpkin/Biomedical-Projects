clear all
close all

[x,y] = meshgrid(-2:0.1:2);
f = x.*exp(-x.^2 - y.^2);
[dx,dy] = gradient(f,0.2,0.2);

figure(1)
contour(x,y,f)
hold on
quiver(x,y,dx,dy)
hold off

figure(2)
mesh(x,y,f)
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
