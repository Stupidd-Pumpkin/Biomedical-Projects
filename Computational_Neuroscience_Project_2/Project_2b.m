clear all

%% Setup Inputs
global para iext
%para=[gca(1), gk(2), gl(3), vca(4), vk(5), vl(6), phi(7), v1(8), v2(9), v3(10), v4(11), v5(12), v6(13), C(14)]'
% para=[4.4, 8.0, 2, 120, -84, -60, 0.04, -1.2, 18, 2, 30, 2, 30, 20];iext = [95,0,2000];
para=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];iext = [0,0,2000];
% para=[4.0, 8.0, 2, 120, -84, -60, 0.0667, -1.2, 18, 12, 17.4, 12, 17.4, 20];iext = [0,0,2000];
tspan = [0, 2000]; %Tstart, Tstop


%% MLE
x1=-80:1:80;
I=1;Err=10;
for i=x1
    
    w_nc(I)=0.5*(1+tanh((i-para(10))/para(11)));
    V_nc(I)=(iext(1)-para(1)*(0.5*(1+tanh((i-para(8))/para(9))))*(i-para(4))-para(3)*(i-para(6)))/(para(2)*(i-para(5)));
    j=abs(w_nc(I)-V_nc(I));
    if j<Err
        Err=j;
        Eq=[i,w_nc(I)];
    end
    I=I+1;
end

x2=-80:10:80;
y2=0:0.05:0.8;
mx2=meshgrid(x2)';
my2=meshgrid(y2);
I=0;J=0;
for i=x2
    tw = 1/cosh((i-para(12))/(2*para(13)));
    ms = 0.5*(1 + tanh((i-para(8))/para(9)));
    ws = 0.5*(1 + tanh((i-para(10))/para(11)));
    I=I+1;
    for j=y2
        J=J+1;
        dydt1(I,J) = (iext(1)-para(1)*ms*(i-para(4)) - para(2)*j*(i-para(5)) - para(3)*(i-para(6)))/para(14);
        dydt2(I,J) = para(7)*(ws-j)/tw;
    end
    J=0;
end

Eq
y0=Eq;   %initial conditions [V0,w0]
options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0, options);  %ODE solver

%% Outputs
figure(1); clf
[ax, h1, h2] = plotyy(t, y(:,1), t, y(:,2));
axes(ax(1)); axis([0 tspan(2) -80 80]); ylabel('V, mV')
axes(ax(2)); axis([0 tspan(2) 0 0.5]); ylabel('w')
xlabel('Time, ms');

figure(2); clf
hold on
plot(x1, w_nc*100,'r--');
plot(x1, V_nc*100,'g--');
quiver(mx2, my2*100, dydt1, dydt2*100,'b');
plot(y(:,1), y(:,2)*100);
quiver(y(:,1), y(:,2)*100, gradient(y(:,1)), gradient(y(:,2))*100)
legend('w nullcline','V nullcline','trajectory')
xlabel('V, mV');ylabel('W')
axis([-80,80,-20,100])
hold off

%% Functions
function dydt = MLE(t, y)

global para iext

tw = 1./cosh((y(1)-para(12))/(2*para(13)));
ms = 0.5*(1 + tanh((y(1)-para(8))/para(9)));
ws = 0.5*(1 + tanh((y(1)-para(10))/para(11)));

dydt = zeros(2,1);
if t>=iext(2) & t<iext(3)
    dydt(1) = (iext(1) - para(1)*ms.*(y(1)-para(4)) - para(2)*y(2).*(y(1)-para(5)) - para(3)*(y(1)-para(6)))/para(14);
else
    dydt(1) = (-para(1)*ms.*(y(1)-para(4)) -  para(2)*y(2).*(y(1)-para(5)) - para(3)*(y(1)-para(6)))/para(14);
end
dydt(2) = para(7)*(ws-y(2))./tw;
end

function jac = JAC(t,y)

global para iext

tw = 1./cosh((y(1)-para(12))/(2*para(13)));
ms = 0.5*(1 + tanh((y(1)-para(8))/para(9)));
ws = 0.5*(1 + tanh((y(1)-para(10))/para(11)));

jac = zeros(2,2);
jac(1,1) = (-para(1)*0.5/(para(9)*cosh((y(1)-para(8))/para(9))^2)*(y(1)-para(4)) - para(1)*ms - para(2)*y(2) - para(3))/para(14);
jac(1,2) = -para(2)*(y(1)-para(5))/para(14);
jac(2,1) = para(7)*(0.5/ (para(11)*cosh((y(1)-para(10))/para(11))^2*tw) + (ws-y(2))*sinh((y(1)-para(12))/(2*para(13)))/(2*para(13)));
jac(2,2) = -para(7)/tw;
end
