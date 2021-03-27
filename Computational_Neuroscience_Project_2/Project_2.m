clear all
close all

%% Setup Inputs
tspan = [0:0.01:200];
y0 = [0,0]; %Initial values of V and w
Iext=0;
C=20;
phi=0.01;

gK=8;
gCa=4.4;
gL=2;

EK=-84;
ECa=120;
EL=-60;

V1=-1.2;
V2=18;
V3=2;
V4=30;
V5=2;
V6=30;

% mss=0.5*(1+tanh((y(1)-V1)/V2));
% wss=0.5*(1+tanh((y(1)-V3)/V4));
% tw=1/cosh((y(1)-V3)/(2*V4));
%dVdt=(Iext-gK*w*(V-EK)-gCa*mss*(V-ECa)-gL*(V-EL))/C;
%dwdt=phi*(wss-w)/tw;
%y=[V,w]

%% MLE

[t,y]=ode45(@(t,y)[(Iext-gK*y(2)*(y(1)-EK)-gCa*(0.5*(1+tanh((y(1)-V1)/V2)))*(y(1)-ECa)-gL*(y(1)-EL))/C; (phi*((0.5*(1+tanh((y(1)-V3)/V4)))-y(2))/(1/cosh((y(1)-V3)/(2*V4))))],tspan,y0);

w_nc=0.5*(1+tanh((y(:,1)-V3)/V4));
V_nc=(Iext-gCa*(0.5*(1+tanh((y(:,1)-V1)/V2))).*(y(:,1)-ECa)-gL*(y(:,1)-EL))./(gK*(y(:,1)-EK));

dVdt = gradient(y(:,1));
dwdt = gradient(y(:,2));

j=0;
Err=1e-4;
for i=1:length(tspan)
    if abs(w_nc(i)-V_nc(i))<=Err
        j=i;
        Eq(j,:)=y(i,:);  
    end
end

d11=gradient((Iext-gK*y(2)*(y(j,1)-EK)-gCa*(0.5*(1+tanh((y(j,1)-V1)/V2)))*(y(j,1)-ECa)-gL*(y(j,1)-EL))/C,y(j,1));
d12=gradient((Iext-gK*y(2)*(y(1)-EK)-gCa*(0.5*(1+tanh((y(1)-V1)/V2)))*(y(1)-ECa)-gL*(y(1)-EL))/C,y(2));
d21=gradient(phi*((0.5*(1+tanh((y(1)-V3)/V4)))-y(2))/(1/cosh((y(1)-V3)/(2*V4))),y(1));
d22=gradient(phi*((0.5*(1+tanh((y(1)-V3)/V4)))-y(2))/(1/cosh((y(1)-V3)/(2*V4))),y(2));
%% Outputs
figure(1)
hold on
plot(y(:,1),w_nc*100);
plot(y(:,1),V_nc*100);
plot(y(:,1),y(:,2)*100); 
quiver(y(:,1),y(:,2)*100,dVdt,dwdt)
legend('w nullcline','V nullcline','trajectory')
hold off