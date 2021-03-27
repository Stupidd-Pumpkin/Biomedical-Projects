clear all
close all
warning off

%% Setup Inputs
global para iext
%para=[gca(1), gk(2), gl(3), vca(4), vk(5), vl(6), phi(7), v1(8), v2(9), v3(10), v4(11), v5(12), v6(13), C(14)]' iext = [i_val, i_Tstart, i_Tstop]
tspan = [0, 2000]; %Tstart, Tstop

% Question 2,3

fprintf('\nQ2\n');
para=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];iext = [0,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=Eq(1:2);   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);
Ref_Eq=Eq(1:2);

%% Question 5

fprintf('\nQ5\n');
para=[4.4, 8.0, 2, 120, -84, -60, 0.04, -1.2, 18, 2, 30, 2, 30, 20];iext = [0,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=Eq(1:2);   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

%% Question 6

fprintf('\nQ6\n');
para=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];iext = [0,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=[-14.96712078,Eq(2)];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

para=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];iext = [0,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=[-14.96712076,Eq(2)];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

for i=1:50
    y0=[-14.9671208+i*2e-9,Eq(2)];
    [t,y] = ode15s(@MLE, tspan, y0);
    y_max(i)=max(y(:,1));
end
figure;
plot(-14.9671208+[2:2:100]*1e-9,y_max);
xlabel('Initial value of V, mV');ylabel('Maximum value of the action potential, mV')

%% Question 7,8

fprintf('\nQ7\n');
para=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];iext = [86,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=Ref_Eq;   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

y0=Eq(1:2);   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

y0=[-27.9, 0.17];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

y0=[-27.9, 0.17];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE_b, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

for i=1:50
    y0=[Eq(1)+i/5, Eq(2)];   %initial conditions [V0,w0]
    % options = odeset('Jacobian', @JAC,'Vectorized', 'on');
    [t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
    y_max(i)=max(y(:,1));
end

figure;
plot(-27.9+[0.2:0.2:10],y_max);
xlabel('Initial value of V, mV');ylabel('Maximum value of the action potential, mV')

%% Question 9

fprintf('\nQ9\n');
para=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];

iext = [80,tspan(1),tspan(2)];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=[Eq(1)+0.1,Eq(2)];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

iext = [86,tspan(1),tspan(2)];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=[Eq(1)+0.1,Eq(2)];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

iext = [90,tspan(1),tspan(2)];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
y0=[Eq(1)+0.1,Eq(2)];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);

for i=1:20
    iext = [80+i,tspan(1),tspan(2)];
    [x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots
    y0=[Eq(1)+0.1,Eq(2)];   %initial conditions [V0,w0]
    % options = odeset('Jacobian', @JAC,'Vectorized', 'on');
    [t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
    for j=1:length(t)
        if y(j,1)<=0
            y(j,1)=0;
        end
    end
    y_rate(i)=length(findpeaks(y(:,1)));
end

figure;
plot(80+[1:20],y_rate);
xlabel('Value of Iext, uA/cm2');ylabel('Rate of firing action potentials')

%% Question 10

fprintf('\nQ10\n');
para=[4.0, 8.0, 2, 120, -84, -60, 0.0667, -1.2, 18, 12, 17.4, 12, 17.4, 20];
iext = [30,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots

Eq2(1)=fzero(@Eq_Find,Eq(1)-20);
Eq2(2)= 0.5*(1+tanh((Eq2(1)-para(10))/para(11)));
Eq3(1)=fzero(@Eq_Find,Eq2(1)-20);
Eq3(2)= 0.5*(1+tanh((Eq3(1)-para(10))/para(11)));
jac=JAC(Eq2);    %fidning the jacobian matrix at equilibrium
lambda2=eig(jac); %finding the eigenvalues at equilibrium
jac=JAC(Eq3);    %fidning the jacobian matrix at equilibrium
lambda3=eig(jac); %finding the eigenvalues at equilibrium

y0=[-40.3309,0.00243578];;   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver
Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y);
Output(x1,w_nc,V_nc,Eq3,lambda3,mx2,my2,Grad_V,Grad_w,tspan,t,y);
close;close;close;close;
Output(x1,w_nc,V_nc,Eq2,lambda2,mx2,my2,Grad_V,Grad_w,tspan,t,y);
hold on
plot(Eq(1),Eq(2)*100,'bo');
plot(Eq3(1),Eq3(2)*100,'bo');
text(Eq(1),Eq(2)*100-2,['(' num2str(Eq(1)) ',' num2str(Eq(2)*100) ')']);
text(Eq3(1),Eq3(2)*100-2,['(' num2str(Eq3(1)) ',' num2str(Eq3(2)*100) ')']);
hold off

%% Question 11

fprintf('\nQ11\n');
para=[4.0, 8.0, 2, 120, -84, -60, 0.0667, -1.2, 18, 12, 17.4, 12, 17.4, 20];

for i=1:20
    iext = [30+i,0,2000];
[x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup(); %Setup function returns the nullclines, equilibrium points and eigenvalues and other parameters for plots

Eq2(1)=fzero(@Eq_Find,Eq(1)-20);
Eq2(2)= 0.5*(1+tanh((Eq2(1)-para(10))/para(11)));
Eq3(1)=fzero(@Eq_Find,Eq2(1)-20);
Eq3(2)= 0.5*(1+tanh((Eq3(1)-para(10))/para(11)));
jac=JAC(Eq2);    %fidning the jacobian matrix at equilibrium
lambda2=eig(jac); %finding the eigenvalues at equilibrium
jac=JAC(Eq3);    %fidning the jacobian matrix at equilibrium
lambda3=eig(jac); %finding the eigenvalues at equilibrium

y0=[-40.3309,0.00243578];   %initial conditions [V0,w0]
% options = odeset('Jacobian', @JAC,'Vectorized', 'on');
[t,y] = ode15s(@MLE, tspan, y0);%, options);  %ODE solver

for j=1:length(t)
        if y(j,1)<=0
            y(j,1)=0;
        end
    end
    y_rate(i)=length(findpeaks(y(:,1)));
end

figure;
plot(30+[1:20],y_rate);
xlabel('Value of Iext, uA/cm2');ylabel('Rate of firing action potentials')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long

%para=[Gkb(1),Ek(2),GNab(3),ENa(4),GLb(5),EL(6),phi(7),C(8),fni(9)];
%y=[V,n,m,h];
%iext = [i_val, i_Tstart, i_Tstop]
global para iext
global Vss
tspan = [0, 500]; %Tstart, Tstop
Ref_Eq=[-60,0.317676914060697,0.052932485257250,0.596120753508460];
Vss=Ref_Eq(1);
iext = [0, tspan(1),tspan(2)];
T=6.3;

%% Q14

fprintf('\nQ14\n');
para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
iext = [0, 100,125];
V=fzero(@Eq_V,-60);
Eq(:,1)=Eq_Find_HH(V)
y0=Ref_Eq+[0,0,0,0];
[t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver

%% Q15

fprintf('\nQ15\n');
para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];

for k=[12.575757,12.575758]
    iext = [k, tspan(1),tspan(2)];
    V=fzero(@Eq_V,-60);
    Eq(:,1)=Eq_Find_HH(V)
    y0=Eq(:,1)+[0.1,0,0,0]';
    [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
    
    figure
    plot(t,y(:,1))
end

%% Q16

fprintf('\nQ16\n');

for k=[0,0.1,0.17,0.2]
    para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,k];
    V=fzero(@Eq_V,-60);
    Eq(:,1)=Eq_Find_HH(V);
    
    j=2;i=1;
    while i>0
        V=fzero(@Eq_V,Eq(1,j-1)+60);
        Eq(:,j)=Eq_Find_HH(V);
        if Eq(1,j)-Eq(1,j-1)<1e-5
            i=0;
        end
        j=j+1;
    end
    
    Eq(:,1:end-1)
    y0=Ref_Eq+[10,0,0,0];
    [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
    
    figure(1)
    set(gca,'fontsize',24);
    subplot(3,1,1);
    plot(t,y(:,2))
    hold on
    subplot(3,1,2);
    plot(t,y(:,3))
    hold on
    subplot(3,1,3);
    plot(t,y(:,4))
    hold on
    figure(2)
    set(gca,'fontsize',24);
    plot(t,y(:,1))
    hold on
    figure(3)
    set(gca,'fontsize',24);
    plot(y(:,1),y(:,2));
    hold on
end

%% Q17

fprintf('\nQ17\n');

for k=[0,0.1,0.17,0.2]
    para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,k];
    V=fzero(@Eq_V,-60);
    Eq(:,1)=Eq_Find_HH(V);
    
    j=2;i=1;
    while i>0
        V=fzero(@Eq_V,Eq(1,j-1)+60);
        Eq(:,j)=Eq_Find_HH(V);
        if Eq(1,j)-Eq(1,j-1)<1e-5
            i=0;
        end
        j=j+1;
    end
    
    Eq(:,1:end-1);
    Vss=Eq(1,1);
    y0=Ref_Eq+[10,0,0,0];
    [t,y] = ode15s(@HH_Vnr, tspan, y0);%, options);  %ODE solver
    
    figure(2)
    set(gca,'fontsize',24);
    plot(t,y(:,1))
    hold on
    figure(3)
    set(gca,'fontsize',24);
    plot(y(:,1),y(:,2));
    hold on
end

%% Q18

fprintf('\nQ18\n');

for k=[0,0.1,0.17,0.2]
    para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,k];
    
    y0=Ref_Eq+[10,0,0,0];
    [t,y] = ode15s(@HH_Vnr, tspan, y0);%, options);  %ODE solver
end


%% Q19

fprintf('\nQ19\n');

para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
iext = [-3, 100,120];
V=fzero(@Eq_V,-60);
Eq(:,1)=Eq_Find_HH(V)
y0=Ref_Eq+[1,0,0,0];
[t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver

figure
plot(t,y(:,1))
set(gca,'fontsize',24);
figure
plot(y(:,1),y(:,2));
set(gca,'fontsize',24);

%% Q20

fprintf('\nQ20\n');

para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
iext = [-3,100,120];
V=fzero(@Eq_V,-60);
Eq(:,1)=Eq_Find_HH(V)
y0=Ref_Eq+[1,0,0,0];
Vss=Eq(1,1);
[t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver

figure
plot(t,y(:,1))
set(gca,'fontsize',24);
figure
plot(y(:,1),y(:,3));
set(gca,'fontsize',24);


%% Functions

function [x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w]=Setup()

global para iext

x1=-80:1:80;    %x-axis of phase plane plot
I=1;
for i=x1
    w_nc(I)=0.5*(1+tanh((i-para(10))/para(11)));    %w nullcline
    V_nc(I)=(iext(1)-para(1)*(0.5*(1+tanh((i-para(8))/para(9))))*(i-para(4))-para(3)*(i-para(6)))/(para(2)*(i-para(5)));    %V nullcline
    I=I+1;
end

Eq(1)=fzero(@Eq_Find,0); %Finding the equilibrium point. Fzero function returns the value where the function is zero.
Eq(2)= 0.5*(1+tanh((Eq(1)-para(10))/para(11))); %Finding the value of w at equilibrium

jac=JAC(Eq(1:2));    %fidning the jacobian matrix at equilibrium
lambda=eig(jac); %finding the eigenvalues at equilibrium

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
        Grad_V(I,J) = (iext(1)-para(1)*ms*(i-para(4)) - para(2)*j*(i-para(5)) - para(3)*(i-para(6)))/para(14);   %finding the gradinet of V in the mesh defined
        Grad_w(I,J) = para(7)*(ws-j)/tw; %finding the gradinet of w in the mesh defined
    end
    J=0;
end
end

function Output(x1,w_nc,V_nc,Eq,lambda,mx2,my2,Grad_V,Grad_w,tspan,t,y)

fprintf('The Equilibrium Point is at (V,w)=(%g,%g)\n',Eq(1),Eq(2));
fprintf('The Eigenvalues are:%g+%gi,%g+%gi)\n',real(lambda(1)),imag(lambda(1)),real(lambda(2)),imag(lambda(2)));
SN=0;
if imag(lambda)==[0,0]
    if sign(lambda)==[1,1]
        disp('The Equilibrium point is unstable')
    elseif sign(lambda)==[-1,-1]
        disp('The Equilibrium point is stable')
    else
        disp('The Equilibrium point is a saddle node') %Also includes the crital cases. Be careful
        MF10=[Eq(1)-0.01,Eq(2)];
        [t_MF,MF1] = ode15s(@MLE_b, tspan, MF10);%, options);  %ODE solver
        MF20=[Eq(1)+0.01,Eq(2)];
        [t_MF,MF2] = ode15s(@MLE_b, tspan, MF20);%, options);  %ODE solver
        MF30=[Eq(1)-0.01,Eq(2)];
        [t_MF,MF3] = ode15s(@MLE, tspan, MF30);%, options);  %ODE solver
        MF40=[Eq(1)+0.01,Eq(2)];
        [t_MF,MF4] = ode15s(@MLE, tspan, MF40);%, options);  %ODE solver
        SN=1;
    end
elseif sign(real(lambda))==[1,1]
    disp('The Equilibrium point is spiral outwards')
else
    disp('The Equilibrium point is spiral inwards') %Also includes the crital cases. Be careful
end

figure;
[ax, h1, h2] = plotyy(t, y(:,1), t, y(:,2));
axes(ax(1)); axis([0 tspan(2) -80 80]); ylabel('V, mV')
axes(ax(2)); axis([0 tspan(2) 0 0.5]); ylabel('w')
xlabel('Time, ms');

figure; %the values of w are multiplied by 100 in the plots, to give a visible output from the quiver function.
hold on
plot(x1, w_nc*100,'r--');
plot(x1, V_nc*100,'g--');
quiver(mx2, my2*100, Grad_V, Grad_w*100,'b');
% plot(y(:,1), y(:,2)*100);
quiver(y(:,1), y(:,2)*100, gradient(y(:,1)), gradient(y(:,2))*100,'color','m')
if SN==1
    quiver(MF1(:,1), MF1(:,2)*100, -gradient(MF1(:,1)), -gradient(MF1(:,2))*100,'color','k')
    quiver(MF2(:,1), MF2(:,2)*100, -gradient(MF2(:,1)), -gradient(MF2(:,2))*100,'color','k')
    quiver(MF3(:,1), MF3(:,2)*100, gradient(MF3(:,1)), gradient(MF3(:,2))*100,'color','k')
    quiver(MF4(:,1), MF4(:,2)*100, gradient(MF4(:,1)), gradient(MF4(:,2))*100,'color','k')
end
plot(Eq(1),Eq(2)*100,'bo');
text(Eq(1),Eq(2)*100-2,['(' num2str(Eq(1)) ',' num2str(Eq(2)*100) ')']);
% legend('w nullcline','V nullcline','Gradient','Trajectory','Equilibrium')
xlabel('V, mV');ylabel('w x100')
axis([-80,80,-20,100])
hold off

end

function Eq = Eq_Find(x)

global para iext

w_nc=0.5*(1+tanh((x-para(10))/para(11)));
V_nc=(iext(1)-para(1)*(0.5*(1+tanh((x-para(8))/para(9))))*(x-para(4))-para(3)*(x-para(6)))/(para(2)*(x-para(5)));
Eq=w_nc-V_nc;   %finding the intersection of both the nullclines

end

function dydt = MLE(t,y)

global para iext

tw = 1./cosh((y(1)-para(12))/(2*para(13)));
ms = 0.5*(1 + tanh((y(1)-para(8))/para(9)));
ws = 0.5*(1 + tanh((y(1)-para(10))/para(11)));

dydt = zeros(2,1);
if t>=iext(2) && t<iext(3)
    dydt(1) = (iext(1) - para(1)*ms*(y(1)-para(4)) - para(2)*y(2)*(y(1)-para(5)) - para(3)*(y(1)-para(6)))/para(14);
else
    dydt(1) = (-para(1)*ms*(y(1)-para(4)) -  para(2)*y(2)*(y(1)-para(5)) - para(3)*(y(1)-para(6)))/para(14);
end
dydt(2) = para(7)*(ws-y(2))./tw;
end

function dydt = MLE_b(t,y)

global para iext

tw = 1./cosh((y(1)-para(12))/(2*para(13)));
ms = 0.5*(1 + tanh((y(1)-para(8))/para(9)));
ws = 0.5*(1 + tanh((y(1)-para(10))/para(11)));

dydt = zeros(2,1);
if t>=iext(2) && t<iext(3)
    dydt(1) = -(iext(1) - para(1)*ms*(y(1)-para(4)) - para(2)*y(2)*(y(1)-para(5)) - para(3)*(y(1)-para(6)))/para(14);
else
    dydt(1) = -(-para(1)*ms*(y(1)-para(4)) -  para(2)*y(2)*(y(1)-para(5)) - para(3)*(y(1)-para(6)))/para(14);
end
dydt(2) = -para(7)*(ws-y(2))./tw;
end

function jac = JAC(y)

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



function Eq=Eq_Find_HH(V)

global para

if abs(V+50)<1e-5
    an=para(7)/10;
else
    an=-0.01*para(7)*(V+50)/(exp(-(V+50)/10)-1);
end
bn=0.125*para(7)*exp(-(V+60)/80);

if abs(V+35)<1e-5
    am=para(7);
else
    am=-0.1*para(7)*(V+35)/(exp(-(V+35)/10)-1);
end
bm=4*para(7)*exp(-(V+60)/18);

ah=0.07*para(7)*exp(-(V+60)/20);
bh=para(7)/(exp(-(V+30)/10)+1);

Eq(1)=V;
Eq(2)=an/(an+bn);
Eq(3)=am/(am+bm);
Eq(4)=ah/(ah+bh);

end

function Eq=Eq_V(V)

global para iext

if abs(V+50)<1e-5
    an=para(7)/10;
else
    an=-0.01*para(7)*(V+50)/(exp(-(V+50)/10)-1);
end
bn=0.125*para(7)*exp(-(V+60)/80);

if abs(V+35)<1e-5
    am=para(7);
else
    am=-0.1*para(7)*(V+35)/(exp(-(V+35)/10)-1);
end
bm=4*para(7)*exp(-(V+60)/18);

ah=0.07*para(7)*exp(-(V+60)/20);
bh=para(7)/(exp(-(V+30)/10)+1);

n=an/(an+bn);
m=am/(am+bm);
h=ah/(ah+bh);

Eq=(iext(1) - para(1)*(n^4)*(V-para(2)) - para(3)*((1-para(9))*(m^3)*h+para(9)*(m^3))*(V-para(4)) - para(5)*(V-para(6)))/para(8);
end

function dydt = HH(t,y)

global para iext

if abs(y(1)+50)<1e-5
    an=para(7)/10;
else
    an=-0.01*para(7)*(y(1)+50)/(exp(-(y(1)+50)/10)-1);
end
bn=0.125*para(7)*exp(-(y(1)+60)/80);

if abs(y(1)+35)<1e-5
    am=para(7);
else
    am=-0.1*para(7)*(y(1)+35)/(exp(-(y(1)+35)/10)-1);
end
bm=4*para(7)*exp(-(y(1)+60)/18);

ah=0.07*para(7)*exp(-(y(1)+60)/20);
bh=para(7)/(exp(-(y(1)+30)/10)+1);

dydt = zeros(4,1);

if t>=iext(2) && t<iext(3)
    dydt(1) = (iext(1) - para(1)*(y(2)^4)*(y(1)-para(2)) - para(3)*((1-para(9))*(y(3)^3)*y(4)+para(9)*(y(3)^3))*(y(1)-para(4)) - para(5)*(y(1)-para(6)))/para(8);
else
    dydt(1) = (- para(1)*(y(2)^4)*(y(1)-para(2)) - para(3)*((1-para(9))*(y(3)^3)*y(4)+para(9)*(y(3)^3))*(y(1)-para(4)) - para(5)*(y(1)-para(6)))/para(8);
end
dydt(2) = an*(1-y(2))-bn*y(2);
dydt(3) = am*(1-y(3))-bm*y(3);
dydt(4) = ah*(1-y(4))-bh*y(4);
end


function dydt = HH_Vnr(t,y)

global para iext Vss

if abs(y(1)+50)<1e-5
    an=para(7)/10;
else
    an=-0.01*para(7)*(y(1)+50)/(exp(-(y(1)+50)/10)-1);
end
bn=0.125*para(7)*exp(-(y(1)+60)/80);

if abs(y(1)+35)<1e-5
    am=para(7);
else
    am=-0.1*para(7)*(y(1)+35)/(exp(-(y(1)+35)/10)-1);
end
bm=4*para(7)*exp(-(y(1)+60)/18);

ah=0.07*para(7)*exp(-(Vss+60)/20);
bh=para(7)/(exp(-(Vss+30)/10)+1);

y(3)=am/(am+bm);
y(4)=ah/(ah+bh);

dydt = zeros(4,1);

if t>=iext(2) && t<iext(3)
    dydt(1) = (iext(1) - para(1)*(y(2)^4)*(y(1)-para(2)) - para(3)*((1-para(9))*(y(3)^3)*y(4)+para(9)*(y(3)^3))*(y(1)-para(4)) - para(5)*(y(1)-para(6)))/para(8);
else
    dydt(1) = (- para(1)*(y(2)^4)*(y(1)-para(2)) - para(3)*((1-para(9))*(y(3)^3)*y(4)+para(9)*(y(3)^3))*(y(1)-para(4)) - para(5)*(y(1)-para(6)))/para(8);
end
dydt(2) = an*(1-y(2))-bn*y(2);
dydt(3) = 0;
dydt(4) = 0;
end


function dydt = HH_Vmr(t,y)

global para iext Vss

if abs(Vss+50)<1e-5
    an=para(7)/10;
else
    an=-0.01*para(7)*(Vss+50)/(exp(-(Vss+50)/10)-1);
end
bn=0.125*para(7)*exp(-(Vss+60)/80);

if abs(y(1)+35)<1e-5
    am=para(7);
else
    am=-0.1*para(7)*(y(1)+35)/(exp(-(y(1)+35)/10)-1);
end
bm=4*para(7)*exp(-(y(1)+60)/18);

ah=0.07*para(7)*exp(-(Vss+60)/20);
bh=para(7)/(exp(-(Vss+30)/10)+1);

y(2)=an/(an+bn);
y(4)=ah/(ah+bh);

dydt = zeros(4,1);

if t>=iext(2) && t<iext(3)
    dydt(1) = (iext(1) - para(1)*(y(2)^4)*(y(1)-para(2)) - para(3)*((1-para(9))*(y(3)^3)*y(4)+para(9)*(y(3)^3))*(y(1)-para(4)) - para(5)*(y(1)-para(6)))/para(8);
else
    dydt(1) = (- para(1)*(y(2)^4)*(y(1)-para(2)) - para(3)*((1-para(9))*(y(3)^3)*y(4)+para(9)*(y(3)^3))*(y(1)-para(4)) - para(5)*(y(1)-para(6)))/para(8);
end
dydt(2) = 0;
dydt(3) = am*(1-y(2))-bm*y(2);
dydt(4) = 0;
end
