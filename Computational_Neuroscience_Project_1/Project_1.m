clear all
close all

%Please Wait for about a minute for the program to run completely

%t_ODE45 contains the time taken to run ODE45 
%and similarly t_ODE15s for ODE15s for u=0.1,1,100 in order

%%
i=1;
u=0.1;      %value of u is set to 0.1
tspan = [0:0.1:100]; %tspan is selected appropriately for the u
y0 =[1;0];  %y0 contains the initial values of the state variable y(1) and y2)

tic      %tic toc is a function used to measure time taken by ODE45 and ODE15s functions
[t,y] = ode45(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0); %Van Der Pol Equation
t_ODE45(i)=toc;

figure(1)
plot(t,y);
title('y1,y2 vs t (with ODE45) for u=0.1');
legend('y1','y2');
figure(2)
plot(y(:,1),y(:,2))
title('y2 vs y1 (with ODE45) for u=0.1');

tic
[t,y] = ode15s(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0); 
t_ODE15s(i)=toc;

%%Figures of ODE15s can be observed below but they are the same as ODE45
% figure(3)
% plot(t,y);
% title('y1,y2 vs t (with ODE15s) for u=0.1');
% legend('y1', 'y2');
% figure(4)
% plot(y(:,1),y(:,2))
% title('y2 vs y1 (with ODE15s) for u=0.1');

%%
i=2;
u=1;
tspan = [0:0.1:100];
y0 =[1;0];

tic
[t,y] = ode45(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0); 
t_ODE45(i)=toc;

figure(5)
plot(t,y);
title('y1,y2 vs t (with ODE45) for u=1');
legend('y1','y2');
figure(6)
plot(y(:,1),y(:,2))
title('y2 vs y1 (with ODE45) for u=1');

tic
[t,y] = ode15s(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0); 
t_ODE15s(i)=toc;

%%Figures of ODE15s can be observed below but they are the same as ODE45
% figure(7)
% plot(t,y);
% title('y1,y2 vs t (with ODE15s) for u=1');
% legend('y1', 'y2');
% figure(8)
% plot(y(:,1),y(:,2))
% title('y2 vs y1 (with ODE15s) for u=1');


%%
i=3;
u=100;
tspan = [0:100:20000];
y0 =[1;0];

tic
[t,y] = ode45(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0); 
t_ODE45(i)=toc;

figure(9)
plot(t,y);
title('y1,y2 vs t (with ODE45) for u=100');
legend('y1','y2');
figure(10)
plot(y(:,1),y(:,2))
title('y2 vs y1 (with ODE45) for u=100');

tic
[t,y] = ode15s(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0); 
t_ODE15s(i)=toc;

%%Figures of ODE15s can be observed below but they are the same as ODE45
% figure(11)
% plot(t,y);
% title('y1,y2 vs t (with ODE15s) for u=100');
% legend('y1', 'y2');
% figure(12)
% plot(y(:,1),y(:,2))
% title('y2 vs y1 (with ODE15s) for u=100');
% i=i+1;

u=[0.1,1,100]
t_ODE45
t_ODE15s

%% Setup_2
% clear all
% close all
% 
% y0 =[1;0];
% 
% for u=[1]
%     tspan = [0:u/10:u*1000];
%     [t,y] = ode15s(@(t,y) [y(2); u*(1-y(1)^2)*y(2)-y(1)/u] ,tspan,y0);
% %     subplot(2,1,1);
% %     plot(t,y(:,1));
%     Y1=fftshift(fft(y(:,1)));
%     Y2=fftshift(fft(y(:,2)));
% %     subplot(2,1,2);
% %     plot(t,(y(:,1)),'o');
%     figure
%     plot(t,y(:,1))
% end
