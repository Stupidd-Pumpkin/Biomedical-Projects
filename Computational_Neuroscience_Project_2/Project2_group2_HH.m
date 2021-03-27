close all; clear all; clc
global Iext
global Gk
global Gl 
global Gna
global Vna
global Vk
global Vl
global C
global fni
global h_eqbm
global n_eqbm
global m_eqbm

Gk = 36; %mS/cm^2
Vk = -72; %mV
Gna = 120; %mS/cm^2
Vna = 55; %mV
Gl = 0.3; %mS/cm^2
C = 1; %muF/cm^2
Iext=0;
fni =0;


%% Question 13

x(1) = -60;
an=-0.01*(x(1)+50)/(exp(-(x(1)+50)/10)-1);
bn=0.125*exp(-(x(1)+60)/80);
am=-0.1*(x(1)+35)/(exp(-(x(1)+35)/10)-1);
bm=4*exp(-(x(1)+60)/18);
ah=0.07*exp(-(x(1)+60)/20);
bh=1/(exp(-(x(1)+30)/10)+1);
x(2) = an/(an+bn);
x(3) = am/(am+bm);
x(4) = ah/(ah+bh);
g_K = Gk*(x(2)^4);
g_Na = Gna*(x(3))^3*(x(4));
g_L = Gl;
Vl = x(1)-(Iext - g_K*(x(1)-Vk)-g_Na*(x(1)-Vna))/g_L

V_eqbm = x(1) 
n_eqbm = x(2) 
m_eqbm = x(3) 
h_eqbm = x(4)

Iext = 10;
tspan = [1 200];
y0 = [V_eqbm n_eqbm m_eqbm h_eqbm];
[t,y] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
figure
plot(t,y(:,1))
title('Spiking at Iext = 10muA/cm^2')
xlabel('time')
ylabel('Voltage in mV')



%% Question 14
%Stability at Iext =0
for i = 0
    Iext = i;
    y_init = [-60; 0;0;0];
    equiPoint3 = fsolve(@hh,y_init);
    %plot(equiPoint3(1),100*equiPoint3(3),'b*')
    disp(equiPoint3)

    syms V m
    J = jacobian([ (Iext -(Gna*m^3*((0.07*(exp(-(V+60)/20)))/((0.07*(exp(-(V+60)/20)))+ (1/(exp(-(V+30)/10)+1))))*(V-Vna)) - Gk*((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))/((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))+ (0.125*exp(-(V+60)/80))))^4*(V-Vk) - Gl*(V-Vl))/C, (-0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12))*(1-m) - (4*exp(-(V+60)/18))*m],[V,m]);
    equi = zeros(1,2);
    equi(1) = equiPoint3(1);
    equi(2) = equiPoint3(3);
    A = double(subs(J,[V,m],equi));
    [V,D] = eig(A)
    if D(1,1) < 0 && D(2,2) < 0
        display(i)
        display("Stable")
    else
        display(i)
        display("Unstable")
    end
    
y0 = [equiPoint3(1)+0.1 equiPoint3(2) equiPoint3(3) equiPoint3(4)];
tspan = [0 50];
Iext = i;
[t,y] = ode15s(@(t,y) hh_ode(t,y), tspan, y0);
figure
subplot(2,1,1)
plot(t,y(:,1));
xlabel('Time t');
ylabel('Voltage V');
axis([0 50 -100 100]);
subplot(2,1,2)
plot(t,y(:,2),t,y(:,3),t,y(:,4));
title('Variation of m,n,h with time');
xlabel('Time t');
ylabel('x');
legend('n','m','h');
        
end    
% Checking the threshold
Iext = 0;
tspan = [1 50];
y0 = [-50 n_eqbm m_eqbm h_eqbm];
[t1,y1] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
y0 = [-45 n_eqbm m_eqbm h_eqbm];
[t2,y2] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
y0 = [-40 n_eqbm m_eqbm h_eqbm];
[t3,y3] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
y0 = [-49 n_eqbm m_eqbm h_eqbm];
[t4,y4] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
y0 = [-49.9 n_eqbm m_eqbm h_eqbm];
[t5,y5] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
y0 = [-49.999 n_eqbm m_eqbm h_eqbm];
[t6,y6] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);

figure
plot(t1,y1(:,1)); hold on; plot(t2,y2(:,1)); hold on; plot(t3,y3(:,1)); hold on; plot(t4,y4(:,1));hold on; plot(t5,y5(:,1));hold on; plot(t6,y6(:,1))
title('Checking the threshold at Iext = 0muA/cm^2')
xlabel('time')
ylabel('Voltage in mV')
legend('Vi = -50', 'Vi = -45', 'Vi = -40','Vi = -49','Vi = -49.9','Vi = -49.999')

%% Question 15
for i = 8:1:12
    Iext = i;
    y_init = [-60; 0;0;0];
    equiPoint3 = fsolve(@hh,y_init);
    %plot(equiPoint3(1),100*equiPoint3(3),'b*')
    disp(equiPoint3)

    syms V m
    J = jacobian([ (Iext -(Gna*m^3*((0.07*(exp(-(V+60)/20)))/((0.07*(exp(-(V+60)/20)))+ (1/(exp(-(V+30)/10)+1))))*(V-Vna)) - Gk*((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))/((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))+ (0.125*exp(-(V+60)/80))))^4*(V-Vk) - Gl*(V-Vl))/C, (-0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12))*(1-m) - (4*exp(-(V+60)/18))*m],[V,m]);
    equi = zeros(1,2);
    equi(1) = equiPoint3(1);
    equi(2) = equiPoint3(3);
    A = double(subs(J,[V,m],equi));
    [V,D] = eig(A)
    if D(1,1) < 0 && D(2,2) < 0
        display(i)
        display("Stable")
    else
        display(i)
        display("Unstable")
    end
    y0 = [equiPoint3(1)+3 equiPoint3(2) equiPoint3(3) equiPoint3(4)];
tspan = [0 50];
Iext = i;
[t,y] = ode15s(@(t,y) hh_ode(t,y), tspan, y0);
figure
subplot(2,1,1)
plot(t,y(:,1));
xlabel('Time t');
ylabel('Voltage V');
axis([0 50 -100 100]);
subplot(2,1,2)
plot(t,y(:,2),t,y(:,3),t,y(:,4));
title('Variation of m,n,h with time');
xlabel('Time t');
ylabel('x');
legend('n','m','h');
        
end    
clearvars t1 t2 t3 t4 t5 t6 y1 y2 y3 y4 y5 y6 
%% Question 16
%Here, we have an additional parameter fi
%fni = 0, 0.1, 0.17, and 0.2
Iext=0;
fni = 0;
y_init = [-60; 0;0;0];
equiPointf = fsolve(@hh,y_init);
tspan = [1 30];
y0 = [7 equiPointf(2) equiPointf(3) equiPointf(4)];
[t1,y1] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
figure
subplot(2,1,1)
plot(t1,y1(:,1));
xlabel('Time t');
ylabel('Voltage V');
axis([0 50 -100 100]);
subplot(2,1,2)
plot(t1,y1(:,2),t1,y1(:,3),t1,y1(:,4));
title('Variation of m,n,h with time for fni = 0');
xlabel('Time t');
ylabel('x');
legend('n','m','h');

fni = 0.1;
y_init = [-60; 0;0;0];
equiPointf = fsolve(@hh,y_init);
tspan = [1 30];
y0 = [7 equiPointf(2) equiPointf(3) equiPointf(4)];
[t2,y2] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
figure
subplot(2,1,1)
plot(t2,y2(:,1));
xlabel('Time t');
ylabel('Voltage V');
axis([0 50 -100 100]);
subplot(2,1,2)
plot(t2,y2(:,2),t2,y2(:,3),t2,y2(:,4));
title('Variation of m,n,h with time for fni = 0.1');
xlabel('Time t');
ylabel('x');
legend('n','m','h');

fni = 0.17;
y_init = [-60; 0;0;0];
equiPointf = fsolve(@hh,y_init);
tspan = [1 30];
y0 = [7 equiPointf(2) equiPointf(3) equiPointf(4)];
[t3,y3] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
figure
subplot(2,1,1)
plot(t3,y3(:,1));
xlabel('Time t');
ylabel('Voltage V');
axis([0 50 -100 100]);
subplot(2,1,2)
plot(t3,y3(:,2),t3,y3(:,3),t3,y3(:,4));
title('Variation of m,n,h with time for fni = 0.17');
xlabel('Time t');
ylabel('x');
legend('n','m','h');


fni = 0.2;
y_init = [-60; 0;0;0];
equiPointf = fsolve(@hh,y_init);
tspan = [1 30];
y0 = [7 equiPointf(2) equiPointf(3) equiPointf(4)];
[t4,y4] = ode15s(@(t,y) hh_ode(t,y),tspan, y0);
figure
subplot(2,1,1)
plot(t4,y4(:,1));
xlabel('Time t');
ylabel('Voltage V');
axis([0 50 -100 100]);
subplot(2,1,2)
plot(t4,y4(:,2),t4,y4(:,3),t4,y4(:,4));
title('Variation of m,n,h with time for fni = 0.2');
xlabel('Time t');
ylabel('x');
legend('n','m','h');

figure
plot(t1,y1(:,1)); hold on; plot(t2,y2(:,1)); hold on; plot(t3,y3(:,1)); hold on; plot(t4,y4(:,1))
title('Spiking behaviour for different values of fni')
xlabel('time in mS')
ylabel('Voltage in mV')
legend('fni = 0','fni= 0.1','fni= 0.17','fni= 0.2')

%% Question 17
clearvars t1 t2 t3 t4 t5 t6 y1 y2 y3 y4 y5 y6
Iext = 0;
tspan = [1 50];
y0 = [-50 n_eqbm m_eqbm h_eqbm];
[t1,y1] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-56 n_eqbm m_eqbm h_eqbm];
[t2,y2] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-57 n_eqbm m_eqbm h_eqbm];
[t3,y3] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-56.6 n_eqbm m_eqbm h_eqbm];
[t4,y4] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-56.65 n_eqbm m_eqbm h_eqbm];
[t5,y5] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-56.7 n_eqbm m_eqbm h_eqbm];
[t6,y6] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-56.669 n_eqbm m_eqbm h_eqbm];
[t7,y7] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
figure
plot(t1,y1(:,1)); hold on; plot(t2,y2(:,1)); hold on; plot(t3,y3(:,1)); hold on; plot(t4,y4(:,1));hold on; plot(t5,y5(:,1));hold on; plot(t6,y6(:,1));hold on; plot(t7,y7(:,1))
title('Checking the threshold at Iext = 0muA/cm^2 for V-n reduced system')
xlabel('time')
ylabel('Voltage in mV')
legend('Vi = -50', 'Vi = -56','Vi = -57','Vi = -56.6','Vi = -56.65','Vi = -56.7','Vi = -56.669')


y0 = [-40 n_eqbm m_eqbm h_eqbm];
[t8,y8] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;
y0 = [-20 n_eqbm m_eqbm h_eqbm];
[t9,y9] = ode15s(@(t,y) hh_odeVn(t,y,h_eqbm),tspan, y0);clearvars y0;Iext = 0;

figure
plot(y1(:,1),y1(:,2)); hold on; plot(y2(:,1),y2(:,2)); plot(y3(:,1),y3(:,2)); hold on;plot(y4(:,1),y4(:,2)); hold on;plot(y5(:,1),y5(:,2)); hold on;plot(y6(:,1),y6(:,2)); hold on;plot(y7(:,1),y7(:,2));hold on;plot(y8(:,1),y8(:,2));hold on;plot(y9(:,1),y9(:,2));
title('Various trajectories in V-n plane at I_ext =0')
xlabel('V in mV')
ylabel('n')
legend('Vi = -50', 'Vi = -56','Vi = -57','Vi = -56.6','Vi = -56.65','Vi = -56.7','Vi = -56.669','Vi = -40','Vi=-20')

%%% plotting the nullclines
% V= linspace(-79,40,100); %in muV
% for i=1:length(V)
%     alphan(i) = -0.01*(V(i)+50)/(exp(-(V(i)+50)/10)-1);
%     if isnan(alphan(i))
%         alphan(i) = 1;
%     end
%     betan(i) = 0.125*exp(-(V(i)+60)/80);
%     alpham(i) = -0.1*(V(i)+35)/(exp(-(V(i)+35)/10)-1);
%     if isnan(alpham(i))
%         alpham(i) = 1;
%     end
%     betam(i) = 4*exp(-(V(i)+60)/18);
%     alphah(i) = 0.07*(exp(-(V(i)+60)/20));
%     betah(i) = 1/(exp(-(V(i)+30)/10)+1);
%     
%     m(i) = alpham(i)/(alpham(i)+betam(i));
%     h(i) = h_eqbm;
%     dndt = alphan(i)-alphan(i)*n(i) - betan(i)*n(i);
%     dndt = 
%     
%  figure
% plot(V_stim,w_ncw,V_stim,w_ncv);
% title(' W and V nullcline plot');
% xlabel('V');
% ylabel('W');
% legend('V nullcline','W nullcline')

%for some Iext
Iext = 10;
y_init = [-60; 0;0;0];
equiPoint10 = fsolve(@hh,y_init);
display(equiPoint10)
tspan = [1 50];
nnn = equiPoint10(2); mmm = equiPoint10(3); hhh=equiPoint10(4);
y0 = [-40 nnn mmm hhh];
[t10,y10] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
figure
plot(t10, y10(:,1))
title(' Spiking at Iext = 10 for V-n reduced system')
xlabel('time in ms')
ylabel('Voltage in mV')

Iext = 10;
y0 = [-55 nnn mmm hhh];
[t11,y11] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-60 nnn mmm hhh];
[t12,y12] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-20 nnn mmm hhh];
[t13,y13] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 


figure 
plot(y10(:,1), y10(:,2)); hold on; plot(y11(:,1),y11(:,2)); hold on; plot(y12(:,1),y12(:,2));hold on; plot(y13(:,1),y13(:,2))
title('Phase plane analysis of V-n reduced system for Iext = 10')
xlabel('V in mV')
ylabel('n')
legend('Vi=-40', 'Vi=-55','Vi = -60','Vi=-20')

%%% For Iext =8
Iext = 8;
y_init = [-60; 0;0;0];
equiPoint10 = fsolve(@hh,y_init);
display(equiPoint10)
tspan = [1 50];
nnn = equiPoint10(2); mmm = equiPoint10(3); hhh=equiPoint10(4);
y0 = [-40 nnn mmm hhh];
[t10,y10] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
figure
plot(t10, y10(:,1))
title(' Spiking at Iext = 8 for V-n reduced system')
xlabel('time in ms')
ylabel('Voltage in mV')
Iext = 8;
y0 = [-55 nnn mmm hhh];
[t11,y11] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-60 nnn mmm hhh];
[t12,y12] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-20 nnn mmm hhh];
[t13,y13] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
figure 
plot(y10(:,1), y10(:,2)); hold on; plot(y11(:,1),y11(:,2)); hold on; plot(y12(:,1),y12(:,2));hold on; plot(y13(:,1),y13(:,2))
title('Phase plane analysis of V-n reduced system for Iext = 8')
xlabel('V in mV')
ylabel('n')
legend('Vi=-40', 'Vi=-55','Vi = -60','Vi=-20')

%%% For Iext = 12
Iext = 12;
y_init = [-60; 0;0;0];
equiPoint10 = fsolve(@hh,y_init);
display(equiPoint10)
tspan = [1 50];
nnn = equiPoint10(2); mmm = equiPoint10(3); hhh=equiPoint10(4);
y0 = [-40 .7 mmm hhh];
[t10,y10] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
figure
plot(t10, y10(:,1))
title(' Spiking at Iext = 12 for V-n reduced system')
xlabel('time in ms')
ylabel('Voltage in mV')
Iext = 12;
y0 = [-55 nnn mmm hhh];
[t11,y11] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-60 nnn mmm hhh];
[t12,y12] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-20 0.4 mmm hhh];
[t13,y13] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0; 
y0 = [-60 0.5 mmm hhh];
[t14,y14] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0;
y0 = [-30 0.3 mmm hhh];
[t15,y15] = ode15s(@(t,y) hh_odeVn(t,y,hhh),tspan, y0);clearvars y0;
figure 
plot(y10(:,1), y10(:,2)); hold on; plot(y11(:,1),y11(:,2)); hold on; plot(y12(:,1),y12(:,2));hold on; plot(y13(:,1),y13(:,2));hold on; plot(y14(:,1),y14(:,2));hold on; plot(y15(:,1),y15(:,2))
title('Phase plane analysis of V-n reduced system for Iext = 12')
xlabel('V in mV')
ylabel('n')
legend('Vi=-40', 'Vi=-55','Vi = -60','Vi=-20','vi=-60', 'Vi=-30')

%% Question 18
alpham = @(V) -0.1 * (35+V) ./ (exp(-0.1*(35+V)) - 1);
betam = @(V) 4 * exp(-(60+V)/18);

alphah = @(V) 0.07 * exp(-(60+V)/20);
betah = @(V) 1 ./ (exp(-(30+V)/10) + 1);
 

alphan = @(V) 0.01 * (-(50+V)) ./ (exp(-0.1*(50+V)) - 1);
betan = @(V) 0.125 * exp(-(60+V)/80);

figure;
     Iext=0;
fni_array=[0.02 0.07 0.1 0.18 0.2 0.25 0.3 0.35 0.4];
%Ina=gNa*(1-fni)*y(2).^3.*y(4)*(y(1)-VNa)+gNa*(fni)*y(2).^3.*(y(1)-VNa);
for i=1:9
    fni=fni_array(i);
    myfun1 = @(V,n) (-Gk*n.^4.*(V - Vk) - Gna*(alpham(V)./(alpham(V)+betam(V))).^3.*fni.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
    myfun2 = @(V,n) alphan(V).*(1-n)-betan(V).*n;
    
    a1=ezplot(@(V,n) myfun1(V,n), [-72 120 0 1])
    set(a1,'Color','blue', 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    set(gca, 'fontsize', 10)
    a2=ezplot(@(V,n) myfun2(V,n), [-72 120 0 1])
    title('Phase plane for varying fni between 0.02 to 0.4');
    hold on;
%     f_nv = @(t,y) [ (-gK*y(2).^4.*(y(1) - VK) - gNa*(fni)*(am(y(1))./(am(y(1))+bm(y(1)))).^3.*(y(1)-VNa)-gL.*(y(1)-VL) + I)/C; ...
%                an(y(1)).*(1-y(2))-bn(y(1)).*y(2)]; ...  
%     hold on;
%     [Tnv, Ynv] = ode15s(f_nv, [0 10], [-53 0.317]);
%     plot(Ynv(:,1),Ynv(:,2));
%     axis([-72 120 0 1])
%     hold on;
%     title('phase plane plot');
end


%% Question 19
Iext = 0;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               alpham(y(1)).*(1-y(2))-betam(y(1)).*y(2); ...
               alphan(y(1)).*(1-y(3))-betan(y(1)).*y(3); ...
               alphah(y(1)).*(1-y(4))-betah(y(1)).*y(4) ];
[T1,Y1]=ode15s(f, [0 100], [V_eqbm n_eqbm m_eqbm h_eqbm]);

Iext = -3;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               alpham(y(1)).*(1-y(2))-betam(y(1)).*y(2); ...
               alphan(y(1)).*(1-y(3))-betan(y(1)).*y(3); ...
               alphah(y(1)).*(1-y(4))-betah(y(1)).*y(4) ];
[T2,Y2]=ode15s(f, [T1(length(T1)) 20 + T1(length(T1))], [Y1(length(Y1),1),Y1(length(Y1),2),Y1(length(Y1),3),Y1(length(Y1),4)]);

Iext = 0;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               alpham(y(1)).*(1-y(2))-betam(y(1)).*y(2); ...
               alphan(y(1)).*(1-y(3))-betan(y(1)).*y(3); ...
               alphah(y(1)).*(1-y(4))-betah(y(1)).*y(4) ];
[T3,Y3]=ode15s(f, [T2(length(T2)) 150 + T2(length(T2))], [Y2(length(Y2),1),Y2(length(Y2),2),Y2(length(Y2),3),Y2(length(Y2),4)]);
figure
plot([T1;T2;T3], [Y1(:,1); Y2(:,1); Y3(:,1)])
xlabel('time in ms'); ylabel('V in mV');
title('Membrane potential vs time');

figure
subplot(3,1,1)
plot([T1;T2;T3], [Y1(:,2); Y2(:,2); Y3(:,2)])
xlabel('time in ms'); ylabel('V in mV');
title('n vs time');
subplot(3,1,2)
plot([T1;T2;T3], [Y1(:,3); Y2(:,3); Y3(:,3)])
xlabel('time in ms'); ylabel('V in mV');
title('m vs time');
subplot(3,1,3)
plot([T1;T2;T3], [Y1(:,4); Y2(:,4); Y3(:,4)])
xlabel('time in ms'); ylabel('V in mV');
title('h vs time');


%% Question 20
Iext=0;
n = n_eqbm;
h = h_eqbm;
myfun1 = @(V,m) (-Gk*n.^4.*(V - Vk) - Gna*m.^3.*h.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
myfun2 = @(V,m) alpham(V).*(1-m)-betam(V).*m;
figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
title('Phase plane for n and h fixed at rest values');
xlabel('V in mv')
ylabel('m')

y0 = [-55 n_eqbm m_eqbm h_eqbm];
[t1,y1] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [-60 n_eqbm m_eqbm h_eqbm];
[t2,y2] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [-20 n_eqbm m_eqbm h_eqbm];
[t3,y3] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [-40 n_eqbm m_eqbm h_eqbm];
[t4,y4] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [0 n_eqbm m_eqbm h_eqbm];
[t5,y5] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);

figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
hold on;
plot(y1(:,1),y1(:,3)); hold on;
plot(y2(:,1),y2(:,3)); hold on;
plot(y3(:,1),y3(:,3)); hold on;
plot(y4(:,1),y4(:,3)); hold on;
plot(y5(:,1),y5(:,3)); 
title('Phase plane for n and h fixed at rest values');
legend('V nullcline','m nullcline', 'Vi = -55', 'Vi = -60' , 'Vi = -20','Vi = -40', 'Vi = 0')
xlabel('V in mv')
ylabel('m')



n = Y2(length(Y2),3);
h = Y2(length(Y2),4);
myfun1 = @(V,m) (-Gk*n.^4.*(V - Vk) - Gna*m.^3.*h.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
myfun2 = @(V,m) alpham(V).*(1-m)-betam(V).*m;
figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
title('Phase plane with trajectories for n and h fixed at end of stimulus values');
xlabel('V in mv')
ylabel('m')

y0 = [-55 n m_eqbm h];
[t1,y1] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [-60 n m_eqbm h];
[t2,y2] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [-20 n m_eqbm h];
[t3,y3] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [-40 n m_eqbm h];
[t4,y4] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);
y0 = [0 n m_eqbm h];
[t5,y5] = ode15s(@(t,y) hh_odeVm(t,y),tspan, y0);

figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
hold on;
plot(y1(:,1),y1(:,3)); hold on;
plot(y2(:,1),y2(:,3)); hold on;
plot(y3(:,1),y3(:,3)); hold on;
plot(y4(:,1),y4(:,3)); hold on;
plot(y5(:,1),y5(:,3)); 
title('Phase plane with trajectories for n and h fixed at end of stimulus values');
legend('V nullcline','m nullcline', 'Vi = -55', 'Vi = -60' , 'Vi = -20','Vi = -40', 'Vi = 0')
xlabel('V in mv')
ylabel('m')



%% Functions used
function dYdt = hh(y)
    global Iext
    global Gk
    global Gl 
    global Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*(1-fni)*m^3*h*(V-Vna)-Gna*fni*m^3*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt];
end

function dYdt = hh_ode(t,y)
    global Iext
    global Gk
    global Gl 
    global Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*(1-fni)*m^3*h*(V-Vna)-Gna*fni*m^3*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end




function dYdt = hh_odeVn(t,y,hh)
    global Iext
    global Gk
    global Gl 
    global Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    global h_eqbm
    
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    m = alpham/(alpham+betam);
    h = hh;
    dVdt = (Iext -Gna*(1-fni)*m^3*h*(V-Vna)-Gna*fni*m^3*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = 0;
    dndt = alphan*(1-n) - betan*n;
    dhdt = 0;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end


function dYdt = hh_odeVm(t,y)
    global Iext
    global Gk
    global Gl 
    global Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    global n_eqbm
    global h_eqbm
    
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    n = n_eqbm;
    h = h_eqbm;
    dVdt = (Iext -Gna*(1-fni)*m^3*h*(V-Vna)-Gna*fni*m^3*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = 0;
    dhdt = 0;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end

