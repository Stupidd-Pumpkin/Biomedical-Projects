clear all
close all
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

% %% Q14
% 
% fprintf('\nQ14\n');
% para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
% iext = [0, 100,125];
% V=fzero(@Eq_V,-60);
% Eq(:,1)=Eq_Find(V)
% y0=Ref_Eq+[0,0,0,0];
% [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
% 
% %% Q15
% 
% fprintf('\nQ15\n');
% para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
% 
% for k=[12.575757,12.575758]
%     iext = [k, tspan(1),tspan(2)];
%     V=fzero(@Eq_V,-60);
%     Eq(:,1)=Eq_Find(V)
%     y0=Eq(:,1)+[0.1,0,0,0]';
%     [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
%     
%     figure
%     plot(t,y(:,1))
% end
% 
% %% Q16
% 
% fprintf('\nQ16\n');
% 
% for k=[0,0.1,0.17,0.2]
%     para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,k];
%     V=fzero(@Eq_V,-60);
%     Eq(:,1)=Eq_Find(V);
%     
%     j=2;i=1;
%     while i>0
%         V=fzero(@Eq_V,Eq(1,j-1)+60);
%         Eq(:,j)=Eq_Find(V);
%         if Eq(1,j)-Eq(1,j-1)<1e-5
%             i=0;
%         end
%         j=j+1;
%     end
%     
%     Eq(:,1:end-1)
%     y0=Ref_Eq+[10,0,0,0];
%     [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
%     
%     figure(1)
%     set(gca,'fontsize',24);
%     subplot(3,1,1);
%     plot(t,y(:,2))
%     hold on
%     subplot(3,1,2);
%     plot(t,y(:,3))
%     hold on
%     subplot(3,1,3);
%     plot(t,y(:,4))
%     hold on
%     figure(2)
%     set(gca,'fontsize',24);
%     plot(t,y(:,1))
%     hold on
%     figure(3)
%     set(gca,'fontsize',24);
%     plot(y(:,1),y(:,2));
%     hold on
% end
% 
% %% Q17
% 
% fprintf('\nQ17\n');
% 
% for k=[0,0.1,0.17,0.2]
%     para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,k];
%     V=fzero(@Eq_V,-60);
%     Eq(:,1)=Eq_Find(V);
%     
%     j=2;i=1;
%     while i>0
%         V=fzero(@Eq_V,Eq(1,j-1)+60);
%         Eq(:,j)=Eq_Find(V);
%         if Eq(1,j)-Eq(1,j-1)<1e-5
%             i=0;
%         end
%         j=j+1;
%     end
%     
%     Eq(:,1:end-1);
%     Vss=Eq(1,1);
%     y0=Ref_Eq+[10,0,0,0];
%     [t,y] = ode15s(@HH_Vnr, tspan, y0);%, options);  %ODE solver
%     
%     figure(2)
%     set(gca,'fontsize',24);
%     plot(t,y(:,1))
%     hold on
%     figure(3)
%     set(gca,'fontsize',24);
%     plot(y(:,1),y(:,2));
%     hold on
% end
% 
% %% Q18
% 
% fprintf('\nQ18\n');
% 
% for k=[0,0.1,0.17,0.2]
%     para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,k];
%     
%     y0=Ref_Eq+[10,0,0,0];
%     [t,y] = ode15s(@HH_Vnr, tspan, y0);%, options);  %ODE solver
%     [x1,n_nc,V_nc,Eq,mx2,my2,Grad_V,Grad_n]=Setup_HH();
% end
%
% 
% %% Q19
% 
% fprintf('\nQ19\n');
% 
% para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
% iext = [-3, 100,120];
% V=fzero(@Eq_V,-60);
% Eq(:,1)=Eq_Find(V)
% y0=Ref_Eq+[1,0,0,0];
% [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
% 
% figure
% plot(t,y(:,1))
% set(gca,'fontsize',24);
% figure
% plot(y(:,1),y(:,2));
% set(gca,'fontsize',24);
% 
% %% Q20
% 
% fprintf('\nQ20\n');
% 
% para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
% iext = [-3,100,120];
% V=fzero(@Eq_V,-60);
% Eq(:,1)=Eq_Find(V)
% y0=Ref_Eq+[1,0,0,0];
% Vss=Eq(1,1);
% [t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver
% 
% figure
% plot(t,y(:,1))
% set(gca,'fontsize',24);
% figure
% plot(y(:,1),y(:,3));
% set(gca,'fontsize',24);

%% Functions

function [x1,n_nc,V_nc,Eq,mx2,my2,Grad_V,Grad_n]=Setup_HH()

global para iext

V=fzero(@Eq_V,-60);

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
%   
% j=2;i=1;
% while i>0
%     V=fzero(@Eq_V,Eq(1,j-1)+60);
%     Eq(:,j)=Eq_Find(V);
%     if Eq(1,j)-Eq(1,j-1)<1e-5
%         i=0;
%     end
%     j=j+1;
% end
% 
% Eq(:,1:end-1)
Vss=Eq(1,1);
    
x1=-70:1:70;    %x-axis of phase plane plot
I=1;
for i=x1
    if abs(i+50)<1e-5
        an=para(7)/10;
    else
        an=-0.01*para(7)*(i+50)/(exp(-(i+50)/10)-1);
    end
    bn=0.125*para(7)*exp(-(i+60)/80);
    n_nc(I)=an/(an+bn);   
    ms=am/(am+bm);
    V_nc(I)=abs(((iext(1) - para(3)*((1-para(9))*(ms^3)*Eq(4)+para(9)*(ms^3))*(i-para(4)) - para(5)*(i-para(6)))/(para(1)*(i-para(2))))^0.25);
    I=I+1;
end

x2=-70:10:70;
y2=0:0.05:0.7;
mx2=meshgrid(x2)';
my2=meshgrid(y2);
I=0;J=0;
for i=x2
    if abs(i+50)<1e-5
        an=para(7)/10;
    else
        an=-0.01*para(7)*(i+50)/(exp(-(i+50)/10)-1);
    end
    bn=0.125*para(7)*exp(-(i+60)/80);
    ms=am/(am+bm);
    I=I+1;
    for j=y2
        J=J+1;
        Grad_V(I,J) = (iext(1)-para(1)*(j^4)*(i-para(2)) - para(3)*((1-para(9))*(ms^3)*Eq(4)+para(9)*(ms^3))*(i-para(4)) -  para(5)*(i-para(6)))/para(8);
        Grad_n(I,J) = an*(1-j)-bn*j;
    end
    J=0;
end

figure; %the values of w are multiplied by 100 in the plots, to give a visible output from the quiver function.
hold on
plot(x1, n_nc*100,'r--');
plot(x1, V_nc*100,'g--');
quiver(mx2, my2*100, Grad_V, Grad_n*100,'b');
plot(Eq(1),Eq(2)*100,'bo');
plot_txt=text(Eq(1),Eq(2)*100-2,['(' num2str(Eq(1)) ',' num2str(Eq(2)*100) ')']);
plot_txt(1).FontSize = 20;
% legend('n nullcline','V nullcline','Gradient','Trajectory','Equilibrium')
xlabel('V, mV');ylabel('n x100')
axis([-70,70,-20,100])
set(gca,'fontsize',24);
hold off

end

function Eq=Eq_Find(V)

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
y(4)=1-y(2)/0.8;

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
