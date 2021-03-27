clear all
close all
format long
%% Setup Inputs
global para iext
%para=[Gkb(1),Ek(2),GNab(3),ENa(4),GLb(5),EL(6),phi(7),C(8),fni(9)];
%y=[V,n,m,h];
%iext = [i_val, i_Tstart, i_Tstop]
tspan = [0, 500]; %Tstart, Tstop
T=6.3;
para=[36,-72,120,55,0.3,-49.4010790306083,3^((T-6.3)/10),1,0];
iext = [0, 0,500];
Ref_Eq=[-60,0.317676914060697,0.052932485257250,0.596120753508460];

% para(9)=0.17;

V=fzero(@Eq_V,-60);
Eq(:,1)=Eq_Find(V);

j=2;i=1;
while i>0
    V=fzero(@Eq_V,Eq(1,j-1)+60);
    Eq(:,j)=Eq_Find(V);
    if Eq(1,j)-Eq(1,j-1)<1e-5
        i=0;
    end
    j=j+1;
end

global Vss
Vss=Eq(1,1);
Eq(:,1:end-1)
y0=Eq(:,1)+[0,0,0,0]';
[t,y] = ode15s(@HH, tspan, y0);%, options);  %ODE solver

figure
subplot(3,1,1);
plot(t,y(:,2))
subplot(3,1,2);
plot(t,y(:,3))
subplot(3,1,3);
plot(t,y(:,4))
figure
plot(t,y(:,1))

figure
plot(y(:,1),y(:,2));    
%% Functions

function Eq=Eq_Find(V)

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

if abs(y(1)+50)<1e-5
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
