function [B, Vx]=SAR(Vin,nob, A)

l=length(Vin);
B=zeros(1,nob);
Vx=zeros(1,nob);

lsb=2*A/2^nob; %V full swing = 2*A
Qn=lsb^2/12;
Ctot = 1.38e-23*300*10/Qn;
Cu=Ctot/2^nob;
Cp=1*Cu*(1+0.001*randn(1,1));
C_arr=Cu*(2.^(0:nob-1)+0.001*randn(1,nob));
B(nob)=1;
for i = 1:nob
    Vx(i)= (2*A*(dot(B,C_arr))-(sum(C_arr)+Cu+Cp)*Vin)/(sum(C_arr)+Cu+Cp);
    B(nob-i+1) = Comparator(0,Vx(i),0);
    if i <nob
        B(nob-i)=1;
    end
end
