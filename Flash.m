function [T,B]=Flash(Vin,Vrp,Vrn,nob,flag_ideality)
% flag_ideality   =0 =>ideal
%                 =1 =>only comparator is non ideal
%                 =2 =>only resistor ladder has random mismatch
%                 =3 =>only resistor ladder gradient effect
%                 =4 =>all non idealities present

l=length(Vin);
lsb=(Vrp-Vrn)/2^nob;
rand_err=0.1;
grad_err=0.1;
switch flag_ideality
    case 0
        Vos=zeros(1,2^nob-1);
        res=ones(1,2^nob);
        step_size=(Vrp-Vrn)/sum(res);
        Vth=Vrn+step_size*(0:1:2^nob-1);
    case 1
        Vos=lsb*randn(1,2^nob-1);
        res=ones(1,2^nob);
        step_size=(Vrp-Vrn)/sum(res);
        Vth=Vrn+step_size*(0:1:2^nob-1);
    case 2
        Vos=lsb*randn(1,2^nob-1);
        res=ones(1,2^nob)-rand_err/2+rand_err*rand(1,2^nob);
        Vth=Vrn+(Vrp-Vrn)*cumsum(res)/sum(res);
    case 3
        Vos=lsb*randn(1,2^nob-1);
        res=ones(1,2^nob)+linspace(-grad_err,grad_err,2^nob);
        Vth=Vrn+(Vrp-Vrn)*cumsum(res)/sum(res);
    otherwise
        Vos=lsb*randn(1,2^nob-1);
        res=ones(1,2^nob)-rand_err/2+rand_err*rand(1,2^nob)+linspace(-grad_err,grad_err,2^nob);
        Vth=Vrn+(Vrp-Vrn)*cumsum(res)/sum(res);
end

T=zeros(l,2^nob-1);
B=zeros(l,1);
for i=1:l
    for j=1:2^nob-1
        T(i,j)= Comparator(Vin(i),Vth(j),Vos(j));
    end
    B(i)=sum(T(i,:));
end