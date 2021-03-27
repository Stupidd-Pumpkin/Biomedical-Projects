function [Vout]=Cap_DAC(Din,A,nob)

lsb=2*A/2^nob;
Qn=lsb^2/12;
Ctot = 1.38e-23*300*10/Qn;
Cu=Ctot/2^nob;
Cp=1*Cu;
Cu_arr=Cu*(1+0.001*randn(1,2^nob));
for i=1:length(Din)
    Vout(i)=2*A*sum(Cu_arr(1:Din(i)))/(sum(Cu_arr)+Cp);
end