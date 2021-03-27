function [D]=Comparator(Vin, Vth, Vos)
D=0;
if Vin > Vth+Vos
    D=1;
end