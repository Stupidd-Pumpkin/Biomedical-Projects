disp('London Baby!!');
a=0.1;
truevalue=exp(a);
for i=1:5
    h=10^(-i);
    aprvalue(i)=(1+h)^(a/h);
    err(i)=abs(truevalue-aprvalue(i));
    hall(i)=h;
end

loglog(hall,err,':g*')

