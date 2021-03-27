 disp('London Babyy');
 a=1;
 h=0.01;
 truevalue=1/(1+a^2);
 aprvalue=(atan(a+h)-atan(a))/h;
 err=abs(truevalue-aprvalue);
 
 