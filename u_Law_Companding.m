clear all;
u=7;
Fs=100e3;
t=[0:1/Fs:255*1/Fs]; %63,127,255 for different number of samples
L=length(t);
x=1.5*cos(2*pi*t*5e3)+sin(2*pi*t*4e3);
x=x/max(x);
B=4;
p=2^B;
MSE=0;
MSE_q=0;

if x*2^(B-1)-floor(x*2^(B-1))<0.5
    x_q=floor(x*2^(B-1))/2^(B-1);
else
    x_q=ceil(x*2^(B-1))/2^(B-1);
end

MSE=sum((x_q-x).^2);

y=sign(x).*log(1+u.*abs(x))/log(1+u);

if y*2^(B-1)-floor(y*2^(B-1))<0.5
    y_q=floor(y*2^(B-1))/2^(B-1);
else
    y_q=ceil(y*2^(B-1))/2^(B-1);
end

z=sign(y_q).*((1+u).^abs(y_q)-1)/u;

MSE_q=MSE_q+sum((z-x).^2);

figure(1);
subplot(3,1,1);
plot(t,x);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,2);
plot(t,y);
title('Compressed Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,3);
plot(t,z);
title('Unequally quantized Signal');
xlabel('time(secs)');
ylabel('Amplitute');


figure(2);
subplot(3,1,1);
plot(t,x);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,2);
plot(t,x_q);
title('Equally Quantized Signal');
xlabel('time(secs)');
ylabel('Amplitute');