clear all;
A=87.6;
Fs=100e3;
t=[0:1/Fs:255*1/Fs];
L=length(t);
x=1.5*cos(2*pi*t*5e3)+sin(2*pi*t*4e3);
B=6;
p=2^B;
x=x/max(x);
MSE=0;  
MSE_q=0;

x_q=round(x*2^(B-1))/2^(B-1);

MSE=sum((x_q-x).^2);

for i=1:1:length(x)
    if abs(x(i))<1/A
        y(i)=sign(x(i))*A*abs(x(i))/(1+log(A));
    else 
        y(i)=sign(x(i))*(1+log(A*abs(x(i))))/(1+log(A));
    end
end

y_q=round(y*2^(B-1))/2^(B-1);

for i=1:1:length(x)
    if abs(y_q(i))<1/(1+log(A))
        z(i)=sign(y_q(i))*abs(y_q(i))*(1+log(A))/A;
    else
        z(i)=sign(y_q(i))*exp(abs(y_q(i))*(1+log(A))-1)/A;
    end

end

MSE_q=MSE_q+sum((z-x).^2);
figure(1);
subplot(3,1,1);
plot(t,x);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,2);
plot(t,y);
title('Compressed & Quantized Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,3);
plot(t,z);
title('Unequally quantized Signal');
xlabel('time(secs)');
ylabel('Amplitute');


% figure(2);
% subplot(3,1,1);
% plot(t,x);
% title('Input Signal');
% xlabel('time(secs)');
% ylabel('Amplitute');
% 
% subplot(3,1,2);
% plot(t,x_q);
% title('Equally Quantized Signal');
% xlabel('time(secs)');
% ylabel('Amplitute');