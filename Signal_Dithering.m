clear all;

Fs=100e3;
t=(0:1/Fs:511*1/Fs); %63,127,255 for different number of samples
L=length(t);

x=10*cos(2*pi*t*1e3)+6*cos(2*pi*t*2e3)+2*cos(2*pi*t*4e3);

xm=x/max(abs(x));
B=12;
B_Final=8;

if xm*2^(B_Final-1)-floor(xm*2^(B_Final-1))<0.5
    x_q=floor(xm*2^(B_Final-1))/2^(B_Final-1);
else
    x_q=ceil(xm*2^(B-1))/2^(B-1);
end

if xm*2^(B-1)-floor(xm*2^(B-1))<0.5
    y_x=floor(xm*2^(B-1))/2^(B-1);
else
    y_x=ceil(xm*2^(B-1))/2^(B-1);
end
% y=y_x;
% y=awgn(y_x,20*(B_Final-1)*0.301);
n=2^(B-B_Final)*rand(size(y_x))-2^(B-B_Final-1);

if n*2^(B-B_Final)-floor(n*2^(B-B_Final))<0.5
    n=floor(n*2^(B-B_Final))/2^(B-B_Final);
else
    n=ceil(n*2^(B-B_Final))/2^(B-B_Final);
end
n=(2^(B-B_Final)*n)/2^(B);
y=y_x+n;

if y*2^(B_Final-1)-floor(y*2^(B_Final-1))<0.5
    z=floor(y*2^(B_Final-1))/2^(B_Final-1);
else
    z=ceil(y*2^(B_Final-1))/2^(B_Final-1);
end
audio=audioplayer(z,Fs);
play(audio);

figure(1);
subplot(3,1,1);
plot(t,n);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,2);
plot(t,y);
title('Noise Added Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,3);
plot(t,z);
title('Quantized Signal');
xlabel('time(secs)');
ylabel('Amplitute');

figure(2);
subplot(3,1,1);
plot(t,y_x);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');

subplot(3,1,2);
plot(t,x_q);
title('Direct Quantized Signal');
xlabel('time(secs)');
ylabel('Amplitute');


% n=wgn(length(x),1,0);
% nm=n/max(abs(n));
% B_noise=B-3;
% if nm*2^(B_noise-1)-floor(nm*2^(B_noise-1))<0.5
%     y_noise=floor(nm*2^(B_noise-1))/2^(B_noise-1);
% else
%     y_noise=ceil(nm*2^(B_noise-1))/2^(B_noise-1);
% end
% 
% for i=1:1:length(x)
%     y(i)=y_x(i)+y_noise(i);
% end