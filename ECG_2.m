clear all;
signal=load('C:\Users\Mohammed Sameer\Documents\MATLAB\ECG\n06m');

x_1=signal.val(1,:);
x_1=x_1(round(length(x_1)/10):round(9*length(x_1)/10));
x_2=signal.val(2,:);
x_2=x_2(round(length(x_2)/10):round(9*length(x_2)/10));
x_1=x_1-mean(x_1);
x_1=x_1/max(abs(x_1));
x_2=x_2-mean(x_2);
x_2=x_2/max(abs(x_2));

Fs=128;     %128 for normal signals and 250 for abnormal signals
t=(0:1/Fs:(length(x_1)-1)/Fs);
L=length(t);
F=-Fs/2:Fs/L:Fs/2-Fs/L;

[B,A]=butter(3,0.5/128,'high');
x_1=filter(B,A,x_1);
x_2=filter(B,A,x_2);
[B,A]=butter(3,[49/128,51/128],'stop','s');
x_1=filter(B,A,x_1);
x_2=filter(B,A,x_2);
[B,A]=butter(1,40/128);
x_1=filter(B,A,x_1);
x_2=filter(B,A,x_2);

X_1=fftshift(fft(x_1));
X_2=fftshift(fft(x_2));

En1(1:length(x_1))=0;
En2(1:length(x_2))=0;

for i=2:1:length(x_1)
    if F(i)>=0
        En1(i)=En1(i-1)+abs(X_1(i))^2;
        En2(i)=En2(i-1)+abs(X_2(i))^2;
    end
end

for i=2:1:length(x_1)
    if En1(i)>=max(En1)/2
        Fc1=F(i);
        break;
    end
end
for i=2:1:length(x_1)
    if En2(i)>=max(En2)/2
        Fc2=F(i);
        break;
    end
end

P=0;
Q=0;
k=1;
i=1;
Fd1(1)=0;
while(i<length(x_1)*999/1000)
    if F(i)>=0;
        Fd1(k+1)=0;
        for j=i:i+length(x_1)/1000
            P=P+abs(X_1(j))^2;
            q=Fd1(k);
            Fd1(k)=max(Fd1(k),abs(X_1(j))^2);
            if q~=Fd1(k)
                fbin1(k)=F(j);
            end
        end
        if P>max(En1)/20;
            Ebin1(k)=P;
            k=k+1;
            i=j;
        end
        P=0;
    end
    i=i+1;
end
fbin1(length(fbin1))=[];

P=0;
Q=0;
k=1;
i=1;
Fd2(1)=0;
while(i<length(x_2)*999/1000)
    if F(i)>=0;
        Fd2(k+1)=0;
        for j=i:i+length(x_2)/1000
            P=P+abs(X_2(j))^2;
            q=Fd2(k);
            Fd2(k)=max(Fd2(k),abs(X_2(j))^2);
            if q~=Fd2(k)
                fbin2(k)=F(j);
            end
        end
        if P>max(En2)/20;
            Ebin2(k)=P;
            k=k+1;
            i=j;
        end
        P=0;
    end
    i=i+1;
end
fbin2(length(fbin2))=[];

figure(1);
subplot(2,1,1);
plot(t,x_1);
title('Time Domain Signal 1');
xlabel('Time(secs)');
ylabel('Normalized Amplitute');

subplot(2,1,2);
plot(t,x_2);
title('Time Domain Signal 2');
xlabel('Time(secs)');
ylabel('Normalized Amplitute');

figure(2);
subplot(2,1,1);
plot(F,2*abs(X_1).^2/L);
title('Magnitute Spectrum of signal 1');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

subplot(2,1,2);
plot(F,2*abs(X_2).^2/L);
title('Magnitute Spectrum of signal 2');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

Threshold1=length(fbin1)
Threshold2=Fc1
Threshold3=length(fbin2)
Threshold4=Fc2