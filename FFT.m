Signal=load('signal.csv');
Time=Signal(:,1);                        %the data has 1st coloum as time samples
x=Signal(:,2);                           %the data has 2nd coloum as signal
fs=1/(Time(2)-Time(1));                  %sample rate


% % %some arbitriary example
% % fs=100000;                          %defining sampling rate
% % Time=(0:1/fs:5);                    %setting samples
% % w1=2*pi*40000;                    
% % w2=2*pi*7500;            
% % x=3*sin(w1*Time)+ sin(w2*Time);
% % %plot(Time,x,'r');

L=length(Time);                     %Number of samples

X=fft(x);               %applying Fast Fourier Transform
X=X(1:1:L/2);                       %for producing the single sided spectrum
f=0:fs/L:fs/2-fs/L;                 %defining bins 
plot(f,2*abs(X)/L);


for i=2:1:(L/2)-1  
    if abs(X(i))>abs(X(i-1)) && abs(X(i))>abs(X(i+1)) && abs(X(i))>mean(abs(X))
      fprintf('the given signal has a sinusoid with frequency: %d and amplitude: %d \n',f(i),2*abs(X(i)/L));
    end
end