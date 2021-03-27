clear all;
signal=load('C:\Users\Mohammed Sameer\Documents\MATLAB\ECG_400.csv');
x=signal(:,1);
fs=400;
Time=(0:1/fs:length(x)/fs-1/fs)';
f=1./Time;
X=fft(x);
hold on;
for i=1:length(f)
    if f(i)==50
        X(i)=X(i-1);
    end
    if f(i)>100
        X(i)=0;
    end
end
Voltage=real(ifft(X));
plot(Time,Voltage,'g');
xlabel('Time');ylabel('Volts');grid;title('ECG Signal');hold on;
beats=0;
%Extraction of R,S,T
R=[];RTime=[];
S=[];STime=[];
T=[];TTime=[];
i=11;
while(i<length(Voltage)-10)
    i=i+1;
    if Voltage(i)==max(Voltage(i-10:1:i+10)) && Voltage(i)>0.6
        R=[R,Voltage(i)];
        RTime=[RTime,Time(i)];
        i=i+10;
    end
    if Voltage(i)==min(Voltage(i-10:1:i+10)) && Voltage(i)<-0.4
        S=[S,Voltage(i)];
        STime=[STime,Time(i)];
        i=i+10;
    end
    if Voltage(i)==max(Voltage(i-10:1:i+10)) && Voltage(i)<0.6 && Voltage(i)>0.2
        T=[T,Voltage(i)];
        TTime=[TTime,Time(i)];
        i=i+10;
    end
end
plot(RTime,R,'bo');plot(STime,S,'ko');plot(TTime,T,'ro');

%Displaying the Info
RRTimeInterval=(RTime(end)-RTime(1))/(length(R)-1);
disp('R-R Time interval is:');
disp(RRTimeInterval);

STsample=[];
for i=1:min(length(S),length(T))
    if TTime(1)>STime(1)
        STsample=[STsample,TTime(i)-STime(i)];
    else
        STsample=[STsample,TTime(i+1)-STime(i)];
    end
end
STTimeInterval=mean(STsample);
disp('S-T Time interval is:');
disp(STTimeInterval);
    
HeartBeatRate=60/RRTimeInterval;
disp('Heart Beat Rate is:');
disp(HeartBeatRate);
hold off;