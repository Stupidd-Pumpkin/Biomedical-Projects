clc
clear all

% instrument connection
%% Tektronix MDO4104B
% Find a VISA-USB object.
objOS = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0450::C011393::0::INSTR', 'Tag', '');

% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(objOS)
    objOS = visa('Tek', 'USB0::0x0699::0x0450::C011393::0::INSTR');
else
    fclose(objOS);
    objOS = objOS(1);
end
% Connect to instrument object, obj1.
fopen(objOS);
fclose(objOS);
objOS.InputBufferSize = 100000;
fopen(objOS);
%set up oscillascope
%acqusion mode
fprintf(objOS, 'ACQuire:STOPAfter RUNSTop');%set in runstop mode
fprintf(objOS, 'ACQuire:STATE RUN');%set it run

%% source meter
visa_SMU2450 = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x05E6::0x2450::04330068::0::INSTR', 'Tag','');

	% Create the VISA-USB object if it does not exist
	% otherwise use the object that was found.
	if isempty(visa_SMU2450)
		visa_SMU2450 = visa('Tek', 'USB0::0x05E6::0x2450::04330068::0::INSTR');
	else
		fclose(visa_SMU2450);
		visa_SMU2450 = visa_SMU2450(1);
	end

fopen(visa_SMU2450);
    
%% setup arduino    
a = arduino('COM8','due','Libraries','SPI');
dev = spidev(a,'D10','Mode',0,'BitOrder','lsbfirst','Bitrate',10000);
% %small gain
% data1 = 47;
% out1  = writeRead(dev,data1)
% data2 = 113;
% out2  = writeRead(dev,data2)
%big gain
data1 = 45;
out1  = writeRead(dev,data1);
data2 = 113;
out2  = writeRead(dev,data2);
figure()

for i = 1:31
 
fprintf(visa_SMU2450, ':OUTP:STAT ON');        % activate output source   
% inputC = [':SOUR:CURRENT ' num2str((i-1)*0.5e-6)];
inputC = [':SOUR:VOTL ' num2str((i-1)*0.5e-6)];
fprintf(visa_SMU2450, inputC);
pause(2)
fprintf(objOS, 'ACQuire:STATE STOP');%set scope stop
pause(2)
% % read channel data
%% channel 1
clear data1 p1 dat1
fprintf(objOS, ':DATa:SOUrce CH1');
fprintf(objOS, ':DATa:STARt 1');
fprintf(objOS, ':DATa:STOP 10000');
fprintf(objOS, ':DATa:ENCdg ASCIi');
fprintf(objOS, ':DATa:WIDth 1');
fprintf(objOS, ':HEADer 1');
fprintf(objOS, ':VERBose');
mesg1{i} = query(objOS, ':WFMOutpre?');
fprintf(objOS, ':HEADer 0');
data1 = query(objOS, ':CURVe?');
p1 = regexp(data1,',','split');
m = size(p1,2);
for j=1:m
    dat1(j) = str2double(p1{j});
end

ch1(:,i) = dat1+128;

thr(i) = mean(ch1(:,i))/128;
plot(i, thr(i),'o');hold on;
% if thr(i) < 0.4
% data1 = 47;
% out1  = writeRead(dev,data1);
% data2 = 113;
% out2  = writeRead(dev,data2);
% pause(2)
% end
pause(2)
fprintf(objOS, 'ACQuire:STOPAfter RUNSTop');%set in runstop mode
fprintf(objOS, 'ACQuire:STATE RUN');%set it run
info = ['End of #' num2str(i) '...\n'];
fprintf(info);
i = i+1;
end

save('DataAGC3', 'ch1','thr');

load('DataAGC1.mat');
i=(1:30)/2;
figure()
plot(-thr(2:end)+thr(1),'ro');hold on
load('DataAGC2.mat');
plot(-thr(2:end)+thr(1),'bo');hold on
load('DataAGC3.mat');
plot(-thr(2:24)+thr(1),'ko');hold on
