clear all;
Fs=8000;
record=audiorecorder(Fs,24,1); 
disp('recording started');
recordblocking(record,5); %%record for 5 seconds and stop recording
disp('recording stopped');
pause(2);
disp('playing record');
play(record);
audio=getaudiodata(record);
t=[1/Fs:1/Fs:length(audio)/Fs];

plot(t,audio);
disp('done :)');

