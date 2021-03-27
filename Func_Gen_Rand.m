%% MATLAB script to generate custom waveform from Tektronix AFG
% This example illustrates the use of MATLAB to set up a Tektronix
% AFG3000 series to generate a custom waveform created in MATLAB. Tektronix
% VISA will be used for communicating with your instrument, please ensure
% that it is installed. 
% Xinyao Tang
% 05/01/2017

%% Instrument Connection (Dual Function Generator)
% Clear MATLAB workspace of any previous instrument connections
instrreset;

% Create a VISA-USB object.
objDF = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0352::C010624::0::INSTR', 'Tag', '');

% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(objDF)
    objDF = visa('Tek', 'USB0::0x0699::0x0352::C010624::0::INSTR');
else
    fclose(objDF);
    objDF = objDF(1);
end


% Change the |OutputBufferSize| depending on the size of the custom
% waveform being transfered.
objDF.OutputBufferSize = 10000; 

% Set the |ByteOrder| to match the requirement of the instrument
objDF.ByteOrder = 'littleEndian';

% Connect to instrument object, obj1.
fopen(objDF);
% Reset the function generator to a know state
fprintf(objDF, '*RST');
fprintf(objDF, '*CLS;'); 

%Sine wave
fsine=250.5e3;
fs = 6*fsine;
t=1/fs:1/fs:1000/fs;
L=length(t);
i=1;
while(i<=L)
    x(i)=sin(2*pi*fsine*t(i));
    i=i+1;
end

D = randi([0 1],100,1);
j=0;
k=1;
for i=1:L
    j=j+1;
    waveform(i)=x(i)*D(k);
    if j>=L/100
        k=k+1;
        j=0;
    end
end


% Convert the double values integer values between 0 and 16382 (as required
% by the instrument)
waveform =  round((waveform + 1.0)*8191);
waveformLength = length(waveform);

% Encode variable 'waveform' into binary waveform data for AFG.  This is the
% same as AWG5000B but marker bits are ignored. Refer to AWG5000B series
% programmer manual for bit definitions.
binblock = zeros(2 * waveformLength, 1);
binblock(2:2:end) = bitand(waveform, 255);
binblock(1:2:end) = bitshift(waveform, -8);
binblock = binblock';

% Build binary block header
bytes = num2str(length(binblock));
header = ['#' num2str(length(bytes)) bytes];

% Resets the contents of edit memory and define the length of signal
fprintf(objDF, ['DATA:DEF EMEM, ' num2str(length(t)) ';']); %1001

% Transfer the custom waveform from MATLAB to edit memory of instrument
fwrite(objDF, [':TRACE EMEM, ' header binblock ';'], 'uint8');

% Associate the waveform in edit memory to channel 1
fprintf(objDF, 'SOUR1:FUNC EMEM');

% Set the output frequency to 1 Hz. This is important since the custom
% waveform's frequency is further upconverted to the value set below.
fprintf(objDF, 'SOUR1:FREQ:FIXed 1Hz');

% Turn on Channel 1 output
fprintf(objDF, ':OUTP1 ON');
