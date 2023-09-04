%% ex4: MVDR beamforming (use library)
clear
clc

[source_signal1,fs]=readwav('S_01_01.wav');
L1=length(source_signal1);

%% microphone parameters
num_mics=8;
c = 340;  % speed of sound, in m/s
f=1000;
lambda=c/f;
microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 40000]);
ula = phased.ULA('NumElements',num_mics,'ElementSpacing', lambda/2,'Element',microphone);

angle_s1 = [20; 0];   % direction of target signal

%% generate target signal
collector = phased.Collector('Sensor',ula,'PropagationSpeed',c, 'OperatingFrequency',fs);
x1=collector(source_signal1, angle_s1);

figure;
subplot(311); 
plot(real(x1));axis tight;
title('Source signal');xlabel('Time (s)');ylabel('Magnitude (V)');


%% generate noise
% random seed
rs = RandStream.create('mt19937ar','Seed',200);

noisePwr = 0.001;   % noise power, 30dB SNR 
noise= sqrt(noisePwr/2)*randn(L1, num_mics);
rxSignal = x1 + noise; 

subplot(312); 
plot(real(rxSignal(:,1)));axis tight;
title('Received signal at Mic 1');xlabel('Time (s)');ylabel('Magnitude (V)');

% MVDR beamformer
mvdrbeamformer = phased.MVDRBeamformer('SensorArray',ula,...
    'Direction',angle_s1,'PropagationSpeed',c,'OperatingFrequency',fs,...
    'WeightsOutputPort',true);

%mvdrbeamformer.TrainingInputPort = true;
[yMVDR,wMVDR] = mvdrbeamformer(rxSignal);

%% plot 
% output signal
subplot(313)
plot(real(yMVDR)); axis tight;
title('Output of MVDR Beamformer');
xlabel('Time (s)');ylabel('Magnitude (V)');
% sound(real(yMVDR),fs)

% beam pattern of MVDR
figure;
pattern(ula,fs,-180:180,0,'Weights',wMVDR,'Type','powerdb',...
    'PropagationSpeed',c,'Normalize',false,...
    'CoordinateSystem','rectangular');


% compare with DSB
psbeamformer = phased.PhaseShiftBeamformer('SensorArray',ula,...
    'OperatingFrequency',fs,'PropagationSpeed',c,'Direction',angle_s1,...
    'WeightsOutputPort', true);
[yCbf,w] = psbeamformer(rxSignal);
hold on
pattern(ula,fs,-180:180,0,'Weights',w,'Type','powerdb',...
    'PropagationSpeed',c,'Normalize',false,...
    'CoordinateSystem','rectangular');
hold off
legend('MVDR','DSB')
axis([-180 180 -50 20]);

%% consider interference signal
angle_i1 = [30; 0];   % direction of interference signal 1
angle_i2 = [60; 0];   % direction of interference signal 2

i1 = randn(rs,L1,1);
i2 = randn(rs,L1,1);
% [source_signal2,fs2]=readwav('S_01_02.wav');
% l1=length(source_signal1);
% l2=length(source_signal2);
% i1=[source_signal2;zeros(l1-l2,1)];
% 
% [source_signal3,fs3]=readwav('S_01_03.wav');
% l3=length(source_signal3);
% i2=[source_signal3;zeros(l1-l3,1)];

% interference at 30 degrees and 60 degrees
interference = collector([i1,i2],[angle_i1,angle_i2]);
rxInt = interference + noise;                 % total interference + noise
rxSignal = x1 + rxInt;                % total received Signal


figure;
subplot(311); 
plot(real(x1));axis tight;
title('Source signal');xlabel('Time (s)');ylabel('Magnitude (V)');

subplot(312); 
plot(real(rxSignal(:,1)));axis tight;
title('Received signal at Mic 1');xlabel('Time (s)');ylabel('Magnitude (V)');

%% MVDR beamformer
mvdrbeamformer = phased.MVDRBeamformer('SensorArray',ula,...
    'Direction',angle_s1,'PropagationSpeed',c,'OperatingFrequency',fs,...
    'WeightsOutputPort',true);

%mvdrbeamformer.TrainingInputPort = true;
[yMVDR,wMVDR] = mvdrbeamformer(rxSignal);

%% plot 
% output signal
subplot(313)
plot(real(yMVDR)); axis tight;
title('Output of MVDR Beamformer With Presence of Interference');
xlabel('Time (s)');ylabel('Magnitude (V)');
sound(real(yMVDR),fs)

% beam pattern of MVDR
figure;
pattern(ula,fs,-180:180,0,'Weights',wMVDR,'Type','powerdb',...
    'PropagationSpeed',c,'Normalize',false,...
    'CoordinateSystem','rectangular');
title('Beam pattern of MVDR Beamformer with interference');
axis([-180 180 -100 20]);


% LCMV
steeringvector = phased.SteeringVector('SensorArray',ula,...
    'PropagationSpeed',c);
LCMVbeamformer = phased.LCMVBeamformer('DesiredResponse',1,...
    'TrainingInputPort',true,'WeightsOutputPort',true);
LCMVbeamformer.Constraint = steeringvector(fs,angle_s1);
LCMVbeamformer.DesiredResponse = 1;
[yLCMV,wLCMV] = LCMVbeamformer(rxSignal,rxInt);

% beam pattern of LCMV
figure;
pattern(ula,fs,-180:180,0,'Weights',wLCMV,'Type','powerdb',...
    'PropagationSpeed',c,'Normalize',false,...
    'CoordinateSystem','rectangular');
title('Beam pattern of LCMV Beamformer with interference');
hold on
pattern(ula,fs,-180:180,0,'Weights',wMVDR,'Type','powerdb',...
    'PropagationSpeed',c,'Normalize',false,...
    'CoordinateSystem','rectangular');
hold off
legend('LCMV','MVDR')
axis([-180 180 -100 20]);
