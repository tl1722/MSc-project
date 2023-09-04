%% DSB beamforming in time domain (with RIR)

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master';

[source_signal,fs]=readwav('S_01_01.wav');
L=length(source_signal);


%% parameters
c = 340;                  % Sound velocity (m/s)
num_mics=4;               % micphone number
% distance between microphone (linear array)
d=0.035;
% Microphone array position [x y z] (m)
position_mics = [2 1.5 2; 
                 2+d 1.5 2;
                 2+2*d 1.5 2;
                 2+3*d 1.5 2];   
% Source position [x y z] (m)
position_source = [3 3.5 2];    
% DOA
theta_source=atan2(position_source(1,2)-position_mics(1,2), position_source(1,1)-position_mics(1,1));
theta_angle=theta_source*180/pi;


%% RIR generator
room = [5 4 6];                % Room dimensions [x y z] (m)
rev_time = 0.6;             % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir = rir_generator(c, fs, position_mics, position_source, room, rev_time, n, mtype, order, dim, orientation, hp_filter);

figure;
plot(rir(1,:));
title('RIR')
%% combine the rir to the speech signal

% get the received signal x in each microphone
y=zeros(num_mics,length(source_signal));
for i=1:num_mics
  y(i,:)=filter(rir(i,:),1,source_signal).*10;
end

% plot the received signal in each microphone
figure;
subplot(5,1,1)
plot(source_signal);
title('source signal')
subplot(5,1,2)
plot(y(1,:));
title('received signal in Mic 1')
subplot(5,1,3)
plot(y(2,:));
title('received signal in Mic 2')
subplot(5,1,4)
plot(y(3,:));
title('received signal in Mic 3')
subplot(5,1,5)
plot(y(4,:));
title('received signal in Mic 4')
% sound(y(1,:),fs)


%% delay and sum
% compute delay for the received signal in each microphone
c=340;
distances=[0:num_mics-1]*d*cos(theta_source);
tao=distances/c;    
delay_samples=round(tao * fs);

% conpensate delay to each microphone
max_delay_samples = max(delay_samples);
y_beam=zeros(num_mics, L+max_delay_samples);
for i=1:num_mics
    delay_mic=delay_samples(i);
    y_beam(i,delay_mic+1:delay_mic+L) = (1/num_mics)*y(i,:);
end

% sum the y in each microphone
y_sum=sum(y_beam,1);
figure;
plot(source_signal);
hold on
plot(y_sum)
title('DSB in time domain')
legend('input signal','output signal')
sound(y_sum,fs)







