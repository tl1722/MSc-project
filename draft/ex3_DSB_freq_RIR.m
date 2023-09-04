%% ex2: DSB  beamforming in frequency domain (with RIR)
clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% load the audio signal
[source_signal,fs]=readwav('S_01_01.wav');
L=length(source_signal);
f=3000;

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

%% combine the rir to the speech signal

% get the received signal x in each microphone
received_signal=zeros(num_mics,length(source_signal));
for i=1:num_mics
  received_signal(i,:)=filter(rir(i,:),1,source_signal).*10;
end

% plot the received signal in each microphone
figure;
subplot(5,1,1)
plot(source_signal);
title('source signal')
subplot(5,1,2)
plot(received_signal(1,:));
title('received signal in Mic 1')
subplot(5,1,3)
plot(received_signal(2,:));
title('received signal in Mic 2')
subplot(5,1,4)
plot(received_signal(3,:));
title('received signal in Mic 3')
subplot(5,1,5)
plot(received_signal(4,:));
title('received signal in Mic 4')
% sound(y(1,:),fs)

%% frequency domain frame-based processing for each microphone
overlap_factor = 4;   
% frame increment in samples
inc = 40;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = hamming(N_window,'periodic');
nfft=N_window;

for i=1:num_mics
  [frame,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa');
  n_frame=size(frame,1);

  % compute the delay in ith sensor
  distance=(i-1)*d*cos(theta_source);
  tao=distance/c; 
  % DSB weight
  w=[1,1,1,1]/num_mics;

  % for each frame, conpensate the delay by the factor exp(-j*w*tao)
  for j=1:n_frame
     beam_frame(j,:,i)=w(i)*(frame(j,:)*exp(-i*2*pi*f*tao));
  end
end

% sum Y in every sensor
all_frame=sum(beam_frame,3);
output_freq=reshape(all_frame,size(frame,1),size(frame,2));

% use overlapadd to reconstruct the signal y
output_signal=v_overlapadd(v_irfft(output_freq,N_window,2),window,inc); 


%% output the signal
writewav(output_signal,fs,'output_signal_01.wav')
figure;
plot(received_signal(1,:));
hold on
plot(output_signal)
title('DSB in frequency domain')
legend('input signal','output signal')
sound(output_signal,fs)

