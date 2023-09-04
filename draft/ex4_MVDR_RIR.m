%% ex5: LCMV in time domain (with RIR)
clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% load the audio signal
[source,fs]=readwav('S_01_01.wav');
L=length(source);
vn1 = 3*randn(L,1);  % interference signal 1
vn2 = 3*randn(L,1);  % interference signal 2

%% parameters
c = 340;                  % Sound velocity (m/s)
num_mics=4;               % micphone number
% distance between microphone (linear array)
d=0.2;
% Microphone array position [x y z] (m)
position_mics = [1 1.5 2; 
                 1+d 1.5 2;
                 1+2*d 1.5 2;
                 1+3*d 1.5 2];   
% Source position [x y z] (m)
position_source = [2 3.5 2]; 

% Source position [x y z] (m)
position_interference1 = [3 3 2];  
position_interference2 = [3 4 2];  

% DOA
theta_source=atan2(position_source(1,2)-position_mics(1,2), position_source(1,1)-position_mics(1,1));
theta_angle=theta_source*180/pi;

% direction of intereference
theta_interference1=atan2(position_interference1(1,2)-position_mics(1,2), position_interference1(1,1)-position_mics(1,1));
theta_angle_i1=theta_interference1*180/pi;

% direction of intereference
theta_interference2=atan2(position_interference2(1,2)-position_mics(1,2), position_interference2(1,1)-position_mics(1,1));
theta_angle_i2=theta_interference2*180/pi;

%% RIR generator
room = [5 4 6];                % Room dimensions [x y z] (m)
rev_time = 0;             % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=theta_source;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, n, mtype, order, dim, orientation, hp_filter);
rir_undesired1 = rir_generator(c, fs, position_mics, position_interference1, room, rev_time, n, mtype, order, dim, orientation, hp_filter);
rir_undesired2 = rir_generator(c, fs, position_mics, position_interference2, room, rev_time, n, mtype, order, dim, orientation, hp_filter);

%% combine the rir to the speech signal

% get the received signal x in each microphone
desired_signal=zeros(num_mics,length(source));
interference_signal_1=zeros(num_mics,length(source));
interference_signal_2=zeros(num_mics,length(source));

for i=1:num_mics
  desired_signal(i,:)=filter(rir_desired(i,:),1,source).*10;
  interference_signal_1(i,:)=filter(rir_undesired1(i,:),1,vn1).*10;
  interference_signal_2(i,:)=filter(rir_undesired2(i,:),1,vn2).*10;
end

received_signal=desired_signal+interference_signal_1+interference_signal_2;

% plot the received signal in each microphone
figure;
subplot(5,1,1)
plot(source);
title('source signal')
subplot(5,1,2)
plot(desired_signal(1,:));
title('received signal in Mic 1')
subplot(5,1,3)
plot(interference_signal_1(1,:));
title('interference signal in Mic 1')
% sound(y(1,:),fs)

% LCMV
R_x = 1/L*(received_signal*received_signal');
R_x_inv = inv(R_x);

n_fft=160;
f=linspace(0, fs/2, n_fft);
w_lcmv=zeros(num_mics,n_fft);
w_mvdr=zeros(num_mics,n_fft);
for i=1:length(f)
  vs = exp(-1i*2*pi*[0:num_mics-1].'*d*cos(theta_source)*f(i)/c)/num_mics; % signal direction vector
  vn1 = exp(-1i*2*pi*[0:num_mics-1].'*d*cos(theta_interference1)*f(i)/c)/num_mics; % interference direction vector
  vn2 = exp(-1i*2*pi*[0:num_mics-1].'*d*cos(theta_interference2)*f(i)/c)/num_mics; % interference direction vector
  C=[vs,vn1,vn2];
  F=[1,0,0].';
  w_lcmv(:,i) = (R_x_inv*C) / (C'*R_x_inv*C)*F;
  w_mvdr(:,i) = (R_x_inv*vs) / (vs'*R_x_inv*vs);
end

% plot the beam pattern
theta = 180 * [-1:0.01:1];
B=zeros(n_fft,length(theta));
figure;
for i=1:length(f)
  v = exp(-1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
  B(i,:) = abs(w_mvdr(:,i)' * v);
end
k=100;
plot(theta,20*log10(B(k,:)/max(B(k,:))))
hold on
title('Beam pattern of LCMV')
xlabel('degree')
ylabel('dB')
grid on

figure;
P=20*log10(B);
mesh(f,theta,P.');
