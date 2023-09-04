%% ex5: LCMV and MVDR in frequency domain (with RIR)
clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% load the audio signal
[source,fs]=readwav('S_01_01.wav');
L=length(source);
vn1 = 0.01*randn(L,1);  % interference signal 1
vn2 = randn(L,1);  % interference signal 2

%% parameters
c = 340;                  % Sound velocity (m/s)
num_mics=4;               % micphone number
% distance between microphone (linear array)
d=0.035;
% Microphone array position [x y z] (m)
position_mics = [1 1.5 2; 
                 1+d 1.5 2;
                 1+2*d 1.5 2;
                 1+3*d 1.5 2];   
% Source position [x y z] (m)
position_source = [2 3.5 2]; 
% interference position [x y z] (m)
position_interference = [3 3 2];  
% DOA
theta_source=atan2(position_source(1,2)-position_mics(1,2), position_source(1,1)-position_mics(1,1));
angle_source=theta_source*180/pi;

% direction of intereference
theta_interference=atan2(position_interference(1,2)-position_mics(1,2), position_interference(1,1)-position_mics(1,1));
angle_interference=theta_interference*180/pi;

%% RIR generator
room = [8 6 5];                % Room dimensions [x y z] (m)
rev_time = 0;             % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, n, mtype, order, dim, orientation, hp_filter);
rir_undesired = rir_generator(c, fs, position_mics, position_interference, room, rev_time, n, mtype, order, dim, orientation, hp_filter);

%% combine the rir to the speech signal

% get the received signal x in each microphone
desired_signal=zeros(num_mics,length(source));
interference_signal_1=zeros(num_mics,length(source));

for i=1:num_mics
  desired_signal(i,:)=filter(rir_desired(i,:),1,source)*10;
  interference_signal_1(i,:)=filter(rir_undesired(i,:),1,vn1)*10;
end

received_signal = desired_signal + interference_signal_1;


%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc = 512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = hamming(N_window,'periodic');

% segment signal into frames and do DFT
for i=1:num_mics
  [frame,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa'); 
  n_frame=size(frame,1);
  n_fft=size(frame,2);
  total_frame(i,:,:)=frame;
end

% set steering vector
freqVec=linspace(0, fs/2, n_fft);
vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source)/c); % signal direction vector
vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interference)/c); % interference direction vector

w_mvdr = zeros(num_mics, length(freqVec));
w_lcmv = zeros(num_mics, length(freqVec));
eps=1e-4;
for f = 1:length(freqVec)  
   Phi=1/n_frame * (total_frame(:,:,f)*total_frame(:,:,f)');
   
   % LCMV
   C = [vs(:, f),vn(:,f)];
   denominator = C'/Phi*C;
   if rcond(denominator)<eps
     denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
   end
   F=[1,0].';
   w_lcmv(:, f) = (Phi\C) / (denominator)*F;

   % MVDR
   vf=vs(:, f);
   w_mvdr(:, f) = (Phi\vf) / (vf'/Phi*vf);

end

for i=1:num_mics
  [frame,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa'); 
  for n=1:n_frame
      %beam_frame(n,:,i)=w_mvdr(i,:).*(frame(n,:));
      beam_frame(n,:,i)=w_lcmv(i,:).*(frame(n,:));
  end
end

% sum Y in every sensor
sum_frame=sum(beam_frame,3);

% use overlapadd to reconstruct the signal y
output_signal=v_overlapadd(v_irfft(sum_frame,N_window,2),window,inc); 


%% plot
% % plot waveform
% figure;
% subplot(4,1,1)
% plot(source);
% title('source signal')
% subplot(4,1,2)
% plot(desired_signal(1,:));
% title('desired signal in Mic 1')
% subplot(4,1,3)
% plot(received_signal(1,:));
% title('received signal in Mic 1')
% % sound(desired_signal(1,:),fs)

% writewav(output_signal,fs,'output_signal_01.wav')
% subplot(4,1,4);
% plot(received_signal(1,:));
% hold on
% plot(output_signal)
% title('Beamformed signal (MVDR)')
% legend('received signal in mics','beamformed signal')
% sound(output_signal,fs)

% spectrogram
figure;
subplot(1,2,1)
[s,f,t] = spectrogram(received_signal(1,:),window,inc,N_window,fs) ;
waterfall(f,t,abs(s)'.^2)
set(gca,XDir="reverse",View=[30 50])
subplot(1,2,2)
[s,f,t] = spectrogram(output_signal,window,inc,N_window,fs) ;
waterfall(f,t,abs(s)'.^2)
set(gca,XDir="reverse",View=[30 50])
xlabel("Frequency (Hz)")
ylabel("Time (s)")


% plot the beam pattern
theta = 180 * [-1:0.01:1];
B_mvdr=zeros(n_fft,length(theta));
B_lcmv=zeros(n_fft,length(theta));
for i=1:length(f)
  v = exp(-1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
  B_mvdr(i,:) = abs(w_mvdr(:,i)' * v);
  B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);

end
% figure;
% k=100;
% plot(theta,20*log10(B_mvdr(k,:)/max(B_mvdr(k,:))))
% hold on
% title('Beam pattern of LCMV')
% xlabel('degree')
% ylabel('dB')
% grid on

figure;
subplot(1,2,1)
P_mvdr=20*log10(B_mvdr);
mesh(f,theta,P_mvdr.');
set(gca,YDir="reverse")
title('MVDR beam pattern')
xlabel('Frequency')
ylabel('angle')
zlabel('magnitude(dB)')

subplot(1,2,2)
P_lcmv=20*log10(B_lcmv);
mesh(f,theta,P_lcmv.');
set(gca,YDir="reverse")
title('LCMV beam pattern')
xlabel('Frequency')
ylabel('angle')
zlabel('magnitude(dB)')
