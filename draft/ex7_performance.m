%% ex5: LCMV and MVDR in frequency domain (with RIR)
%% no noise, analysis performance using SIR

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% load the audio signal
% source
[source,fs]=readwav('S_01_01.wav');
L=length(source);
% interference
rng(100);
% interf1 = 0.01*randn(L,1); 
[interf1,fs1]=readwav('babble noise.wav');
[P,Q] = rat(fs/fs1);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);

%% parameters
c = 340;                  % Sound velocity (m/s)
num_mics=8;               % micphone number
% distance between microphone (linear array)
d=0.05;
% Microphone array position [x y z] (m)
position_mics = [1 1.5 2; 
                 1+d 1.5 2;
                 1+2*d 1.5 2;
                 1+3*d 1.5 2;
                 1+4*d 1.5 2;
                 1+5*d 1.5 2;
                 1+6*d 1.5 2;
                 1+7*d 1.5 2]; 

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
room = [5 4 6];                % Room dimensions [x y z] (m)
rev_time = 0;             % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=theta_source;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, n, mtype, order, dim, orientation, hp_filter);
rir_undesired = rir_generator(c, fs, position_mics, position_interference, room, rev_time, n, mtype, order, dim, orientation, hp_filter);

%% combine the rir to the speech signal
% get the received signal x in each microphone
desired_signal=zeros(num_mics,length(source));
interference_signal_1=zeros(num_mics,length(source));

for i=1:num_mics
  desired_signal(i,:)=filter(rir_desired(i,:),1,source)*10;
  interference_signal_1(i,:)=filter(rir_undesired(i,:),1,interf1)*10;
end

received_signal = desired_signal + interference_signal_1;
% sound(received_signal_1(1,:),fs)

%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc = 512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = hamming(N_window,'periodic');

% segment signal into frames and do DFT
for i=1:num_mics
  [frame_desired,T,WS]=v_enframe(desired_signal(i,:),window,1/overlap_factor,'fa'); 
  [frame_interf,T,WS]=v_enframe(interference_signal_1(i,:),window,1/overlap_factor,'fa');
  total_frame_desired(i,:,:)=frame_desired;
  total_frame_interf(i,:,:)=frame_interf;
  n_frame=size(frame_desired,1);
  n_fft=size(frame_desired,2);
end

% steering vector
freqVec=linspace(0, fs/2, n_fft);
vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source)/c); % signal direction vector
vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interference)/c); % interference direction vector

% compute MVDR/LCMV weight
w_mvdr = zeros(num_mics, length(freqVec));
w_lcmv = zeros(num_mics, length(freqVec));
eps=1e-4;
for f = 1:length(freqVec)  
   Phi_n=1/n_frame * (total_frame_interf(:,:,f)*total_frame_interf(:,:,f)');
   
   % LCMV
   C = [vs(:, f),vn(:,f)];
   denominator = C'/Phi_n*C;
   if rcond(denominator)<eps
     denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
   end
   F=[1,0].';
   w_lcmv(:, f) = (Phi_n\C) / (denominator)*F;

   % MVDR
   vf=vs(:, f);
   w_mvdr(:, f) = (Phi_n\vf) / (vf'/Phi_n*vf);

end

for i=1:num_mics
  [frame_desired,T,WS]=v_enframe(desired_signal(i,:),window,1/overlap_factor,'fa'); 
  [frame_interf,T,WS]=v_enframe(interference_signal_1(i,:),window,1/overlap_factor,'fa'); 
  for n=1:n_frame
      % mvdr output frame
%       beam_frame_desired(n,:,i)=w_mvdr(i,:).*(frame_desired(n,:));
%       beam_frame_interf(n,:,i)=w_mvdr(i,:).*(frame_interf(n,:));

      % lcmv output frame
      beam_frame_desired(n,:,i)=w_lcmv(i,:).*(frame_desired(n,:));
      beam_frame_interf(n,:,i)=w_lcmv(i,:).*(frame_interf(n,:));
% 
  end
end

% sum Y in every mics
sum_frame_desired=sum(beam_frame_desired,3);
sum_frame_interf=sum(beam_frame_interf,3);


% use overlapadd to reconstruct the signal y
output_desired=v_overlapadd(v_irfft(sum_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(sum_frame_interf,N_window,2),window,inc); 

output_signal=output_desired+output_interf;
soundsc(output_signal,fs)
% sound(output_interf,fs)

%% calculate SIR improvement
% P_sd=(source'*source)/length(source);
% P_si=(interf1'*interf1)/length(interf1);
P_sd=(desired_signal(1,:)*desired_signal(1,:)')/length(desired_signal(1,:));
P_si=(interference_signal_1(1,:)*interference_signal_1(1,:)')/length(interference_signal_1(1,:));
sir_in=10*log10(P_sd/P_si);

P_yd=(output_desired'*output_desired)/length(output_desired);
P_yi=(output_interf'*output_interf)/length(output_interf);
sir_o=10*log10(P_yd/P_yi);

sir_imp=sir_o-sir_in;

%% PESQ
writewav(desired_signal(1,:),16000,'desired_signal.wav')
writewav(output_signal,16000,'output_signal.wav')
score=pesq('desired_signal.wav','output_signal.wav')

%% plot
% plot waveform
figure;
subplot(4,1,1)
plot(source);
title('source signal')
subplot(4,1,2)
plot(desired_signal(1,:));
title('desired signal in Mic 1')
subplot(4,1,3)
plot(interference_signal_1(1,:));
title('interference signal in Mic 1')
% sound(desired_signal(1,:),fs)

% writewav(output_signal,fs,'output_signal_01.wav')
subplot(4,1,4);
plot(desired_signal(1,:));
hold on
plot(output_signal)
title('Beamformed signal (LCMV)')
legend('desired signal in mics 1','output signal')


% spectrogram
figure;
subplot(1,2,1)
[s,f,t] = spectrogram(received_signal(1,:),window,inc,N_window,fs) ;
spectrogram(desired_signal(1,:),window,inc,N_window,fs,'yaxis');
title('desired signal in mics 1')
subplot(1,2,2)
spectrogram(output_signal,window,inc,N_window,fs,'yaxis') ;
title('output signal')


% plot the 3D beam pattern
theta = 180 * [-1:0.01:1];
B_mvdr=zeros(n_fft,length(theta));
B_lcmv=zeros(n_fft,length(theta));
for i=1:length(f)
  v = exp(-1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
  B_mvdr(i,:) = abs(w_mvdr(:,i)' * v);
  B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);

end

figure;
subplot(1,2,1)
P_mvdr=20*log10(B_mvdr./max(B_mvdr,[],2));
mesh(f,theta,P_mvdr.');
set(gca,YDir="reverse")
title('MVDR beam pattern')
xlabel('Frequency')
ylabel('angle')
zlabel('magnitude(dB)')

subplot(1,2,2)
P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));
mesh(f,theta,P_lcmv.');
set(gca,YDir="reverse")
title('LCMV beam pattern')
xlabel('Frequency')
ylabel('angle')
zlabel('magnitude(dB)')

% figure;
% k=70;
% plot(theta,20*log10(B_mvdr(k,:)/max(B_mvdr(k,:))))
% hold on
% title('Beam pattern of LCMV')
% xlabel('degree')
% ylabel('dB')
% grid on
