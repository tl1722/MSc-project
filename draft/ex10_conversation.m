%% ex5: LCMV and MVDR in frequency domain (with RIR)
% consider interference signal (babble noise) and wgn noise

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% load the audio signal
% source
[source1,fs]=readwav('S_01_01.wav');
[source2,fs]=readwav('S_01_02.wav');

L1=length(source1);
L2=length(source2);

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

% source location 
D=1; 
theta_source1 = 60*pi/180;   % DOA 
position_source1 = D*[cos(theta_source1) sin(theta_source1) 0] + position_mics(1,:); 

theta_source2 = 120*pi/180;   % DOA 
position_source2 = D*[cos(theta_source2) sin(theta_source2) 0] + position_mics(1,:); 

% interference location 
D_interf=3; 
theta_interf= 20*pi/180; 
position_interf = D_interf*[cos(theta_interf) sin(theta_interf) 0] + position_mics(1,:); 

%% RIR generator
room = [8 6 5];                % Room dimensions [x y z] (m)
rev_time = 0;              % Reverberation time (s)
m = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir_desired1 = rir_generator(c, fs, position_mics, position_source1, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
rir_desired2 = rir_generator(c, fs, position_mics, position_source2, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
rir_undesired = rir_generator(c, fs, position_mics, position_interf, room, rev_time, m, mtype, order, dim, orientation, hp_filter);

%% combine the rir to the speech signal
% get the received signal x in each microphone
desired_signal1=zeros(num_mics,length(source1));
desired_signal2=zeros(num_mics,length(source2));
t_internal=0.5;
sample_internal=t_internal*fs;

for i=1:num_mics
  desired_signal1(i,:)=filter(rir_desired1(i,:),1,source1)*10;
  desired_signal2(i,:)=filter(rir_desired2(i,:),1,source2)*10;
end

conversation=zeros(num_mics,L1+L2+sample_internal*2);
conversation(:,1:L1)=desired_signal1;
conversation(:,L1+sample_internal+1:L1+sample_internal+L2)=desired_signal2;
L=length(conversation);

% add interference
SIR=0;
rng(100);
% interf1 =0.01*randn(1,L); 
[interf1,fs2]=readwav('babble noise.wav');
[P,Q] = rat(fs/fs2);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);
% interf1 = sqrt(sum(sum(conversation.^2)) / sum(sum(interf1.^2))/(10^(SIR/10))) * interf1;


for i=1:num_mics
  interf_signal(i,:)=filter(rir_undesired(i,:),1,interf1)*10;
end

noise=0.0001*randn(num_mics,L);
received_signal=conversation + interf_signal + noise;

% soundsc(received_signal(1,:),fs)


%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc =256;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

% segment signal into frames and do DFT
for i=1:num_mics
  [frame_total,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa'); 
  total_frame(i,:,:)=frame_total;
  n_frame=size(frame_total,1);
  n_fft=size(frame_total,2);
end

%% calculate w using data independent model
theta_noise=(0:10:360)'; 
R=zeros(num_mics,num_mics);
for k=1:length(theta_noise)
  position_noise = D*[cos(theta_noise(k)/180*pi) sin(theta_noise(k)/180*pi) 0] + position_mics(1,:);   
  order=0;
  rir_diffuse = rir_generator(c, fs, position_mics, position_noise, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
  atf=fft(rir_diffuse,n_fft,2);
  atf=atf./atf(1,:);
%   R=R+rir_diffuse*rir_diffuse';
  R=R+atf*atf';
end
R=R/length(theta_noise);

freqVec=linspace(0, fs/2, n_fft);
% steering vector
vs1 = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source1)/c); % signal direction vector
vs2 = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source2)/c); % signal direction vector
vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interf)/c); % interference direction vector

% compute MVDR/LCMV weight
w1_lcmv = zeros(num_mics, length(freqVec));
w2_lcmv = zeros(num_mics, length(freqVec));
eps=1e-4;

for f = 1:length(freqVec)
%    R=sin(2*pi*f*l_mn/c/n_fft)/(2*pi*f*l_mn/c/n_fft);

   % LCMV
   C = [vs1(:, f), vs2(:, f), vn(:,f)];
   denominator = C'/(R)*C;
   if rcond(denominator)<100
      denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
   end

   F1=[1,0,0].';
   w1_lcmv(:, f) = ((R)\C) / (denominator)*F1;

   F2=[0,1,0].';
   w2_lcmv(:, f) = ((R)\C) / (denominator)*F2;

end 

frame_s1=floor((L1+sample_internal)/inc);
frame_s2=n_frame-frame_s1;

% compute output frame
for i=1:num_mics
  [frame_desired,T,WS]=v_enframe(conversation(i,:),window,1/overlap_factor,'fa'); 
  [frame_interf,T,WS]=v_enframe(interf_signal(i,:),window,1/overlap_factor,'fa'); 
  [frame_noise,T,WS]=v_enframe(noise(i,:),window,1/overlap_factor,'fa'); 
  for n=1:n_frame
      if n<frame_s1
        beam_frame_desired(n,:,i)=w1_lcmv(i,:).*(frame_desired(n,:));
        beam_frame_interf(n,:,i)=w1_lcmv(i,:).*(frame_interf(n,:));
        beam_frame_noise(n,:,i)=w1_lcmv(i,:).*(frame_noise(n,:));
      else
        beam_frame_desired(n,:,i)=w2_lcmv(i,:).*(frame_desired(n,:));
        beam_frame_interf(n,:,i)=w2_lcmv(i,:).*(frame_interf(n,:));
        beam_frame_noise(n,:,i)=w2_lcmv(i,:).*(frame_noise(n,:));
      end
  end
end

% sum Y in every mics
sum_frame_desired=sum(beam_frame_desired,3);
sum_frame_interf=sum(beam_frame_interf,3);
sum_frame_noise=sum(beam_frame_noise,3);


% use overlapadd to reconstruct the signal y
output_desired=v_overlapadd(v_irfft(sum_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(sum_frame_interf,N_window,2),window,inc); 
output_noise=v_overlapadd(v_irfft(sum_frame_noise,N_window,2),window,inc); 

output_signal= output_desired + output_interf + output_noise;

% soundsc(received_signal(1,:),fs)
soundsc(output_signal,fs)

%% calculate SIR improvement
sir_in=sir_test(conversation(1,:),interf_signal(1,:),fs)
sir_out=sir_test(output_desired,output_interf,fs)

sir_imp=sir_out-sir_in

%% PESQ
% resample in 16000Hz
fs1=16000;
[P,Q] = rat(fs1/fs);
desired = resample(conversation(1,:),P,Q);
received = resample(received_signal(1,:),P,Q);
output = resample(output_signal,P,Q);

writewav(received,16000,'received_signal.wav')
writewav(desired,16000,'desired_signal.wav')
writewav(output,16000,'output_signal.wav')
score1=pesq('desired_signal.wav','received_signal.wav')
score2=pesq('desired_signal.wav','output_signal.wav')
del_score=score2-score1

figure;
[s,f,t] = spectrogram(received_signal(1,:),window,inc,N_window,fs,'yaxis');
subplot(2,4,1)
spectrogram(received_signal(1,:),window,inc,N_window,fs,'yaxis');
title('Received signal')
lim=caxis;
subplot(2,4,2)
spectrogram(output_signal,window,inc,N_window,fs,'yaxis');
title('beamformed signal')
caxis(lim);
subplot(2,4,3)
spectrogram(interf_signal(1,:),window,inc,N_window,fs,'yaxis');
title('interference signal')
caxis(lim);
subplot(2,4,4)
spectrogram(output_interf,window,inc,N_window,fs,'yaxis');
title('interference output')
caxis(lim);
subplot(2,4,5)
plot(received_signal(1,:))
ylim([-0.5 0.5])
subplot(2,4,6)
plot(output_signal)
ylim([-0.5 0.5])
subplot(2,4,7)
plot(interf_signal(1,:))
ylim([-0.5 0.5])
subplot(2,4,8)
plot(output_interf)
ylim([-0.5 0.5])



% % spectrogram
% figure;
% subplot(1,3,1)
% [s,f,t] = spectrogram(interf_signal(1,:),window,inc,N_window,fs) ;
% spectrogram(conversation(1,:),window,inc,N_window,fs,'yaxis');
% title('received signal in mic 1')
% lim=caxis;
% subplot(1,3,2)
% spectrogram(interf_signal(1,:),window,inc,N_window,fs,'yaxis');
% title('total received signal in mic 1')
% caxis(lim);
% subplot(1,3,3)
% spectrogram(output_interf,window,inc,N_window,fs,'yaxis') ;
% title('beamformed signal')
% caxis(lim);
% 
% 
% % plot the 3D beam pattern
% figure;
% theta = 180 * [-1:0.001:1];
% B_lcmv=zeros(n_fft,length(theta));
% for i=1:length(f)
%   v = exp(-1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
%   B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);
% 
% end
% 
% P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));
% mesh(f,theta,P_lcmv.');
% set(gca,YDir="reverse")
% title('LCMV beam pattern')
% xlabel('Frequency')
% ylabel('angle')
% zlabel('magnitude(dB)')
% 
% figure;
% k=1335;
% plot(freqVec,P_lcmv(:,k))
% hold on
% k1=1112;
% plot(freqVec,P_lcmv(:,k1))
% title('Beam pattern of LCMV')
% xlabel('Frequency')
% ylabel('dB')
% grid on
% 
% % 2D plot
