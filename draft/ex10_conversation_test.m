%% draft: simulation of two talker conversation
% Tongkai Li, 2023.08
clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'
addpath 'speech\male'
addpath 'speech\female'

%% load the audio signal
% source
[source1,fs]=readwav('ieee01m01.wav');
[source2,fs]=readwav('ieee01f02.wav');

L1=length(source1);
L2=length(source2);

%% parameters 
c = 340;                  % Sound velocity (m/s) 
num_mics=8;               % micphone number 
% distance between microphone (linear array) 
d=0.05; 
% Microphone array position [x y z] (m) 
position_mics = [4 3 2;  
                 4+d 3 2; 
                 4+2*d 3 2; 
                 4+3*d 3 2; 
                 4+4*d 3 2; 
                 4+5*d 3 2; 
                 4+6*d 3 2; 
                 4+7*d 3 2];  

% source location 
D=1; 
theta_source1 = 60;   % DOA 
position_source1 = D*[cos(theta_source1*pi/180) sin(theta_source1*pi/180) 0] + position_mics(1,:); 

theta_source2 = 120;   % DOA 
position_source2 = D*[cos(theta_source2*pi/180) sin(theta_source2*pi/180) 0] + position_mics(1,:); 

% interference location 
D_interf=1.5; 
theta_interf= 20; 
position_interf = D_interf*[cos(theta_interf*pi/180) sin(theta_interf*pi/180) 0] + position_mics(1,:); 

% distance matrix of mics
l_mn = zeros(num_mics, num_mics);
for i = 1:num_mics
    for j = 1:num_mics
        l_mn(i, j) = abs(i - j) * d;
    end
end


%% RIR generator
room = [8 6 5];                % Room dimensions [x y z] (m)
rev_time = 0.6;              % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir_desired1 = rir_generator(c, fs, position_mics, position_source1, room, rev_time, n);
rir_desired2 = rir_generator(c, fs, position_mics, position_source2, room, rev_time, n);
rir_undesired = rir_generator(c, fs, position_mics, position_interf, room, rev_time, n);

%% combine the rir to the speech signal
% get the received signal x in each microphone
desired_signal1=zeros(num_mics,length(source1));
desired_signal2=zeros(num_mics,length(source2));
t_internal=0;
sample_internal=t_internal*fs;

for i=1:num_mics
  desired_signal1(i,:)=filter(rir_desired1(i,:),1,source1)*10;
  desired_signal2(i,:)=filter(rir_desired2(i,:),1,source2)*10;
end

desired_signal=zeros(num_mics,L1+L2+sample_internal*2);
desired_signal(:,1:L1)=desired_signal1;
desired_signal(:,L1+sample_internal+1:L1+sample_internal+L2)=desired_signal2;
L=length(desired_signal);

% add interference
SIR=0;
rng(100);
% interf1 =0.01*randn(1,L); 
[interf1,fs2]=readwav('fan noise.wav');
[P,Q] = rat(fs/fs2);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);
interf1 = sqrt(sum(sum(desired_signal.^2)) / sum(sum(interf1.^2))/(10^(SIR/10))) * interf1;


for i=1:num_mics
  interf_signal(i,:)=filter(rir_undesired(i,:),1,interf1)*10;
end

% babble noise
SBNR=0;
[babble_noise,fs3]=readwav('babble noise.wav');
[P,Q] = rat(fs/fs3);
babble_noise = resample(babble_noise,P,Q);
theta_noise=[0,90,180,270];
D_noise=2;
% rir_noise=cell(length(theta_noise),1);
bnoise_all=zeros(num_mics,L);
for p=1:4
   position_noise = D_noise*[cos(theta_noise(p)/180*pi) sin(theta_noise(p)/180*pi) 0] + position_mics(1,:); 
   rir_noise = rir_generator(c, fs, position_mics, position_noise, room, rev_time, n);
   bnoise=babble_noise(p*100000+1:p*100000+L);
   bnoise =sqrt(sum(sum(desired_signal(1,:).^2)) / sum(sum(bnoise.^2))/(10^(SBNR/10)) ) * bnoise;

   for i=1:num_mics
      received_noise(i,:)=filter(rir_noise(i,:),1,bnoise)*10;
   end
   bnoise_all=bnoise_all+received_noise;
end

snoise=0*randn(num_mics,L);
noise=snoise + bnoise_all;

received_signal=desired_signal + interf_signal + noise;
% received_signal=desired_signal + interf_signal;

% soundsc(received_signal(1,:),fs)


%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc =512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

% segment signal into frames and do DFT
for i=1:num_mics
  [frame_total,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa'); 
  frame_desired=v_rfft(v_enframe(desired_signal(i,:),window,inc),N_window,2); 
  frame_interf=v_rfft(v_enframe(interf_signal(i,:),window,inc),N_window,2); 
  frame_noise=v_rfft(v_enframe(noise(i,:),window,inc),N_window,2); 

  total_frame(i,:,:)=frame_total;
  total_frame_desired(i,:,:)=frame_desired;
  total_frame_interf(i,:,:)=frame_interf;
  total_frame_noise(i,:,:)=frame_noise;

  n_frame=size(frame_total,1);
  n_fft=size(frame_total,2);
end

  % steering vector
h_interf = rir_generator(c, fs, position_mics, position_interf, room, 0, n);
temp=fft(h_interf,N_window,2);
vn=temp(:,1:n_fft);
vn=vn./vn(1,:);

h_desired1 = rir_generator(c, fs, position_mics, position_source1, room, 0, n);
h_desired2 = rir_generator(c, fs, position_mics, position_source2, room, 0, n);


freqVec=linspace(0, fs/2, n_fft);

% compute MVDR/LCMV weight
frame_speaker1=floor((L1+sample_internal)/inc);
for frm=1:n_frame
    % source DOA (a function of time t)
   if frm<=frame_speaker1
       theta_source=60;
       h_desired=h_desired1;
   else
       theta_source=120;
       h_desired=h_desired2;

   end
   temp=fft(h_desired,N_window,2);
   vs=temp(:,1:n_fft);
   vs=vs./vs(1,:);

  % compute MVDR/LCMV weight
  w_lcmv = zeros(num_mics, length(freqVec));
  for f = 1:length(freqVec)
     R = sinc(2*(f-1)/N_window*fs*l_mn/c);
     eps=1e-4;
     if rcond(R)<100
        R = R+eye(num_mics)*min(diag(R))*eps;
     end

     % LCMV
     C = [vs(:, f),vn(:,f)];
     denominator = C'/R*C;
     if rcond(denominator)<100
      denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
     end
     F=[1,0].';
     w_lcmv(:, f) = (R\C) / (denominator)*F;

  % compute output frame
     beam_frame_desired(frm,f)=w_lcmv(:,f)'*total_frame_desired(:,frm,f);
     beam_frame_interf(frm,f)=w_lcmv(:,f)'*total_frame_interf(:,frm,f);
     beam_frame_noise(frm,f)=w_lcmv(:,f)'*total_frame_noise(:,frm,f);

  end

  frm

end

output_desired=v_overlapadd(v_irfft(beam_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(beam_frame_interf,N_window,2),window,inc); 
output_noise=v_overlapadd(v_irfft(beam_frame_noise,N_window,2),window,inc); 

output_signal= output_desired + output_interf + output_noise;

soundsc(output_signal,fs)



%% calculate SIR improvement
sinr_in=sinr_test(desired_signal(1,:),interf_signal(1,:), noise(1,:),fs);
sinr_out=sinr_test(output_desired,output_interf,output_noise,fs);
sinr_imp=sinr_out-sinr_in

sir_in=sir_test(desired_signal(1,:),interf_signal(1,:),fs);
sir_out=sir_test(output_desired,output_interf,fs);
sir_imp=sir_out-sir_in


%% PESQ
% resample in 16000Hz
fs1=16000;
[P,Q] = rat(fs1/fs);
desired = resample(desired_signal(1,:),P,Q);
received = resample(received_signal(1,:),P,Q);
output = resample(output_signal,P,Q);

writewav(received,16000,'received_signal.wav')
writewav(desired,16000,'desired_signal.wav')
writewav(output,16000,'output_signal.wav')
score1=pesq('desired_signal.wav','received_signal.wav')
score2=pesq('desired_signal.wav','output_signal.wav')
del_score=score2-score1

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(received_signal(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(received_signal(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);


title('Received signal at reference Mic','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis([-120 -40])
lim=caxis;

subplot(3,1,3)
plot(t0,received_signal(1,:));
ylim([-0.8,0.8])
xlim([0,length(received_signal(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_signal)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_signal,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);
title('Output signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')

subplot(3,1,3)
plot(t0,output_signal);
ylim([-0.8,0.8])
xlim([0,length(output_signal)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

t0=[1:1:length(desired_signal(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(desired_signal(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Desired component in received signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,desired_signal(1,:));
ylim([-0.8,0.8])
xlim([0,length(desired_signal(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_desired)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_desired,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Desired component in output','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,output_desired);
ylim([-0.8,0.8])
xlim([0,length(output_desired)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(interf_signal(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(interf_signal(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Interference component in received signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,interf_signal(1,:));
ylim([-0.8,0.8])
xlim([0,length(interf_signal(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_interf)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_interf,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Interference component in output','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,output_interf);
ylim([-0.8,0.8])
xlim([0,length(output_interf)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(noise(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(noise(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Noise component in received signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,noise(1,:));
ylim([-0.8,0.8])
xlim([0,length(noise(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_noise)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_noise,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Noise component in LCMV output','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,output_noise);
ylim([-0.8,0.8])
xlim([0,length(output_noise)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')


%%
t_frame=N_window/fs;
t1=t_frame/2*[1:1:n_frame];
figure;
subplot(2,1,1)
t0=[1:1:length(desired_signal(1,:))]/fs;
plot(t0,desired_signal(1,:))
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')
title('Speech waveform','FontSize',17,'FontName','Times New Roman')

subplot(2,1,2)
plot(t1,Pt1,'LineWidth',1.2);
hold on
plot(t1,Pt2,'LineStyle','-.','LineWidth',1.2);
hold on
plot(t1,Pn1,'LineStyle',':','LineWidth',1.2);
hold on
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Gain (dB)','FontSize',15,'FontName','Times New Roman')
title('Gains at the direction of beam and null','FontSize',17,'FontName','Times New Roman')
