%% chap 2: comparison of mvdr and lcmv (draft)
% output: spectrograms and waveforms, beam pattern of MVDR and LCMV
% TONGKAI LI, 2023.08

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'
addpath 'speech\male'
addpath 'speech\female'

set(0,'DefaultAxesFontSize',13)
set(0,'defaultAxesFontName','Times New Roman')
set(0,'DefaultLineLineWidth',1);

%% load the audio signal
% source
[source,fs]=audioread('ieee01m10.wav');
L=length(source);
% interference
SIR=0;
rng(100);
% interf1 =0*randn(1,L);
[interf1,fs2]=readwav('ieee58f03.wav');
[P,Q] = rat(fs/fs2);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);
interf1 = sqrt( sum(sum(source.^2)) / sum(sum(interf1.^2))/(10^(SIR/10)) ) * interf1;

%% parameters 
c = 340;                  % Sound velocity (m/s) 
num_mics=8;               % micphone number 
% distance between microphone (linear array) 
d=0.04; 
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
theta_source = 60;   % DOA 
position_source = D*[cos(theta_source/180*pi) sin(theta_source/180*pi) 0] + position_mics(1,:); 
TDOA = sqrt(sum((bsxfun(@minus, position_source, position_mics)).^2, 2))/c; 

% interference location 
D_interf=1; 
theta_interf=20; 
position_interf = D_interf*[cos(theta_interf/180*pi) sin(theta_interf/180*pi) 0] + position_mics(1,:); 

% distance matrix of mics
l_mn = zeros(num_mics, num_mics);
for i = 1:num_mics
    for j = 1:num_mics
        l_mn(i, j) = abs(i - j) * d;
    end
end


%% RIR generator
room = [8 6 5];                % Room dimensions [x y z] (m)
rev_time = 0.6;             % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, n, mtype, order, dim, orientation, hp_filter);
rir_interf = rir_generator(c, fs, position_mics, position_interf, room, rev_time, n, mtype, order, dim, orientation, hp_filter);

%% combine the rir to the speech signal
% get the received signal x in each microphone
desired_signal=zeros(num_mics,length(source));
interference_signal=zeros(num_mics,length(source));

for i=1:num_mics
  desired_signal(i,:)=filter(rir_desired(i,:),1,source)*10;
  interference_signal(i,:)=filter(rir_interf(i,:),1,interf1)*10;
end

% babble noise
SBNR=20;
SSNR=40;
[babble_noise,fs3]=readwav('babble noise.wav');
[P,Q] = rat(fs/fs3);
babble_noise = resample(babble_noise,P,Q);
theta_noise=[0,90,180,270];
D_noise=3;
% rir_noise=cell(length(theta_noise),1);
bnoise_all=zeros(num_mics,L);
for p=1:4
   position_noise = D_noise*[cos(theta_noise(p)/180*pi) sin(theta_noise(p)/180*pi) 0] + position_mics(1,:); 
   rir_noise = rir_generator(c, fs, position_mics, position_noise, room, rev_time, n);
   bnoise=babble_noise(p*100000+1:p*100000+L);
   for i=1:num_mics
      received_noise(i,:)=filter(rir_noise(i,:),1,bnoise)*10;
   end
   bnoise_all=bnoise_all+received_noise;
end
% bnoise_all =sqrt( sum(sum(source.^2)) / sum(sum(bnoise_all.^2))/(10^(SBNR/10)) ) * bnoise_all;

% sensor noise
% snoise=0.01*randn(num_mics,L);
snoise=wgn(num_mics,L,-40);
% noise=snoise+bnoise_all;
noise=snoise;

received_signal = desired_signal + interference_signal+ noise;
% soundsc(received_signal(1,:),fs)
received_passthrough=sum(received_signal,1)/num_mics;

%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc =512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

% segment signal into frames and do DFT
for i=1:num_mics
  frame_total=v_rfft(v_enframe(received_signal(i,:),window,inc),N_window,2); 
  frame_desired=v_rfft(v_enframe(desired_signal(i,:),window,inc),N_window,2); 
  frame_interf=v_rfft(v_enframe(interference_signal(i,:),window,inc),N_window,2); 
  frame_noise=v_rfft(v_enframe(noise(i,:),window,inc),N_window,2); 

  total_frame(i,:,:)=frame_total;
  total_frame_desired(i,:,:)=frame_desired;
  total_frame_interf(i,:,:)=frame_interf;
  total_frame_noise(i,:,:)=frame_noise;

  n_frame=size(frame_total,1);
  n_fft=size(frame_total,2);
end

%% calculate R using data independent model
d_desired = rir_generator(c, fs, position_mics, position_source, room, 0, n);
d_interf = rir_generator(c, fs, position_mics, position_interf, room, 0, n);
temp=fft(d_desired,1024,2);
vs=temp(:,1:n_fft);
vs=vs./vs(1,:);
temp=fft(d_interf,1024,2);
vn=temp(:,1:n_fft);
vn=vn./vn(1,:);

freqVec=linspace(0, fs/2, n_fft);
alpha=0.98;
% steering vector
% vs = exp(1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source*pi/180)/c); % signal direction vector
% vn = exp(1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interf*pi/180)/c); % interference direction vector

% compute MVDR/LCMV weight
for frm=1:n_frame
  % compute MVDR/LCMV weight
  w_lcmv = zeros(num_mics, length(freqVec));
  for f = 1:length(freqVec)
%    R=1/n_frame * (total_frame(:,:,f)*total_frame(:,:,f)');

     if frm == 1
        R = total_frame(:,frm,f) * total_frame(:,frm,f)';
     else
        R = alpha * R + (1-alpha) * total_frame(:,frm,f) * total_frame(:,frm,f)';
     end
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
     total_w(:,f,frm)=w_lcmv(:, f);

     % MVDR
     vf=vs(:, f);
     w_mvdr(:, f) = ((R)\vf) / (vf'/(R)*vf);
     total_w_mvdr(:,f,frm)=w_mvdr(:, f);

     beam_frame_desired(frm,f)=w_lcmv(:,f)'*total_frame_desired(:,frm,f);
     beam_frame_interf(frm,f)=w_lcmv(:,f)'*total_frame_interf(:,frm,f);
     beam_frame_noise(frm,f)=w_lcmv(:,f)'*total_frame_noise(:,frm,f);

     beam_frame_desired1(frm,f)=w_mvdr(:,f)'*total_frame_desired(:,frm,f);
     beam_frame_interf1(frm,f)=w_mvdr(:,f)'*total_frame_interf(:,frm,f);
     beam_frame_noise1(frm,f)=w_mvdr(:,f)'*total_frame_noise(:,frm,f);

  end

  frm
end
% output=v_overlapadd(v_irfft(beam_frame,N_window,2),window,inc); 

total_w=mean(total_w,3);
total_w_mvdr=mean(total_w_mvdr,3);
% fir_coef=circshift(v_irfft(total_w,N_window,2),inc,2);
% figure;
% plot(fir_coef(1,:));

output_desired=v_overlapadd(v_irfft(beam_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(beam_frame_interf,N_window,2),window,inc); 
output_noise=v_overlapadd(v_irfft(beam_frame_noise,N_window,2),window,inc); 

output_signal= output_desired + output_interf + output_noise;

% mvdr
output_desired1=v_overlapadd(v_irfft(beam_frame_desired1,N_window,2),window,inc); 
output_interf1=v_overlapadd(v_irfft(beam_frame_interf1,N_window,2),window,inc); 
output_noise1=v_overlapadd(v_irfft(beam_frame_noise1,N_window,2),window,inc); 

output_signal1= output_desired1 + output_interf1 + output_noise1;

% soundsc(output_signal,fs)
% soundsc(received_signal(1,:),fs)

%%
sir_in=sir_test(desired_signal(1,:),interference_signal(1,:),fs);
sir_lcmv=sir_test(output_desired, output_interf,fs);
sir_mvdr=sir_test(output_desired1,output_interf1,fs);

snr_in=sir_test(desired_signal(1,:),noise(1,:),fs);
snr_lcmv=sir_test(output_desired, output_noise,fs);
snr_mvdr=sir_test(output_desired1,output_noise1,fs);

 %% plot
% plot waveform

% spectrogram
figure;
subplot(2,4,1)
[s,f,t] = spectrogram(interference_signal(1,:),window,inc,N_window,fs) ;
spectrogram(source,window,inc,N_window,fs,'yaxis');
title('Target signal')
lim=caxis;
subplot(2,4,2)
spectrogram(received_signal(1,:),window,inc,N_window,fs,'yaxis');
title('received signal in reference microphone')
caxis(lim);
subplot(2,4,3)
spectrogram(output_signal1,window,inc,N_window,fs,'yaxis') ;
title('MVDR beamformed signal')
caxis(lim);
subplot(2,4,4)
spectrogram(output_signal,window,inc,N_window,fs,'yaxis') ;
title('LCMV beamformed signal')
caxis(lim);

subplot(2,4,5)
plot(source);
subplot(2,4,6)
plot(received_signal(1,:));
subplot(2,4,7)
plot(output_signal1);
subplot(2,4,8);
plot(output_signal);

figure;
subplot(2,4,1)
[s,f,t] = spectrogram(interference_signal(1,:),window,inc,N_window,fs) ;
spectrogram(source,window,inc,N_window,fs,'yaxis');
title('Target signal')
lim=caxis;
subplot(2,4,2)
spectrogram(received_signal(1,:),window,inc,N_window,fs,'yaxis');
title('received signal in reference microphone')
caxis(lim);
subplot(2,4,3)
spectrogram(output_interf1,window,inc,N_window,fs,'yaxis') ;
title('MVDR beamformed signal')
caxis(lim);
subplot(2,4,4)
spectrogram(output_interf,window,inc,N_window,fs,'yaxis') ;
title('LCMV beamformed signal')
caxis(lim);

subplot(2,4,5)
plot(source);
ylim([-0.4,0.4])
subplot(2,4,6)
plot(interference_signal(1,:));
ylim([-0.4,0.4])
subplot(2,4,7)
plot(output_interf1);
ylim([-0.4,0.4])
subplot(2,4,8);
plot(output_interf);
ylim([-0.4,0.4])

% plot the 3D beam pattern
theta = 180 * [-1:0.01:1];
B_mvdr=zeros(n_fft,length(theta));
B_lcmv=zeros(n_fft,length(theta));
for i=1:length(f)
  v = exp(1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
  B_mvdr(i,:) = abs(total_w_mvdr(:,i)' * v);
  B_lcmv(i,:) = abs(total_w(:,i)' * v);

end

% 2D plot
figure;
P_mvdr=20*log10(B_mvdr./max(B_mvdr,[],2));
imagesc(theta',f',P_mvdr);
set(gca,'YDir','normal')
% set(gca,'XTickLabelRotation',23)
title('MVDR beam pattern')
xlabel('DOA (Degree)')
ylabel('Frequency (Hz)')
colorbar
caxis([-50 0])
lim2=caxis;
xline(66,'blue',LineWidth=2);
xline(23,'red',LineWidth=2);

figure;
P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));
imagesc(theta',f',P_lcmv);
set(gca,'YDir','normal')
title('LCMV beam pattern')
xlabel('DOA (Degree)')
ylabel('Frequency (Hz)')
colorbar
caxis(lim2);
hold on
xline(66,'blue',LineWidth=2);
xline(23,'red',LineWidth=2);


% 3D plot
figure;
mesh(f,theta,P_mvdr.');
set(gca,YDir="reverse")
% set(gca,'XTickLabelRotation',23)

title('MVDR beam pattern')
ylabel('DOA (Degree)','FontSize',15,'FontName','Times New Roman')
xlabel('Frequency (Hz)','FontSize',15,'FontName','Times New Roman')
zlabel('Gain (dB)','FontSize',15,'FontName','Times New Roman')

figure;
mesh(f,theta,P_lcmv.');
set(gca,YDir="reverse")
% set(gca,'XTickLabelRotation',23)

title('LCMV beam pattern','FontSize',17,'FontName','Times New Roman')
ylabel('DOA (Degree)','FontSize',15,'FontName','Times New Roman')
xlabel('Frequency (Hz)','FontSize',15,'FontName','Times New Roman')
zlabel('Gain (dB)','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t0=[1:1:length(interf1)]/fs;
% figure;
% subplot(3,1,[1,2])
% [s,f,t] = spectrogram(interference_signal(1,:),window,inc,N_window,fs) ;
% ylabel(colorbar,[])
% 
% spectrogram(interf1,window,inc,N_window,fs,'yaxis');
% ylabel(colorbar,[])
% 
% lim=caxis;
% title('Interference signal','FontSize',17,'FontName','Times New Roman')
% ylabel('Frequency (kHz)')
% subplot(3,1,3)
% plot(t0,interf1);
% ylim([-0.4,0.4])
% xlim([0,length(interf1)/fs])
% xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
% ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
t0=[1:1:length(source)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(source,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

title('Target Signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis([-120 -40])
lim=caxis;

subplot(3,1,3)
plot(t0,source);
ylim([-0.4,0.4])
xlim([0,length(source)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(source)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(received_signal(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

title('Received signal at reference Mic','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis(lim);
subplot(3,1,3)
plot(t0,received_signal(1,:));
ylim([-0.4,0.4])
xlim([0,length(source)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_signal1)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_signal1,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz');
xlabel([]);
title('MVDR output signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis(lim);
subplot(3,1,3)
plot(t0,output_signal1);
ylim([-0.4,0.4])
xlim([0,length(output_signal1)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_signal)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_signal,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz');
xlabel([]);
title('LCMV output signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis(lim);
subplot(3,1,3)
plot(t0,output_signal);
ylim([-0.4,0.4])
xlim([0,length(output_signal)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(interference_signal(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(interference_signal(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz');
xlabel([]);
title('Interference component in received signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis(lim);
subplot(3,1,3)
plot(t0,interference_signal(1,:));
ylim([-0.4,0.4])
xlim([0,length(interference_signal(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_signal)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_interf1,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz');
xlabel([]);
title('MVDR output interference component','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis(lim);
subplot(3,1,3)
plot(t0,output_interf1);
ylim([-0.4,0.4])
xlim([0,length(output_signal)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[1:1:length(output_signal)]/fs;
figure('units','normalized');
subplot(3,1,[1,2])
spectrogram(output_interf,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz');
xlabel([]);
title('Interference component in LCMV output','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis(lim);
subplot(3,1,3)
plot(t0,output_interf);
ylim([-0.4,0.4])
xlim([0,length(output_signal)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')


