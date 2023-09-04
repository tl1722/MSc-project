%% ex5: LCMV and MVDR in frequency domain (with RIR)
% consider interference signal (babble noise) and wgn noise

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'
addpath 'F:\恺\project\code\new code'

%% load the audio signal
% source
[source,fs]=readwav('S_01_02.wav');
L=length(source);
% interference
SIR=0;
rng(100);
% interf1 =0.01*randn(1,L); 
[interf1,fs2]=readwav('S_01_01.wav');
[P,Q] = rat(fs/fs2);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);
interf1 = sqrt(sum(sum(source.^2)) / sum(sum(interf1.^2))/(10^(SIR/10))) * interf1;

% 
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

% position_mics = [1 1.5 2;  1+d 1.5 2; 1+2*d 1.5 2; 1+3*d 1.5 2; 
%                  1+4*d 1.5 2; 1+5*d 1.5 2; 1+6*d 1.5 2; 1+7*d 1.5 2;
%                  1+8*d 1.5 2; 1+9*d 1.5 2; 1+10*d 1.5 2; 1+11*d 1.5 2;
%                  1+12*d 1.5 2; 1+13*d 1.5 2; 1+14*d 1.5 2; 1+15*d 1.5 2]; 
% source location 
D=1; 
theta_source = 60*pi/180;   % DOA 
position_source = D*[cos(theta_source) sin(theta_source) 0] + position_mics(1,:); 
% TDOA = sqrt(sum((bsxfun(@minus, position_source, position_mics)).^2, 2))/c; 

% interference location 
D_interf=2; 
theta_interf= 20*pi/180; 
position_interf = D_interf*[cos(theta_interf) sin(theta_interf) 0] + position_mics(1,:); 

%% RIR generator
room = [8 6 5];                % Room dimensions [x y z] (m)
rev_time = 0.3;             % Reverberation time (s)
m = 4096;                   % Number of samples

rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, m);
rir_undesired = rir_generator(c, fs, position_mics, position_interf, room, rev_time, m);


%% combine the rir to the speech signal
% get the received signal x in each microphone
desired_signal=zeros(num_mics,length(source));
interference_signal=zeros(num_mics,length(source));

for i=1:num_mics
  desired_signal(i,:)=filter(rir_desired(i,:),1,source)*10;
  interference_signal(i,:)=filter(rir_undesired(i,:),1,interf1)*10;
end

noise=0*randn(num_mics,L);

received_signal = desired_signal + interference_signal + noise;
% soundsc(received_signal(1,:),fs)


%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc = 512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

% segment signal into frames and do DFT
for i=1:num_mics
  frame_total=v_rfft(v_enframe(received_signal(i,:),window,inc),N_window,2); 
%   [frame_total,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa'); 
  total_frame(i,:,:)=frame_total;
  n_frame=size(frame_total,1);
  n_fft=size(frame_total,2);
end

% steering vector
freqVec=linspace(0, fs/2, n_fft);

% rir_desired = rir_generator(c, fs, position_mics, position_source, room, 0, m);
% rir_undesired = rir_generator(c, fs, position_mics, position_interf, room, 0, m);
% temp=fft(rir_desired,1024,2);
% vs=temp(:,1:n_fft);
% vs=vs./vs(1,:);
% temp=fft(rir_undesired,1024,2);
% vn=temp(:,1:n_fft);
% vn=vn./vn(1,:);

theta_change=0;
theta_actual=theta_source+ theta_change*pi/180;
interf_actual=theta_interf+theta_change*pi/180;

vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_actual)/c); % signal direction vector
vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(interf_actual)/c); % interference direction vector


% compute MVDR/LCMV weight
w_mvdr = zeros(num_mics, length(freqVec));
w_lcmv = zeros(num_mics, length(freqVec));
eps=1e-4;
alpha=0.9;

for f = 1:length(freqVec)
   R=1/n_frame * (total_frame(:,:,f)*total_frame(:,:,f)');
%   for frm=1:n_frame   
%      if frm == 1
%         R = total_frame(:,frm,f) * total_frame(:,frm,f)';
%      else
%         R = alpha * R + (1-alpha) * total_frame(:,frm,f) * total_frame(:,frm,f)';
%      end
%   end

   % LCMV
   C = [vs(:, f),vn(:,f)];
   denominator = C'/(R)*C;
   if rcond(denominator)<100
     denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
   end
   F=[1,0].';
   w_lcmv(:, f) = ((R)\C) / (denominator)*F;
%    w_lcmv(:, f) = w_lcmv(:, f)./norm(w_lcmv(:,f));
   
   % MVDR
   vf=vs(:, f);
   w_mvdr(:, f) = ((R)\vf) / (vf'/(R)*vf);
end

% compute output frame
for i=1:num_mics
  frame_desired=v_rfft(v_enframe(desired_signal(i,:),window,inc),N_window,2); 
  frame_interf=v_rfft(v_enframe(interference_signal(i,:),window,inc),N_window,2); 
  frame_noise=v_rfft(v_enframe(noise(i,:),window,inc),N_window,2); 

%   [frame_desired,T,WS]=v_enframe(desired_signal(i,:),window,1/overlap_factor,'fa'); 
%   [frame_interf,T,WS]=v_enframe(interference_signal(i,:),window,1/overlap_factor,'fa'); 
%   [frame_noise,T,WS]=v_enframe(noise(i,:),window,1/overlap_factor,'fa'); 
  for n=1:n_frame
      % mvdr
%       beam_frame_desired(n,:,i)=w_mvdr(i,:).*(frame_desired(n,:));
%       beam_frame_interf(n,:,i)=w_mvdr(i,:).*(frame_interf(n,:));
%       beam_frame_noise(n,:,i)=w_mvdr(i,:).*(frame_noise(n,:));

      % lcmv
      beam_frame_desired(n,:,i)=w_lcmv(i,:).*(frame_desired(n,:));
      beam_frame_interf(n,:,i)=w_lcmv(i,:).*(frame_interf(n,:));
      beam_frame_noise(n,:,i)=w_lcmv(i,:).*(frame_noise(n,:));

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
soundsc(output_signal,fs)

%% evaluation
%% calculate SIR improvement
sir_in=sir_test(desired_signal(1,:),interference_signal(1,:),fs);
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
score1=pesq('desired_signal.wav','received_signal.wav');
score2=pesq('desired_signal.wav','output_signal.wav');
del_score=score2-score1

% %% plot
% plot waveform
figure;
subplot(4,1,1)
plot(source);
title('source signal')
subplot(4,1,2)
plot(desired_signal(1,:));
title('desired signal in Mic 1')
subplot(4,1,3)
plot(interference_signal(1,:));
title('interference signal in Mic 1')
subplot(4,1,4);
plot(desired_signal(1,:));
hold on
plot(output_signal)
title('Beamformed signal (LCMV)')
legend('desired signal in mics','beamformed signal')

% spectrogram
figure;
subplot(1,3,1)
[s,f,t] = spectrogram(interference_signal(1,:),window,inc,N_window,fs) ;
spectrogram(desired_signal(1,:),window,inc,N_window,fs,'yaxis');
title('received signal in mic 1')
lim=caxis;
subplot(1,3,2)
spectrogram(interference_signal(1,:),window,inc,N_window,fs,'yaxis');
title('total received signal in mic 1')
caxis(lim);
subplot(1,3,3)
spectrogram(output_interf,window,inc,N_window,fs,'yaxis') ;
title('beamformed signal')
caxis(lim);


% plot the 3D beam pattern
theta = 180 * [-1:0.001:1];
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

figure;
k=1335;
plot(freqVec,P_lcmv(:,k))
hold on
k1=1112;
plot(freqVec,P_lcmv(:,k1))
title('Beam pattern of LCMV')
xlabel('Frequency')
ylabel('dB')
grid on

% 2D plot
figure;
subplot(1,2,1)
P_mvdr=20*log10(B_mvdr./max(B_mvdr,[],2));
imagesc(theta',f',P_mvdr);
set(gca,'YDir','normal')
title('MVDR beam pattern')
xlabel('Angle')
ylabel('Frequency')
colorbar
lim2=caxis;
subplot(1,2,2)
P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));
imagesc(theta',f',P_lcmv);
set(gca,'YDir','normal')
title('LCMV beam pattern')
xlabel('Angle')
ylabel('Frequency')
colorbar
% caxis(lim2);