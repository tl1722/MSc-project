%% ex5: adaptive data independent model

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

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
interf1 = sqrt( sum(sum(source.^2)) / sum(sum(interf1.^2))/(10^(SIR/10)) ) * interf1;

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
theta_source = 60;   % DOA 
position_source = D*[cos(theta_source/180*pi) sin(theta_source/180*pi) 0] + position_mics(1,:); 
TDOA = sqrt(sum((bsxfun(@minus, position_source, position_mics)).^2, 2))/c; 

% interference location 
D_interf=3; 
theta_interf=20; 
position_interf = D_interf*[cos(theta_interf/180*pi) sin(theta_interf/180*pi) 0] + position_mics(1,:); 

% distance matrix of mics
l_mn = zeros(num_mics, num_mics);
for i = 1:num_mics
    for j = 1:num_mics
        l_mn(i, j) = abs(i - j) * d;
%         l_mn(i,j)=norm(position_mics(i,:)-position_mics(j,:))

    end
end


%% RIR generator
room = [5 4 6];                % Room dimensions [x y z] (m)
rev_time = 0;             % Reverberation time (s)
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

noise=0.001*randn(num_mics,L);
% noise=awgn(desired_signal,30,'measured');
received_signal = desired_signal + interference_signal+ noise;
% soundsc(received_signal(1,:),fs)

%% frequency domain frame-based processing for each microphone
overlap_factor = 2;   
% frame increment in samples
inc = 256;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = hamming(N_window,'periodic');

% segment signal into frames and do DFT
for i=1:num_mics
  [frame_interf,T,WS]=v_enframe(interference_signal(i,:),window,1/overlap_factor,'fa'); 
  [frame_noise,T,WS]=v_enframe(noise(i,:),window,1/overlap_factor,'fa'); 
  total_frame_interf(i,:,:)=frame_interf;
  total_frame_noise(i,:,:)=frame_noise;
  n_frame=size(frame_interf,1);
  n_fft=size(frame_interf,2);
end

%% calculate w using data independent model
R=zeros(num_mics,num_mics);
theta_beam=0:10:180;
theta_null=0:20:180;

freqVec=linspace(0, fs/2, n_fft);

% compute MVDR/LCMV weight
w_mvdr = zeros(num_mics, length(freqVec));
w_lcmv = zeros(num_mics, length(freqVec));
total_w=cell(length(theta_beam),length(theta_null));
eps=1e-4;

for bm=1:length(theta_beam)
    for nl=1:length(theta_null)
      for f = 1:length(freqVec)
        R=sinc(2*(f-1)/N_window *fs*l_mn/c);
        if rcond(R)<100
           R = R+eye(num_mics)*min(diag(R))*eps;
        end
      % steering vector
        vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_beam(bm)*pi/180)/c); % signal direction vector
        vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_null(nl)*pi/180)/c); % interference direction vector

      %   Rnn=sin(2*pi*f*l_mn/c)/(2*pi*f*l_mn/c);

         % LCMV
         C = [vs(:, f),vn(:,f)];
         denominator = C'/R*C;
         if rcond(denominator)<100
           denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
         end
         F=[1,0].';
         w_lcmv(:, f) = (R\C) / (denominator)*F;  
      end
      total_w{bm,nl}=w_lcmv;

      % compute output frame
      for i=1:num_mics
        [frame_desired,T,WS]=v_enframe(desired_signal(i,:),window,1/overlap_factor,'fa'); 
        [frame_interf,T,WS]=v_enframe(interference_signal(i,:),window,1/overlap_factor,'fa'); 
        [frame_noise,T,WS]=v_enframe(noise(i,:),window,1/overlap_factor,'fa');   
        for n=1:n_frame

            % lcmv
            beam_frame1(n,:,i)=w_lcmv(i,:).*(frame_desired(n,:));
            beam_frame2(n,:,i)=w_lcmv(i,:).*(frame_interf(n,:));
            beam_frame3(n,:,i)=w_lcmv(i,:).*(frame_noise(n,:));

        end
      end

% sum Y in every mics
      sum_frame1=sum(beam_frame1,3);
      sum_frame2=sum(beam_frame2,3);
      sum_frame3=sum(beam_frame3,3);

% use overlapadd to reconstruct the signal y
      output_desired=v_overlapadd(v_irfft(sum_frame1,N_window,2),window,inc); 
      output_interf=v_overlapadd(v_irfft(sum_frame2,N_window,2),window,inc); 
      output_noise=v_overlapadd(v_irfft(sum_frame3,N_window,2),window,inc); 

     sir_out=sir_test(output_desired,output_interf,fs);
     sir(bm,nl)=sir_out;
     bm
     nl
    end
end


%% calculate SIR improvement
sir_in=sir_test(desired_signal(1,:),interference_signal(1,:),fs);


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
theta = 180 * [-1:0.01:1];
B_mvdr=zeros(n_fft,length(theta));
B_lcmv=zeros(n_fft,length(theta));
for i=1:length(f)
  v = exp(-1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
  B_mvdr(i,:) = abs(w_mvdr(:,i)' * v);
  B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);

end

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
caxis(lim2);

% 3D plot
figure;
subplot(1,2,1)
mesh(f,theta,P_mvdr.');
set(gca,YDir="reverse")
title('MVDR beam pattern')
xlabel('Frequency')
ylabel('angle')
zlabel('magnitude(dB)')

subplot(1,2,2)
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
