%% draft: adaptive data independent LCMV

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'
addpath 'speech\male'
addpath 'speech\female'
addpath '0714new\'
folder1 = 'F:\恺\project\code\speech\male';


[audio,fs]=readData(folder1);

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
% TDOA = sqrt(sum((bsxfun(@minus, position_source, position_mics)).^2, 2))/c; 

% interference location 
D_interf=2; 
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
rev_time = 0.3;             % Reverberation time (s)
n = 4096;                   % Number of samples

rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, n);
rir_interf = rir_generator(c, fs, position_mics, position_interf, room, rev_time, n);

%% calculate R using data independent model
overlap_factor = 2;   
% frame increment in samples
inc =512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

R=zeros(num_mics,num_mics);
theta_beam=0:30:180; 
theta_null=0:30:180;
w_all=cell(length(theta_beam),length(theta_null));
d_vs=cell(length(theta_beam),length(theta_null));
d_vn=cell(length(theta_beam),length(theta_null));

for bm=1:length(theta_beam)
    for nl=1:length(theta_null)
        position_beam = D*[cos(theta_beam(bm)/180*pi) sin(theta_beam(bm)/180*pi) 0] + position_mics(1,:); 
        position_null = D_interf*[cos(theta_null(nl)/180*pi) sin(theta_null(nl)/180*pi) 0] + position_mics(1,:); 

        h_beam = rir_generator(c, fs, position_mics, position_beam, room, 0, n);
        h_null = rir_generator(c, fs, position_mics, position_null, room, 0, n);
        temp=fft(h_beam,N_window,2);
        vs=temp(:,1:513);
        vs=vs./vs(1,:);
        temp=fft(h_null,N_window,2);
        vn=temp(:,1:513);
        vn=vn./vn(1,:);


freqVec=linspace(0, fs/2, 513);
alpha=0.98;
eps=1e-4;

% steering vector
% vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source*pi/180)/c); % signal direction vector
% vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interf*pi/180)/c); % interference direction vector

  % compute MVDR/LCMV weight
  w_lcmv = zeros(num_mics, length(freqVec));
  for f = 1:length(freqVec)
     R = sinc(2*(f-1)/N_window*fs*l_mn/c);
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
  end
  w_all{bm,nl}=w_lcmv;
  bm
  nl
    end
end


%% combine the rir to the speech signal
%load the audio signal
for k=1:1
% source
source=audio{k};
L=length(source);
% interference
SIR=0;
rng(100);
[interf1,fs2]=readwav('ieee58f03.wav');
[P,Q] = rat(fs/fs2);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);
interf1 = sqrt( sum(sum(source.^2)) / sum(sum(interf1.^2))/(10^(SIR/10)) ) * interf1;

% get the received signal x in each microphone
desired_signal=zeros(num_mics,length(source));
interference_signal=zeros(num_mics,length(source));

for i=1:num_mics
  desired_signal(i,:)=filter(rir_desired(i,:),1,source)*10;
  interference_signal(i,:)=filter(rir_interf(i,:),1,interf1)*10;
end

noise=0*randn(num_mics,L);

received_signal = desired_signal + interference_signal+ noise;
% soundsc(received_signal(1,:),fs)

% frequency domain frame-based processing for each microphone

% segment signal into frames and do DFT
n_fft=N_window/2+1;
n_frame=floor(L/N_window*2) - 1;

total_frame=zeros(num_mics, n_frame, n_fft);
total_frame_desired=zeros(num_mics, n_frame, n_fft);
total_frame_interf=zeros(num_mics, n_frame, n_fft);
total_frame_noise=zeros(num_mics, n_frame, n_fft);
beam_frame_desired=zeros(n_frame,n_fft);
beam_frame_interf=zeros(n_frame,n_fft);
beam_frame_noise=zeros(n_frame,n_fft);

for i=1:num_mics
  [frame_total,T,WS]=v_enframe(received_signal(i,:),window,1/overlap_factor,'fa'); 
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


% apply the beamformer
for bm=1:length(theta_beam)
    for nl=1:length(theta_null)
     w_lcmv=w_all{bm,nl};
     for frm=1:n_frame
        for f=1:length(freqVec)
        beam_frame_desired(frm,f)=w_lcmv(:,f)'*total_frame_desired(:,frm,f);
        beam_frame_interf(frm,f)=w_lcmv(:,f)'*total_frame_interf(:,frm,f);
        beam_frame_noise(frm,f)=w_lcmv(:,f)'*total_frame_noise(:,frm,f);

        end
      end
% output=v_overlapadd(v_irfft(beam_frame,N_window,2),window,inc); 

output_desired=v_overlapadd(v_irfft(beam_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(beam_frame_interf,N_window,2),window,inc); 
output_noise=v_overlapadd(v_irfft(beam_frame_noise,N_window,2),window,inc); 

output_signal= output_desired + output_interf + output_noise;

% soundsc(output_signal,fs)
% soundsc(received_signal(1,:),fs)

% calculate SIR improvement
sir_in=sir_test(desired_signal(1,:),interference_signal(1,:),fs);
sir_out=sir_test(output_desired,output_interf,fs);

sir_imp=sir_out-sir_in;
SIR_all(bm,nl)=sir_imp;

score1=pesq_test(source, received_signal(1,:), fs);
score2=pesq_test(source, output_signal, fs);
SCORE(bm,nl)=score2 - score1;

DRR_imp(bm,nl)=DRR_test(rir_desired, fs, w_lcmv, N_window);

bm
nl

    end
end
SCORE_final(k,1)=max(max(SCORE));
[x,y]=find(SCORE==SCORE_final(k));
SIR_final(k,1)=SIR_all(x,y);
DRR_final(k,1)=DRR_imp(x,y);
k
DRR_final(k,1)
end

% % %% plot
% % plot waveform
% figure;
% subplot(4,1,1)
% plot(source);
% title('source signal')
% subplot(4,1,2)
% plot(desired_signal(1,:));
% title('desired signal in Mic 1')
% subplot(4,1,3)
% plot(interference_signal(1,:));
% title('interference signal in Mic 1')
% subplot(4,1,4);
% plot(desired_signal(1,:));
% hold on
% plot(output_signal)
% title('Beamformed signal (LCMV)')
% legend('desired signal in mics','beamformed signal')
% 
% % spectrogram
% figure;
% subplot(1,3,1)
% [s,f,t] = spectrogram(interference_signal(1,:),window,inc,N_window,fs) ;
% spectrogram(desired_signal(1,:),window,inc,N_window,fs,'yaxis');
% title('received signal in mic 1')
% lim=caxis;
% subplot(1,3,2)
% spectrogram(interference_signal(1,:),window,inc,N_window,fs,'yaxis');
% title('total received signal in mic 1')
% caxis(lim);
% subplot(1,3,3)
% spectrogram(output_interf,window,inc,N_window,fs,'yaxis') ;
% title('beamformed signal')
% caxis(lim);
% 
% 
% % plot the 3D beam pattern
% theta = 180 * [-1:0.01:1];
% % B_mvdr=zeros(n_fft,length(theta));
% B_lcmv=zeros(n_fft,length(theta));
% for i=1:length(f)
%   v = exp(1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*f(i)/c);
% %   B_mvdr(i,:) = abs(w_mvdr(:,i)' * v);
%   B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);
% 
% end
% 
% 
% P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));
% imagesc(theta',f',P_lcmv);
% set(gca,'YDir','normal')
% title('LCMV beam pattern')
% xlabel('Angle')
% ylabel('Frequency')
% colorbar
% 
% 
% 
% subplot(1,2,2)
% mesh(f,theta,P_lcmv.');
% set(gca,YDir="reverse")
% title('LCMV beam pattern')
% xlabel('Frequency')
% ylabel('angle')
% zlabel('magnitude(dB)')
% 
% % figure;
% % k=70;
% % plot(theta,20*log10(B_mvdr(k,:)/max(B_mvdr(k,:))))
% % hold on
% % title('Beam pattern of LCMV')
% % xlabel('degree')
% % ylabel('dB')
% % grid on
