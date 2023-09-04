%% ex9: MOVEMENT

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% load the audio signal
% source
[source,fs]=readwav('S_01_01.wav');
L=length(source);

% stft
overlap_factor = 2;   
% frame increment in samples
inc = 256;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = hamming(N_window,'periodic');
[frame_source,T,WS]=v_enframe(source,window,1/overlap_factor,'fa'); 
n_frame=size(frame_source,1);
n_fft=size(frame_source,2);

% interference
rng(100);
% interf1=0.01*randn(1,L);
[interf1,fs1]=readwav('babble noise.wav');
[P,Q] = rat(fs/fs1);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);

[frame_interf,T,WS]=v_enframe(interf1,window,1/overlap_factor,'fa'); 


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

% distance matrix of mics
l_mn = zeros(num_mics, num_mics);
for i = 1:num_mics
    for j = 1:num_mics
        l_mn(i, j) = abs(i - j) * d;
    end
end


%% RIR parameters
room = [5 4 6];                % Room dimensions [x y z] (m)
rev_time = 0.3;             % Reverberation time (s)
m = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order =-1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

% source location 
D=1; 
theta_source = 60;   % DOA 

D_interf=1; 
theta_interf = 20;   % DOA 
position_interf = D_interf*[cos(theta_interf/180*pi) sin(theta_interf/180*pi) 0] + position_mics(1,:); 
rir_undesired = rir_generator(c, fs, position_mics, position_interf, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
atf=fft(rir_undesired,n_fft,2);

%% pre-calculated covariance matrix
%% calculate R
theta_noise=(0:10:360)'; 
Rnn=zeros(num_mics,num_mics);
for k=1:length(theta_noise)
  position_noise = D*[cos(theta_noise(k)/180*pi) sin(theta_noise(k)/180*pi) 0] + position_mics(1,:);   
  order=0;
  rir_diffuse = rir_generator(c, fs, position_mics, position_noise, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
  Rnn=Rnn+rir_diffuse*rir_diffuse';
end
Rnn=Rnn/length(theta_noise);
freqVec=linspace(0, fs/2, n_fft);

%% start processing
desired_signal=zeros(n_frame,n_fft,num_mics);
interf_signal=zeros(n_frame,n_fft,num_mics);

frame_s1=floor((L1+sample_internal)/inc);
for n=1:n_frame
    % source DOA
    if n<frame_s1
       theta_source=60;
    else
       theta_source=120;
    end
   position_source = D*[cos(theta_source/180*pi) sin(theta_source/180*pi) 0] + position_mics(1,:); 
   rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
   atf=fft(rir_desired,n_fft,2);

  %% combine the rir to the speech signal
  for i=1:num_mics
    desired_signal(n,:,i)=atf(i,:).*frame_source(n,:);
    interf_signal(n,:,i)=atf_i(i,:).*frame_interf(n,:);
  end
end
received_signal=desired_signal+interf_signal;

% create a moving signal
received_frame=sum(received_signal,3);
received=v_overlapadd(v_irfft(received_frame,N_window,2),window,inc); 
% soundsc(received,fs);

%% apply beamforming
beam_frame=zeros(n_frame,n_fft,num_mics);
for n=1:n_frame
  theta_interf=theta_interf+0.1*n;
  % steering vector
  vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source*pi/180)/c); % signal direction vector
  vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interf*pi/180)/c); % interference direction vector

  % compute MVDR/LCMV weight
  w_mvdr = zeros(num_mics, length(freqVec));
  w_lcmv = zeros(num_mics, length(freqVec));
  eps=1e-4;
  for f = 1:length(freqVec)
     % LCMV
     C = [vs(:, f),vn(:,f)];
     denominator = C'/Rnn*C;
     if rcond(denominator)<100
       denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
     end
     F=[1,0].';
     w_lcmv(:, f) = (Rnn\C) / (denominator)*F;
   
     % MVDR
     vf=vs(:, f);
     w_mvdr(:, f) = (Rnn\vf) / (vf'/Rnn*vf);
   end

  % compute output frame
  for i=1:num_mics
     beam_frame(n,:,i)=w_lcmv(i,:).*(received_signal(n,:,i));
  end

    % plot the 3D beam pattern
  theta = 180 * [-1:0.01:1];
  B_lcmv=zeros(n_fft,length(theta));
  for i=1:length(freqVec)
    v = exp(-1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*freqVec(i)/c);
    B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);
  end

  P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));


  clf;
  imagesc(theta',freqVec,P_lcmv);
  hold on
  xline(theta_source,'blue',LineWidth=2);
  xline(theta_interf,'red',LineWidth=2);
  title(['LCMV beam pattern, desired signal at ',num2str(theta_source),' deg and interf signal at ',num2str(theta_interf),' deg.']);
  set(gca,'YDir','normal')
  xlabel('Angle')
  ylabel('Frequency')
  colorbar
  caxis([-70 0]);

  pause(0.1)

end

% sum Y in every mics
sum_frame=sum(beam_frame,3);

% use overlapadd to reconstruct the signal y
output_signal=v_overlapadd(v_irfft(sum_frame,N_window,2),window,inc); 

soundsc(output_signal,fs)
