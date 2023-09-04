%% dynamic beampattern varying with time
Tongkai Li, 2023.08

clear
clc
addpath 'E:\MATLAB2021b\RIR-Generator-master'

%% parameters 
fs=22050;
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

% distance matrix of mics
l_mn = zeros(num_mics, num_mics);
for i = 1:num_mics
    for j = 1:num_mics
        l_mn(i, j) = abs(i - j) * d;
    end
end

% Source position [x y z] (m)
% v_x=0.5;
% v_y=0;
% v_z=0;
% t=[0:0.1:length(source)/fs]';

%% RIR parameters
room = [8 6 5];                % Room dimensions [x y z] (m)
rev_time = 0.3;             % Reverberation time (s)
m = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order =-1;                  % Reflection order
dim = 3;                    % Room dimension
orientation=0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter


% for n=1:length(t)
%   position_source = [2+v_x*t(n) 3.5+v_y*t(n) 2+v_z*t(n)]; 

  % interference position [x y z] (m) 
  % position_interf = [3 3 2];   

  % DOA 
  
  thetas_source=46.6; 
  thetas_interf=6.8;
  n_fft=201;
% % 
% %   rir_desired = rir_generator(c, fs, position_mics, position_source, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
% %   %   rir_interf = rir_generator(c, fs, position_mics, position_interf, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
% % 
% %   %% combine the rir to the speech signal
% %   % get the received signal x in each microphone
% %   desired_signal=zeros(num_mics,length(source)+1);
% %   interference_signal=zeros(num_mics,length(source));
% % 
% %   for i=1:num_mics
% %     desired_signal(i,((n-1)*0.1*fs+1):n*0.1*fs)=filter(rir_desired(i,:),1,source(((n-1)*0.1*fs+1):n*0.1*fs))*10;
% % %   interference_signal(i,:)=filter(rir_interf(i,:),1,interf1)*10;
% %   end
% % end
% % 
% % % noise=0.001*randn(num_mics,L);
% % noise=awgn(desired_signal,30,'measured');
% % received_signal = desired_signal + noise;
% % % soundsc(received_signal(1,:),fs)
% % 
% % %% frequency domain frame-based processing for each microphone
% % overlap_factor = 2;   
% % % frame increment in samples
% % inc = 256;                               
% % % DFT window length
% % N_window = inc*overlap_factor;                        
% % window = hamming(N_window,'periodic');
% % 
% % % segment signal into frames and do DFT
% % for i=1:num_mics
% %   [frame_interf,T,WS]=v_enframe(interference_signal(i,:),window,1/overlap_factor,'fa'); 
% %   [frame_noise,T,WS]=v_enframe(noise(i,:),window,1/overlap_factor,'fa'); 
% %   total_frame_interf(i,:,:)=frame_interf;
% %   total_frame_noise(i,:,:)=frame_noise;
% %   n_frame=size(frame_interf,1);
% %   n_fft=size(frame_interf,2);
% % end
% % 
% % %% calculate w using data independent model
% % D_noise=2; 
% % theta_noise=(0:10:180)'; 
% % Rnn=zeros(num_mics,num_mics);
% % % for k=1:length(theta_noise)
% % %   position_noise = D_noise*[cos(theta_noise(k)/180*pi) sin(theta_noise(k)/180*pi) 0] + position_mics(1,:);   
% % %   order=0;
% % %   rir_diffuse = rir_generator(c, fs, position_mics, position_noise, room, rev_time, n, mtype, order, dim, orientation, hp_filter);
% % %   Rnn=Rnn+rir_diffuse*rir_diffuse';
% % % end

% steering vector

for n=1:length(thetas_source)
  freqVec=linspace(0, fs/2, n_fft);
%   theta_noise=(0:10:180)'; 
%   Rnn=zeros(num_mics,num_mics);
%   D_noise=1;
%   for k=1:length(theta_noise)
%     position_noise = D_noise*[cos(theta_noise(k)/180*pi) sin(theta_noise(k)/180*pi) 0] + position_mics(1,:);   
%     order=0;
%     rir_diffuse = rir_generator(c, fs, position_mics, position_noise, room, rev_time, m, mtype, order, dim, orientation, hp_filter);
%     Rnn=Rnn+rir_diffuse*rir_diffuse';
%   end
%   Rnn=Rnn/length(theta_noise);
%     D=1;
%     position_source=D*[cos(thetas_source(n)/180*pi) sin(thetas_source(n)/180*pi) 0] + position_mics(1,:);
%     position_interf=D*[cos(thetas_interf(n)/180*pi) sin(thetas_interf(n)/180*pi) 0] + position_mics(1,:);
%     
%     d_desired = rir_generator(c, fs, position_mics, position_source, room, 0, m);
%     d_interf = rir_generator(c, fs, position_mics, position_interf, room, 0, m);
%     temp=fft(d_desired,400,2);
%     vs=temp(:,1:n_fft);
%     vs=vs./vs(1,:);
%     temp=fft(d_interf,400,2);
%     vn=temp(:,1:n_fft);
%     vn=vn./vn(1,:);


  vs = exp(1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(thetas_source(n)*pi/180)/c); % signal direction vector
  vn = exp(1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(thetas_interf(n)*pi/180)/c); % interference direction vector

  % compute MVDR/LCMV weight
  w_mvdr = zeros(num_mics, length(freqVec));
  w_lcmv = zeros(num_mics, length(freqVec));
  eps=1e-4;
  for f = 1:length(freqVec)
%      Rnn=sin(2*pi*f*l_mn/c/n_fft)/(2*pi*f*l_mn/c/n_fft);
     Rnn = sinc(2*(f-1)/400*fs*l_mn/c);
     if rcond(Rnn)<100
       Rnn = Rnn+eye(num_mics)*min(diag(Rnn))*eps;
     end
     % LCMV
     C = [vs(:, f),vn(:,f)];
     denominator = C'/Rnn*C;
     if rcond(denominator)<eps
       denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
     end
     F=[1,0].';
     w_lcmv(:, f) = (Rnn\C) / (denominator)*F;
   
     % MVDR
     vf=vs(:, f);
     w_mvdr(:, f) = (Rnn\vf) / (vf'/Rnn*vf);
  end

  % plot the 3D beam pattern
  theta = 180 * [-1:0.01:1];
  B_mvdr=zeros(n_fft,length(theta));
  B_lcmv=zeros(n_fft,length(theta));
  for i=1:length(freqVec)
    v = exp(1i*2*pi*[0:num_mics-1].'* d*cosd(theta)*freqVec(i)/c);
    B_mvdr(i,:) = abs(w_mvdr(:,i)' * v);
    B_lcmv(i,:) = abs(w_lcmv(:,i)' * v);
  end

 P_mvdr=20*log10(B_mvdr./max(B_mvdr,[],2));
 P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));


clf;
imagesc(theta',freqVec,P_lcmv);
hold on
xline(thetas_source(n),'blue',LineWidth=2);
xline(thetas_interf(n),'red',LineWidth=2);
title(['Target at ',num2str(thetas_source(n)),' deg and interferer at ',num2str(thetas_interf(n)),' deg.'],'FontSize',17,'FontName','Times New Roman');
set(gca,'YDir','normal')
xlabel('DOA (Degree)','FontSize',15,'FontName','Times New Roman')
ylabel('Frequency (Hz)','FontSize',15,'FontName','Times New Roman')
colorbar
caxis([-70 0]);

pause(10)

% clf;
% k=100;
% plot(theta,20*log10(B_mvdr(k,:)/max(B_mvdr(k,:))))
% hold on
% grid on
% plot(theta,20*log10(B_lcmv(k,:)/max(B_lcmv(k,:))))
% xlabel('degree')
% ylabel('dB')
% xline(thetas_source(n),'blue');
% xline(thetas_interf(n));
% 
% pause(0.1)
end
