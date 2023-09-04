% LCMV beamformer on speech signals using frame-based signal processing
% author: Tongkai Li, MSc project, supervised by Prof. Patrick Naylor
% 2023.08

%% setup paths
clc
clear

addpath 'E:\MATLAB2021b\RIR-Generator-master'
addpath 'F:\恺\project\code'

%% load IEEE speech signals
folder1 = 'F:\恺\project\code\speech\male';
folder2 = 'F:\恺\project\code\speech\female';
[audio,fs]=readData(folder1);
disp('Data loading...')

%% setup parameters
cfg = [];
cfg.Nmics = 8;               % micphone number 
cfg.d = 0.05;                % distance between microphone (linear array) 
cfg.T60 = 0.3;               % Reverberation time (s)
cfg.D = 1;                   % distance between the desired source and array
cfg.theta_source = 60;       % DOA of the desired source (in degree) 
cfg.D_interf = 2;            % distance between the undesired source and array              
cfg.theta_interf = 20;       % DOA of the undesired source (in degree) 
cfg.k0 = 100;                % the number of signals for processing
setup_mics
setup_room

SIR_imp=zeros(cfg.k0,1);
sinr_imp=zeros(cfg.k0,1);
SCORE=zeros(cfg.k0,1);
score1=zeros(cfg.k0,1);
DRR_imp=zeros(cfg.k0,1);

%% load speech signal in turn
disp('Starting processing...')
for k=1:cfg.k0
   source=audio{k};
   L=length(source);
   desired = zeros(cfg.Nmics, L);
   % generate echoic desired signal
   for ch=1:cfg.Nmics
      desired(ch,:) = 10*filter(rir_desired(ch,:),1,source);
   end
   setup_noise
   % generate total received signal in microphone array
   received = desired + interf + noise;
   
   %% apply LCMV beamformer algorithm
   % DD means Data-dependent approach, DI means Data-independent approach,
   % select the one you would like to apply

%    [output_desired,output_interf,output_noise, total_w]=DD_LCMV(desired,interf,noise,cfg);
    [output_desired,output_interf,output_noise, total_w]=DI_LCMV(desired,interf,noise,cfg);
  
    % total output beamformed signal
   output_signal= output_desired + output_interf + output_noise;

   %% evaluation
   % sir
   sir_in=sir_test(desired(1,:),interf(1,:),fs);
   sir_out=sir_test(output_desired,output_interf,fs);
   SIR_imp(k)=sir_out-sir_in;

   % sinr
   sinr_in=sinr_test(desired(1,:),interf(1,:),noise(1,:),fs);
   sinr_out=sinr_test(output_desired,output_interf,output_noise,fs);
   sinr_imp(k)=sinr_out-sinr_in;

   % pesq improvement
   score1(k)=pesq_test(source, received(1,:), fs);
   score2(k)=pesq_test(source, output_signal, fs);
   SCORE(k)=score2(k) - score1(k);  

   % drr improvement
   DRR_imp(k)=DRR_test(rir_desired, fs, total_w, 1024);
   X = sprintf('Processing No. %d signal.',k);
   disp(X)
end

%% plot beam pattern and spectrograms
plot_pattern;
plot_spectrogram;

disp('Completed')


%% listen to the output signal (optional)
% soundsc(output_signal,fs)
% soundsc(received(1,:),fs)