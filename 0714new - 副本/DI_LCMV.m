% implement Data-independent LCMV beamformer

% input: 
%   received desired signal, interf signal, noise signal: Nmics * length
%   cfg: parameters
% output:
%   desired output signal, interf output signal, noise output signal:
%   length * 1 vector
%   beamformer weight: (Nmics * N_fft) matrix
% TONGKAI LI, 2023.08


function [output_desired, output_interf,output_noise, total_w]=DI_LCMV(desired_signal,interference_signal,noise,cfg)
% frequency domain frame-based processing for each microphone
fs=cfg.fs;
num_mics=cfg.Nmics;
overlap_factor = 2;   
% frame increment in samples
inc =512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

received_signal=desired_signal + interference_signal +noise;

% segment signal into frames and do FFT
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
  cfg.n_fft=size(frame_total,2);
  n_fft=cfg.n_fft;
end

% calculate beamformer weight using data independent model
rir_desired = rir_generator(cfg.c, cfg.fs, cfg.position_mics, cfg.position_source, cfg.room, 0, cfg.n);
rir_interf = rir_generator(cfg.c, cfg.fs, cfg.position_mics, cfg.position_interf, cfg.room, 0, cfg.n);
temp=fft(rir_desired,N_window,2);
vs=temp(:,1:n_fft);
vs=vs./vs(1,:);
temp=fft(rir_interf,N_window,2);
vn=temp(:,1:n_fft);
vn=vn./vn(1,:);

freqVec=linspace(0, fs/2, n_fft);
eps=1e-4;

w_lcmv = zeros(num_mics, length(freqVec));

for f = 1:length(freqVec)
    % covariance matrix
     R = sinc(2*(f-1)/N_window*fs*cfg.l_mn/cfg.c);
     % diagnal loading
     if rcond(R)<100
       R = R+eye(num_mics)*min(diag(R))*eps;
     end

     % calculate LCMV beamformer
     C = [vs(:, f),vn(:,f)];
     denominator = C'/R*C;
     if rcond(denominator)<100
      denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
     end
     F=[1,0].';
     w_lcmv(:, f) = (R\C) / (denominator)*F;
end
total_w=w_lcmv;

% apply beamformer on received signals
for frm=1:n_frame
  for f = 1:length(freqVec)
     beam_frame_desired(frm,f)=w_lcmv(:,f)'*total_frame_desired(:,frm,f);
     beam_frame_interf(frm,f)=w_lcmv(:,f)'*total_frame_interf(:,frm,f);
     beam_frame_noise(frm,f)=w_lcmv(:,f)'*total_frame_noise(:,frm,f);

  end
end
% do IFFT and overlap and add
output_desired=v_overlapadd(v_irfft(beam_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(beam_frame_interf,N_window,2),window,inc); 
output_noise=v_overlapadd(v_irfft(beam_frame_noise,N_window,2),window,inc); 


% fir coefficient of LCMV (optional)
% total_w=mean(total_w,3);
% fir_coef=circshift(v_irfft(total_w,N_window,2),inc,2);
% figure;
% plot(fir_coef(1,:));
end