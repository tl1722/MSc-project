% implement Data-dependent LCMV beamformer

% input: 
%   received desired signal, interf signal, noise signal: Nmics * length
%   cfg: parameters
% output:
%   desired output signal, interf output signal, noise output signal:
%   length * 1 vector
%   beamformer weight: (Nmics * N_fft) matrix
% TONGKAI LI, 2023.08

function [output_desired, output_interf,output_noise, total_w]=DD_LCMV(desired_signal,interference_signal,noise,cfg)

% frequency domain frame-based processing for each microphone
fs=cfg.fs;
overlap_factor = 2;   
% frame increment in samples
inc =512;                               
% DFT window length
N_window = inc*overlap_factor;                        
window = sqrt(hamming(N_window,'periodic'));

received_signal=desired_signal + interference_signal + noise;
% segment signal into frames and do DFT
for i=1:cfg.Nmics
  frame=v_rfft(v_enframe(received_signal(i,:),window,inc),N_window,2); 
  frame_desired=v_rfft(v_enframe(desired_signal(i,:),window,inc),N_window,2); 
  frame_interf=v_rfft(v_enframe(interference_signal(i,:),window,inc),N_window,2); 
  frame_noise=v_rfft(v_enframe(noise(i,:),window,inc),N_window,2); 

  total_frame(i,:,:)=frame;
  total_frame_desired(i,:,:)=frame_desired;
  total_frame_interf(i,:,:)=frame_interf;
  total_frame_noise(i,:,:)=frame_noise;

  n_frame=size(frame_desired,1);
  cfg.n_fft=size(frame_desired,2);
  n_fft=cfg.n_fft;
end

% calculate beamformer weight
rir_desired = rir_generator(cfg.c, cfg.fs, cfg.position_mics, cfg.position_source, cfg.room, 0, cfg.n);
rir_interf = rir_generator(cfg.c, cfg.fs, cfg.position_mics, cfg.position_interf, cfg.room, 0, cfg.n);
temp=fft(rir_desired,1024,2);
% steering vector
vs=temp(:,1:n_fft);
vs=vs./vs(1,:);
temp=fft(rir_interf,1024,2);
vn=temp(:,1:n_fft);
vn=vn./vn(1,:);
% steering vector
% vs = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_source*pi/180)/c); % signal direction vector
% vn = exp(-1i*2*pi*[0:num_mics-1].'*freqVec*d*cos(theta_interf*pi/180)/c); % interference direction vector



freqVec=linspace(0, fs/2, n_fft);
% time constant of the recursive estimator
alpha=0.98;
for frm=1:n_frame
  w_lcmv = zeros(cfg.Nmics, length(freqVec));
  for f = 1:length(freqVec)
     % estimate of covariance matrix
     if frm == 1
        R = total_frame(:,frm,f) * total_frame(:,frm,f)';
     else
        R = alpha * R + (1-alpha) * total_frame(:,frm,f) * total_frame(:,frm,f)';
     end
     eps=1e-4;
     if rcond(R)<100
       R = R+eye(cfg.Nmics)*min(diag(R))*eps;
     end
     % LCMV beamformer weight
     C = [vs(:, f),vn(:,f)];
     denominator = C'/R*C;
     if rcond(denominator)<100
      denominator = denominator+eye(size(C,2))*min(diag(denominator))*eps;
     end
     F=[1,0].';
     w_lcmv(:, f) = (R\C) / (denominator)*F;
     total_w(:,f,frm)=w_lcmv(:, f);
     
     % apply the beamformer on the speech signals
     beam_frame_desired(frm,f)=w_lcmv(:,f)'*total_frame_desired(:,frm,f);
     beam_frame_interf(frm,f)=w_lcmv(:,f)'*total_frame_interf(:,frm,f);
     beam_frame_noise(frm,f)=w_lcmv(:,f)'*total_frame_noise(:,frm,f);

  end

end

total_w=mean(total_w,3);
% fir_coef=circshift(v_irfft(total_w,N_window,2),inc,2);
% figure;
% plot(fir_coef(1,:));

% do ifft and overlap and add
output_desired=v_overlapadd(v_irfft(beam_frame_desired,N_window,2),window,inc); 
output_interf=v_overlapadd(v_irfft(beam_frame_interf,N_window,2),window,inc); 
output_noise=v_overlapadd(v_irfft(beam_frame_noise,N_window,2),window,inc); 

end