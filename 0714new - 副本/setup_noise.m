%% setup noise

%% interference
rng(100)
SIR=0;
[interf1,fs2]=readwav('ieee58f03.wav');
[P,Q] = rat(fs/fs2);
interf1 = resample(interf1,P,Q);
interf1=interf1(1:L);
% scale
interf1 = sqrt( sum(sum(source.^2)) / sum(sum(interf1.^2))/(10^(SIR/10)) ) * interf1;

% interference location 
D_interf=cfg.D_interf; 
theta_interf=cfg.theta_interf; 
cfg.position_interf = D_interf*[cos(theta_interf/180*pi) sin(theta_interf/180*pi) 0] + cfg.position_mics(1,:); 

rir_interf = rir_generator(cfg.c, cfg.fs, cfg.position_mics, cfg.position_interf, cfg.room, cfg.T60, cfg.n);

interf=zeros(cfg.Nmics,L);
for ch=1:cfg.Nmics
  interf(ch,:)=10*filter(rir_interf(ch,:),1,interf1);
end

%% babble noise
[speech,fs3]=readwav('babble noise.wav');
[P,Q] = rat(fs/fs3);
speech = resample(speech,P,Q);

SBNR=0;
theta_noise=[0,90,180,270];
D_noise=2.5;
bnoise_all=zeros(cfg.Nmics,L);
for p=1:4
   position_noise = D_noise*[cos(theta_noise(p)/180*pi) sin(theta_noise(p)/180*pi) 0] + cfg.position_mics(1,:); 
   rir_noise = rir_generator(cfg.c, cfg.fs, cfg.position_mics, position_noise, cfg.room, cfg.T60, cfg.n);
   bnoise=speech(p*100000+1:p*100000+L);
   bnoise =sqrt(sum(sum(source.^2)) / sum(sum(bnoise.^2))/(10^(SBNR/10)) ) * bnoise;
   noise_all=zeros(cfg.Nmics,L);
   for ch=1:cfg.Nmics
      noise_all(ch,:)=filter(rir_noise(ch,:),1,bnoise)*10;
   end
   bnoise_all=bnoise_all+noise_all;
end

%% sensor noise
snoise=0.0001*randn(cfg.Nmics,L);

% noise=snoise+bnoise_all;
noise=snoise;
