%% ex4: MVDR beamforming vs LCMV (far field)

clc;
clear;

num_mics=8;
c=340;
f=300;
lambda = c/f;
d =0.3;
thetas = [20];    % desired angle
thetai = [30 60]; % interference angle
n = [0:num_mics-1]';
vs = exp(-1i*2*pi*n*d*sind(thetas)/lambda); % signal direction vector
vn = exp(-1i*2*pi*n*d*sind(thetai)/lambda); % interference direction vector

 
[di,fs]=readwav('S_01_01.wav');    % desired signal
di=di.';
L=length(di);
t = [0:L-1];

vn1 = 3*randn(length(di),1).';  % interference signal 1
vn2 = 3*randn(length(di),1).';  % interference signal 2

% [source_signal2,fs2]=readwav('S_01_02.wav');
% l1=length(di);
% l2=length(source_signal2);
% vn1=[source_signal2;zeros(l1-l2,1)].';
% 
% [source_signal3,fs3]=readwav('S_01_03.wav');
% l3=length(source_signal3);
% vn2=[source_signal3;zeros(l1-l3,1)].';


A = [vs vn];
St = [di;vn1;vn2];
Xt = A*St + randn(num_mics,L);   % the received signal in microphones

% MVDR
R_x = 1/L*(Xt*Xt');
R_x_inv = inv(R_x);
W_MVDR = R_x_inv*vs/(vs'*R_x_inv*vs);

% LCMV
C=[vs,vn];
f=[1,0,0].';
R_x_inv = inv(R_x);
W_LCMV = R_x_inv*C/(C'*R_x_inv*C)*f;

% plot the beam pattern
theta = 180 * [-1:0.01:1];
v = exp(-1i*2*pi*n* d*sind(theta)/lambda);
B_MVDR = abs(W_MVDR' * v);
plot(theta,20*log10(B_MVDR/max(B_MVDR)))
hold on
v = exp(-1i*2*pi*n* d*sind(theta)/lambda);
B_LCMV = abs(W_LCMV' * v);
plot(theta,20*log10(B_LCMV/max(B_LCMV)))
title('Beam pattern of MVDR and LCMV')
legend('MVDR','LCMV')
xlabel('degree')
ylabel('dB')
grid on

% wideband LCMV
n_fft=160;
f=linspace(0, fs/2, n_fft);
W_opt=zeros(num_mics,n_fft);
for i=1:length(f)
  vs = exp(-1i*2*pi*[0:num_mics-1].'*d*sind(thetas)*f(i)/c); % signal direction vector
  vn1 = exp(-1i*2*pi*[0:num_mics-1].'*d*sind(thetai(1))*f(i)/c); % interference direction vector
  vn2 = exp(-1i*2*pi*[0:num_mics-1].'*d*sind(thetai(2))*f(i)/c); % interference direction vector
  C=[vs,vn1,vn2];
  F=[1,0,0].';
  W_opt(:,i) = (R_x_inv*C) / (C'*R_x_inv*C)*F;
end

% plot the beam pattern
theta = 180 * [-1:0.01:1];
B=zeros(n_fft,length(theta));
figure;
for i=1:length(f)
  v = exp(-1i*2*pi*[0:num_mics-1].'* d*sind(theta)*f(i)/c);
  B(i,:) = abs(W_opt(:,i)' * v);
  plot(theta,20*log10(B(i,:)/max(B(i,:))))
  hold on
end
title('Beam pattern of LCMV')
xlabel('degree')
ylabel('dB')
grid on

figure;
P=20*log10(B);
mesh(f,theta,P.');

