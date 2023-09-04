%% ex1: beam pattern
clear; 
clc;

set(0,'DefaultAxesFontSize',11)
set(0,'defaultAxesFontName','Times New Roman')
set(0,'DefaultLineLineWidth',1);

%% parameters
c=340;
fs=22050;                  % sampling frequency
f0=[500,1000,5000,8000];   % signal frequency
f0=1:30:11000;
lambda=c./f0;              % wave length
d=0.04;                   % mic distance 
num_mics=8;                       % mic number

%% beamformer
theta_desired=60*pi/180;   % desired direction
for i=1:length(f0)
    w(i,:)=exp(-1i*2*pi*f0(i)*([0:num_mics-1]*d*cos(theta_desired)/c))/num_mics;
end

%% compute beam pattern
theta_angle=[-180:0.1:180];
theta=theta_angle*pi/180;
beam = zeros(size(theta));
 
for i=1:length(f0)
  for j=1:length(theta)
    a=exp(-1i*2*pi*f0(i)*([0:num_mics-1]*d*cos(theta(j)))/c);   % direction vector
    beam(i,j)=w(i,:)*a';
  end
end
Beam_db =20*log10(abs(beam));
limit_dB = -50;
index = Beam_db <limit_dB;
Beam_db(index)=limit_dB;
 
% polar response
% figure;
% for i=1:length(f0)
%   subplot(2,2,i)
%   polarplot(theta, Beam_db(i,:),'LineWidth',1.2);
%   rlim([-50 0])
%   thetalim([-180 180])
%   title(['f = ' num2str(f0(i)) ' Hz']);
% end
% grid on;
% legend('f=500Hz','f=1000Hz','f=5000Hz','f=8000Hz')


% standard beam pattern
figure;
for i=1:length(f0)
  plot(theta_angle, Beam_db(i,:),'LineWidth',1.2);
end
grid on;
title('beam pattern');
xlim([0 360])
% legend('f=500Hz','f=1000Hz','f=5000Hz','f=8000Hz')

figure;
mesh(theta_angle,f0,Beam_db);
set(gca,'YDir','normal')
title('3-Dimension DSB beam pattern','FontSize',14,'FontName','Times New Roman')
xlabel('DOA (Degrees)','FontSize',12,'FontName','Times New Roman')
ylabel('Frequency (Hz)','FontSize',12,'FontName','Times New Roman')
zlabel('Gain (dB)','FontSize',12,'FontName','Times New Roman')

