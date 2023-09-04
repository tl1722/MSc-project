%%------------------------------------------------------------------------
% [1] R. Mars, V. G. Reju, and A. W. H Khong. "A Frequency-Invariant Fixed
% Beamformer for Speech Enhancement.", APSIPA-2014
%-------------------------------------------------------------------------
clc;
clearvars;

%% Simulation setup:
micNum          = 8;            % microphone number
d               = 0.04;        % mic spacing
noTap           = 200;           % filter taps
c               = 342;           % speed of sound                               
fs              = 16000;         % sampling rate
fl              = 200;           % lower freq limit
fh              = 3400;          % upper freq limit
f_step          = 50;            % freq step
fr              = 1700;          % reference freq
lkDir           = 90;            % target direction
inter1          = 20;            % interference 1 direction
inter2          = 150;           % interference 2 direction

%%
k               = [inter1,inter2]; 
% k               = [inter1]; 

x_array         = ones(1,micNum); 

%% Design Frequency invariant beamformer

disp('Processing: FIB design for target direction ...')
w_FIB = FIB(noTap,micNum,lkDir,fl,fh,f_step,fr,fs,x_array,d,c);

%% Design Frequency invariant null design
disp('Processing: FIB null design for interference directions ...')
weight_matrix_null = zeros(noTap*micNum,size(k,2));
for i = 1:length(k)
    [weight_matrix_null(:,i)] = FIB_null(noTap,micNum,k(i),fl,fh,f_step,fr,fs,x_array,d,c); % store the weights for null directions
end

%% Store steering vectors for CFIBN computation

steer =steer_vec(noTap,micNum,fl,fh,f_step,fs,x_array,d,c);

%% compute CFIBN weights 
disp('Processing: CFIBN design ...')
w_CFIBN = ((weight_matrix_null(:,1)'*steer).*(weight_matrix_null(:,2)'*steer).*(w_FIB'*steer))*steer'*(pinv(steer*steer')); % final CFIBN weight  eqn (21)
P       = plot_beampattern(w_CFIBN,noTap,micNum,200,3400,50,fs,x_array,d,c); % plot beampattern for CFIBN

%% Polar plots for selected frequencies
disp('Processing: Polar plot ...')
thetaFull = (0:pi/180:2*pi)';
Pfull = [P; flipud(P)]; 
Pfull=Pfull(1:end-1,:);
% figure;
% polar_dB(thetaFull, Pfull, 50, 10, 12)

limit_dB = -50;
index = Pfull <limit_dB;
Pfull(index)=limit_dB;
polarplot(thetaFull, Pfull(:,50),'LineWidth',1.2)
hold on
polarplot(thetaFull, Pfull(:,10),'LineWidth',1.2,'LineStyle','-.')
hold on
polarplot(thetaFull, Pfull(:,12),'LineWidth',1.2,'LineStyle','--')
rlim([-50 0])
thetalim([-180 180])
legend('800Hz','1600Hz','2700Hz','FontSize',15,'FontName','Times New Roman')
title('Polar beam pattern of FINB','FontSize',17,'FontName','Times New Roman')
%%
disp('Completed!')
 



















