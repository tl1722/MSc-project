%% plot spectrograms, waveforms and beam pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_window = 1024;                        
window = sqrt(hamming(N_window,'periodic'));
inc=N_window/2;
% received signal in microphone array
t0=[1:1:length(received(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(received(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);
title('Received signal at reference Mic','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
caxis([-120 -40])
lim=caxis;

subplot(3,1,3)
plot(t0,received(1,:));
ylim([-0.8,0.8])
xlim([0,length(received(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%
% output signal
t0=[1:1:length(output_signal)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_signal,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);
title('output signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')

subplot(3,1,3)
plot(t0,output_signal);
ylim([-0.8,0.8])
xlim([0,length(output_signal)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

% desired component in received signal
t0=[1:1:length(desired(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(desired(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Desired component in received signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,desired(1,:));
ylim([-0.8,0.8])
xlim([0,length(desired(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% desired component in output signal
t0=[1:1:length(output_desired)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_desired,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Desired component in output','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,output_desired);
ylim([-0.8,0.8])
xlim([0,length(output_desired)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interf component in received signal
t0=[1:1:length(interf(1,:))]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(interf(1,:),window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Interference component in received signal','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,interf(1,:));
ylim([-0.8,0.8])
xlim([0,length(interf(1,:))/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interf component in output signal
t0=[1:1:length(output_interf)]/fs;
figure;
subplot(3,1,[1,2])
spectrogram(output_interf,window,inc,N_window,fs,'yaxis');
ylabel(colorbar,'dB/Hz')
xlabel([]);

caxis(lim);

title('Interference component in output','FontSize',17,'FontName','Times New Roman')
ylabel('Frequency (kHz)')
subplot(3,1,3)
plot(t0,output_interf);
ylim([-0.8,0.8])
xlim([0,length(output_interf)/fs])
xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % noise component in received signal
% t0=[1:1:length(noise(1,:))]/fs;
% figure;
% subplot(3,1,[1,2])
% spectrogram(noise(1,:),window,inc,N_window,fs,'yaxis');
% ylabel(colorbar,'dB/Hz')
% xlabel([]);
% 
% caxis(lim);
% 
% title('Noise component in received signal','FontSize',17,'FontName','Times New Roman')
% ylabel('Frequency (kHz)')
% subplot(3,1,3)
% plot(t0,noise(1,:));
% ylim([-0.8,0.8])
% xlim([0,length(noise(1,:))/fs])
% xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
% ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % noise component in output signal
% t0=[1:1:length(output_noise)]/fs;
% figure;
% subplot(3,1,[1,2])
% spectrogram(output_noise,window,inc,N_window,fs,'yaxis');
% ylabel(colorbar,'dB/Hz')
% xlabel([]);
% 
% caxis(lim);
% 
% title('Noise component in LCMV output','FontSize',17,'FontName','Times New Roman')
% ylabel('Frequency (kHz)')
% subplot(3,1,3)
% plot(t0,output_noise);
% ylim([-0.8,0.8])
% xlim([0,length(output_noise)/fs])
% xlabel('Time (s)','FontSize',15,'FontName','Times New Roman')
% ylabel('Amplitude','FontSize',15,'FontName','Times New Roman')
% 

