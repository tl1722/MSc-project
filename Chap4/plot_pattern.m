%% plot beam pattern

N_window = 1024;                        
window = sqrt(hamming(N_window,'periodic'));
n_fft=N_window/2+1;
[s,f,t]=spectrogram(received(1,:),window,N_window/2,N_window,fs,'yaxis');
% plot the 3D beam pattern
theta =  180*[-1:0.01:1];
% B_mvdr=zeros(n_fft,length(theta));
B_lcmv=zeros(n_fft,length(theta));
for i=1:length(f)
  v = exp(1i*2*pi*[0:cfg.Nmics-1].'* d*cosd(theta)*f(i)/cfg.c);
  B_lcmv(i,:) = abs(total_w(:,i)' * v);

end

% 2D plot
figure;
P_lcmv=20*log10(B_lcmv./max(B_lcmv,[],2));
imagesc(theta',f',P_lcmv);
set(gca,'YDir','normal')
title('LCMV beam pattern','FontSize',17,'FontName','Times New Roman')
xlabel('DOA (Degree)','FontSize',15,'FontName','Times New Roman')
ylabel('Frequency (Hz)','FontSize',15,'FontName','Times New Roman')
colorbar
caxis([-50 0])
lim2=caxis;
% caxis(lim2);

% 3D plot
figure;
mesh(f,theta,P_lcmv.');
set(gca,YDir="reverse")
title('LCMV beam pattern')
xlabel('Frequency','FontSize',15,'FontName','Times New Roman')
ylabel('DOA','FontSize',15,'FontName','Times New Roman')
zlabel('Gain (dB)','FontSize',15,'FontName','Times New Roman')