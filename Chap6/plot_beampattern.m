function [power_n] = plot_beampattern(w,noTap,micNum,f_low,f_high,f_step,fs,x_array,d,c)

% Plot CFIBN beampattern

i1 = find(x_array == 1);  
i1 = i1(:);
N = length(i1);           
r = d*(i1(:)-i1(1));  
r = r - r(N)/2;          
L = noTap;

count =0;

for freq = f_low:f_step:f_high;
    index = 0;
    count = count+1;
    for theta = 0:180
        index = index+1; 
        x1 = exp(1i*(2*pi*freq)/c*cos(pi/180*theta)*r(:));
        x2 = [1 exp(1i*(2*pi*freq)/fs*(1:L-1))].';
        steer = kron(eye(micNum),x2)*x1;          
        power(count,index) = ((w*(steer))^2);
    end
        
end
p = 10*log10((abs(power)));
theta = 0:180;
freq = (f_low:f_step:f_high);
figure;
mesh(theta,freq,p);
xlim([0 180]);
xlabel('DOA (degrees)');
ylabel('Frequency (Hz)');
zlabel('Gain (dB)');
 
%% Calculating response for few frequency - later used for polar plots
count = 0; 
freqzz = [800;1800;2700];
      
for i = 1:length(freqzz)
    index = 0;
    count = count+1;
    for theta = 0:180
        index = index+1;
        x1 = exp(1i*(2*pi*freqzz(i))/c*cos(pi/180*theta)*r(:));
        x2 = [1 exp(1i*(2*pi*freqzz(i))/fs*(1:L-1))].';
        steer = kron(eye(micNum),x2)*x1;
        power(count,index) = ((w*(steer))^2);
    end
end
power_n = power/max(max(abs(power)));
power_n = power_n';
power_n = (10*log10(abs(power_n)));
end
                
 
        
    
        
        

    
 

 