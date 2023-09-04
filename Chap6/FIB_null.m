function [weight_norm] = FIB_null(noTap,micNum,lkDir,f_low,f_high,f_step,ref_omega,fs,x_array,d,c)
%% Design of  FIB null

% Allocating arrays
s1=[];sx2=[];s_sum_total =zeros(micNum*noTap,micNum*noTap);

i1 = find(x_array == 1);  
i1 = i1(:);
N = length(i1);           
r = d*(i1(:)-i1(1)); 
r = r - r(N)/2;           

%% steering vector and big S for reference frequency for null direction

lkDir1 = [lkDir-5:lkDir+5];
L = noTap;
x1(:,1:length(lkDir1)) = exp(1i*(2*pi*ref_omega)/c*r(:)*cos(pi/180*lkDir1));     % vector of array input spectra
x2 = [1 exp(1i*(2*pi*ref_omega)/fs*(1:L-1))].'; 
s_sum_ref = kron(eye(micNum),x2)*x1;    % steering vector for look direction at reference freq

S_sum_ref = real((s_sum_ref)*(s_sum_ref')); % calculation of big S for reference freq 
                 
     
%% Calculation for undesired angle steering vectors (sidelobs) for reference freq
 
theta1 =  0:(lkDir-5);
theta2 = (lkDir+5):180;  
theta = [theta1,theta2];

x1(:,1:length(theta)) = exp(1i*(2*pi*ref_omega)/c*r(:)*cos(pi/180*theta));
s1(:,1:length(theta)) = kron(eye(micNum),x2)*x1(:,1:length(theta));                    
S1_omega_theta_ref = (sum(s1,2));   
clear s1;                 

%% steering vector calculation for all

theta = 0:180;
x3(:,1:length(theta)) = exp(1i*(2*pi*ref_omega)/c*r(:)*cos(pi/180*theta));
sx2(:,1:length(theta)) = kron(eye(micNum),x2)*x3(:,1:length(theta)); 
 
for omega = f_low:f_step:f_high; 
    x1 = exp(1i*(2*pi*omega)/c*r(:)*cos(pi/180*theta));
    x_new = [1 exp(1i*(2*pi*omega)/fs*(1:L-1))].'; 
    sx1 = kron(eye(micNum),x_new)*x1;                                             
    s_sum_total =  s_sum_total+ real((sx1-sx2)*(sx1-sx2)');
    clear sx1;
end           

%% Calculation of 'R' matrix and 'p'  and weights

beta = 0.95;  
R    = S_sum_ref + beta.*s_sum_total; % eqn (18)
p    = S1_omega_theta_ref;            % eqn (19)
w    = (pinv(R)*p)/(p'*(pinv(R)*p));  % eqn (17)
  
 %% Calculating response
 
count = 0;
      
for freq = f_low:f_step:f_high;
    index = 0;
    count = count+1;
    for theta = 0:180
        index = index+1;
        x1 = exp(1i*(2*pi*freq)/c*cos(pi/180*theta)*r(:));
        x2 = [1 exp(1i*(2*pi*freq)/fs*(1:L-1))].';
        steer = kron(eye(micNum),x2)*x1;          
        resp(count,index) =  ((w'*(steer)));        
    end
end

value = max(max((resp)));
weight_norm = w./value;

count = 0;    
for freq = f_low:f_step:f_high;
    index = 0;
    count = count+1;
    for theta = 0:180
        index = index+1;
        x1 = exp(1i*(2*pi*freq)/c*cos(pi/180*theta)*r(:));
        x2 = [1 exp(1i*(2*pi*freq)/fs*(1:L-1))].';
        steer = kron(eye(micNum),x2)*x1;          
        power(count,index) = ((weight_norm'*(steer))^2);
    end
end

p = (10*log10(abs(power)));
theta = 0:180;
freq = (f_low:f_step:f_high);
figure;
mesh(theta,freq,p);
xlim([0 180]);
xlabel('DOA (degrees)');
ylabel('Frequency (Hz)');
zlabel('Gain (dB)');

end
                
 
        
    
        
        

    
 

 