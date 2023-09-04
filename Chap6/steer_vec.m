function [s] =steer_vec(noTap,micNum,f_low,f_high,f_step,fs,x_array,d,c)
%%
% Allocating arrays

steer=[];

i1 = find(x_array == 1);  % positions of active sensors in thinned array
i1 = i1(:);
N = length(i1);           % number of active sensors
r = d*(i1(:)-i1(1));  % active sensor locations in m
r = r - r(N)/2;           % x = 0 in middle of sensor array
L = noTap;
  
 %% Calculating Steering vectors and store

s =[];
count = 0;
        
for freq = f_low:f_step:f_high;
    index = 0;
    count = count+1;
    for theta =0:180;
        index = index+1;
        x1 = exp(1i*(2*pi*freq)/c*cos(pi/180*theta)*r(:));
        x2 = [1 exp(1i*(2*pi*freq)/fs*(1:L-1))].';
        steer(:,index) = kron(eye(micNum),x2)*x1;
    end
        s = [s steer];
end

end
                
 
        
    
        
        

    
 

 