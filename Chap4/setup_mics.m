%% setup microphone parameters
cfg.c = 340;                  % Sound velocity (m/s) 
cfg.fs=fs;
d=cfg.d;
% Microphone array position [x y z] (m) 
cfg.position_mics = [4 3 2;  
                 4+d 3 2; 
                 4+2*d 3 2; 
                 4+3*d 3 2; 
                 4+4*d 3 2; 
                 4+5*d 3 2; 
                 4+6*d 3 2; 
                 4+7*d 3 2];  

% distance matrix of mics
cfg.l_mn = zeros(cfg.Nmics, cfg.Nmics);
for i = 1:cfg.Nmics
    for j = 1:cfg.Nmics
        cfg.l_mn(i, j) = abs(i - j) * d;
    end
end

