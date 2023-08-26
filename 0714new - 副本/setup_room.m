%% setup room parameters
cfg.room = [8 6 5];                % Room dimensions [x y z] (m)
cfg.n = 4096;                      % Number of samples
D=cfg.D;
theta_source=cfg.theta_source;
% source location 

cfg.position_source = D*[cos(theta_source/180*pi) sin(theta_source/180*pi) 0] + cfg.position_mics(1,:); 
% TDOA = sqrt(sum((bsxfun(@minus, position_source, position_mics)).^2, 2))/c; 

% generate room impulse response for the desired source using RIR generator
rir_desired = rir_generator(cfg.c, cfg.fs, cfg.position_mics, cfg.position_source, cfg.room, cfg.T60, cfg.n);
