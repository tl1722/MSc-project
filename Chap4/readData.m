% read speech signals from ieee dataset
% input: 
%   folder: folder name of dataset
% output:
%   audioData: a cell which includes 100 speech signals
%   fs: sampling frequency
% TONGKAI LI, 2023.08

function [audioData,fs]=readData(folder)
audioData = cell(100, 1);
k=1;
for i = 1:10
    for j = 1:10
        filename = sprintf('ieee%02dm%02d.wav', i, j);
        fullFilePath = fullfile(folder, filename);
       
        [audio, fs] = audioread(fullFilePath);
       
        audioData{k,1} = audio;
        
%         fprintf('Reading %s\n', fullFilePath);
        k=k+1;
    end
end
end