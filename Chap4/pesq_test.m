function score=pesq_test(desired, received, fs)
  fs1=16000;
  [P,Q] = rat(fs1/fs);
  desired = resample(desired,P,Q);
  received = resample(received,P,Q);

writewav(received,16000,'received_signal.wav');
writewav(desired,16000,'desired_signal.wav');
score=pesq('desired_signal.wav','received_signal.wav');
end