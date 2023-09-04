function [sir]=sir_test(desired_signal,interf_signal,fs)
  [lev,af,fso,vad]=v_activlev(desired_signal,fs,'d');
  P_d=0;
  for i=1:length(vad)
      if vad(i)==1
          P_d=P_d+desired_signal(i)^2;
      end
  end
  P_d=P_d/sum(vad(:)==1);

  [lev,af,fso,vad1]=v_activlev(interf_signal,fs,'d');
  P_i=0;
  for i=1:length(vad1)
      if vad(i)==1
          P_i=P_i+interf_signal(i)^2;
      end
  end
  P_i=P_i/sum(vad1(:)==1);
  
  sir=10*log10(P_d/P_i);


end
