function [rateMax,numLast] = Spikelet_Subseq_peakDiffOpposite(subseq,SCstd)
  subseq_norm = normalize(subseq);
  std_1 = std(subseq_norm); 

  MagInfo = Spikelet_generateInitialSpikelet_ver_02(subseq_norm,Inf);
  Mag = MagInfo.magnitude;
  Index_peak = find(abs(Mag) > std_1*(SCstd*2));
  if isempty(Index_peak)
    rateMax = 0;
    numLast = 0;
    return;
  else
    numMax = 0; rateMax = 0;
    for i = 1:length(Index_peak)
      peak_i = Index_peak(i); 
      sign_i = sign(Mag(peak_i));
      if sign_i == 1
        Index = subseq_norm < -1*std_1*(SCstd-1); % opposite
      elseif sign_i == -1
        Index = subseq_norm > 1*std_1*(SCstd-1);  
      end
      rate = sum(Index)/length(subseq);
      Diff = diff(Index);
      num = sum(Diff<0);
      if Index(end) == true
        num = num + 1;  
      end
      if num > numMax
        numMax = num; 
        rateMax = rate;
      end
    end
    numLast = numMax + 1;
  end
  
end