function  [rate,num] = Spikelet_Subseq_doublePeak(subseq,SCstd)
  subseq_norm = normalize(subseq);
  std_1 = std(subseq_norm);
  Index_positive  = (subseq_norm > std_1*SCstd);
  Index_negative = (subseq_norm < -1*std_1*SCstd);
  Index_all = Index_positive | Index_negative;
  if sum(Index_all) == 0
    rate = 0;
    num = 0;
  else 
    rate = sum(Index_all)/length(subseq); 

    Diff_positive = diff(Index_positive);
    Diff_negative = diff(Index_negative);
    
    num_positive = sum(Diff_positive<0);
    num_negative = sum(Diff_negative<0);
    num = num_positive + num_negative;

    if Index_all(end) == true
      num = num + 1;
    end
  end
end
