function  Range = Spikelet_findNearestZeroCross(center_k,mag_k,D)
  from = center_k; to = center_k; 

  if D(center_k) > 0 && mag_k > 0
    for i = center_k:length(D)
      if D(i) <= 0
        to = i;  
        break;
      end
    end
    for i = center_k:-1:1
      if D(i) <= 0
        from = i;  
        break;
      end
    end
  end
  if D(center_k) < 0  && mag_k < 0
    for i = center_k:length(D)
      if D(i) >= 0
        to = i;  
        break;
      end
    end
    for i = center_k:-1:1
      if D(i) >= 0
        from = i;  
        break;
      end
    end
  end

  Range = [from,to];

end
