function [Clen, Names] = Spikelet_Subseq_ConstantExtractionbyMagBand(D,MagThr)
  %{
  ----------------------
  D:   Data
  MagThr: Magnitude threshold
  ------------------

  %}
  %% 
  FuncName = 'Spikelet_Subseq_constantExtraction'; 
  Method = "leftRight_search"; % "rightOnly_search"
  
  %% output  
  Names = ["from","to","length"];
  Clen = NaN(length(D),length(Names));
  
  %% default value
  for i = 1:length(D)
    % if i == 2662; fprintf('[%s] debug\n',FuncName);  end
    if Method == "leftRight_search"
      FromTo = search_lr(i, D, MagThr);
      from = FromTo(1); to = FromTo(2);
    elseif Method == "rightOnly_search"
      count_right = search_right(i, D, MagThr);    
      from = i;
      to  = i+count_right-1;
    end
    Clen(i,:) = [from, to, to-from+1]; 
  end
  Clen = unique(Clen,'rows'); 
  [Length_sorted, Index_sorted] = sort(Clen(:,Names=="length"),"descend");
  Clen = Clen(Index_sorted,:);

end

function Range = search_lr(i, D, MagThr)
  count_right = search_right(i, D, MagThr);    
  count_left = search_left(i, D, MagThr);    
  from = i - (count_left -1);
  to = i + count_right -1;
  Range = [from,to];
end

function count = search_right(pos_i,D,MagThr)
  MagThr_half = MagThr/2;
  count = 0;
  for i = pos_i:length(D)
      if abs(D(i) - D(pos_i)) <= MagThr_half
        count = count + 1;
      else
        break;  
      end
  end
end

function count = search_left(pos_i,D,MagThr)
  MagThr_half = MagThr/2;
  count = 0;
  for i = pos_i:-1:1
      if abs(D(i) - D(pos_i)) <= MagThr_half
        count = count + 1;
      else
        break;  
      end
  end
end


