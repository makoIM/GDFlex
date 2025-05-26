function [OVERLAP, overlap_length,ov_type] = Spikelet_overlapSegment(range_1,range_2)
  from_1 = range_1(1);
  to_1 = range_1(2);
  from_2 = range_2(1);
  to_2 = range_2(2);

  if (from_2 <= from_1 && to_1 <= to_2)      % range_2 includes range_1
    OVERLAP = true;  
    overlap_length = to_1 - from_1 + 1; 
    ov_type = "including";
  elseif (from_1 <= from_2 && to_2 <= to_1)      % range_1 includes range_2
    OVERLAP = true;  
    overlap_length = to_2 - from_2 + 1; 
    ov_type = "included";
  elseif (from_2 <= from_1 && from_1 <= to_2)     % range_2 proceeds
    OVERLAP = true;  
    overlap_length = to_2 - from_1 + 1;
    ov_type = "front";
  elseif (from_1 <= from_2 && from_2 <= to_1) % range_1 proceeds
    OVERLAP = true;  
    overlap_length = to_1 - from_2 + 1;
    ov_type = "back";
  else
    OVERLAP = false;
    overlap_length = 0;
    ov_type = "";
  end
end