function MagInfo = Spikelet_Op_restrictSupportByLegSupportLength(MagInfo) 
  PLOT = true; % PLOT = false;
  
  FuncName = 'Spikelet_Op_restrictSupportByLegSupportLength';  % Magnitude window
  OpName = 'restrictSupportByLegSupportLength';
  %% parameters
  Param = MagInfo.param.operation.(OpName);
  SupportLength =  Param.length;
    
  %% Input and ouput
  % Data_org = MagInfo.data_org; 
  Data = MagInfo.data;
  Mag = MagInfo.magnitude;
  Left = MagInfo.left;
  Right = MagInfo.right;
  
  %% support restriction by Leg length
  Left_new = Left; 
  Right_new = Right;
  M_pos = find(Mag ~= 0);  
  
  for i = 1:length(M_pos)
    t = M_pos(i);
    mag_t = Mag(t); left_t = Left(t); right_t = Right(t); 
    leftSupp_t = t - left_t + 1; rightSupp_t = right_t - t + 1;
    
    Debug_On = false;
    if Debug_On && t == 274; 
      fprintf('[%s] debug t = %d\n',FuncName,t); 
      figure(); 
      plot((left_t:right_t),Data(left_t:right_t)); 
      title([num2str(left_t),' ',num2str(t),' ',num2str(right_t),' len:',num2str(right_t-left_t+1)]);
    end

    if leftSupp_t > SupportLength
      terminal_left = t-SupportLength+1;
      LeftRange = Data(terminal_left:t);  
      if mag_t > 0
        [val,pos_rel_1st] = min(LeftRange);
      elseif mag_t < 0
        [val,pos_rel_1st] = max(LeftRange);  
      end
      pos_rel = find(LeftRange==val,1,"first");  % As long as possible
      pos = terminal_left + pos_rel - 1;
      Left_new(t) = pos;
    end
    if rightSupp_t > SupportLength
      terminal_right = t+SupportLength-1;
      RightRange = Data(t:terminal_right);  
      if mag_t > 0
        [val,pos_rel_1st] = min(RightRange);
      elseif mag_t < 0
        [val,pos_rel_1st] = max(RightRange);  
      end
      pos_rel = find(RightRange==val,1,"last");  % As long as possible
      pos = t + pos_rel - 1;
      Right_new(t) = pos;
    end  
  end
  %% revise magnitude
  Mag_left_new = Data - Data(Left_new);
  Mag_right_new = Data - Data(Right_new);
  Mag_new = sign(Mag_left_new) .* min(abs(Mag_left_new), abs(Mag_right_new));

  %% set result
  MagInfo.left = Left_new;
  MagInfo.right = Right_new;
  MagInfo.magnitude = Mag_new;
  
  %% post processing
  MagInfo.param_str = ['support_length = ',num2str(Param.length)]; 
  MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param);
 
  %% PLOT
  if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
  end
end