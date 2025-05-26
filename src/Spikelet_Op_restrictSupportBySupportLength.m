function MagInfo = Spikelet_Op_restrictSupportBySupportLength(MagInfo) 
  Debug_On = true; Debug_On = false;
  
  FuncName = 'Spikelet_restrictSupportBySupportLength';  % Magnitude window
  OpName = 'restrictSupportBySupportLength';
  %% parameters
  Param = MagInfo.param.operation.(OpName);
  if isfield(Param,'query')
    M_formula = ParamSRB.query;
  else
    M_formula = 'Mag ~= 0';  % Unconditionally true
  end
  SupportLength =  Param.length;
    
  %% Input and ouput
  Data = MagInfo.data; 
  Mag = MagInfo.magnitude;
  Left = MagInfo.left;
  Right = MagInfo.right;
  
  %% support restriction by window
  % output
  Left_new = Left; 
  Right_new = Right;
  
  % find target times 
  Time = (1:length(Mag))';
  Supp = Right - Left + 1;
  
  Index_M = eval(M_formula);
  Index_MW = Index_M & (Supp > SupportLength);
  
  M_pos = find(Index_MW == true);
  
  % restrict support
  for i = 1:length(M_pos)
    t = M_pos(i);
    left_t = Left(t); right_t = Right(t); 
    if Debug_On && t == 274; 
      fprintf('[%s] debug t = %d\n',FuncName,t); 
      figure(); 
      plot((left_t:right_t),Data(left_t:right_t)); 
      title([num2str(left_t),' ',num2str(t),' ',num2str(right_t),' len:',num2str(right_t-left_t+1)]);
    end

    MaxMag_i = 0; 
    for from_i = max(t - SupportLength + 2,left_t): min(t-1,right_t)
      [Left_ij, Left_ij_rel] = max(abs(Data(t) - Data(from_i:t-1)));
      [Right_ij,Right_ij_rel] = max(abs(Data(t+1:min(from_i+SupportLength-1,right_t)) - Data(t)));
      Mag_ij = min(Left_ij,Right_ij); 
      if Mag_ij > MaxMag_i
        MaxMag_i = Mag_ij;  
        MaxLeft_ij = from_i + Left_ij_rel - 1;
        MaxRight_ij = t+1 + Right_ij_rel - 1;
      end
    end
    Left_new(t) = MaxLeft_ij; 
    Right_new(t) = MaxRight_ij; 
   
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
  MagInfo.param_str = [OpName, ' ',num2str(Param.length)];
  MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param);
  
 %% debug
 if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
 end
end