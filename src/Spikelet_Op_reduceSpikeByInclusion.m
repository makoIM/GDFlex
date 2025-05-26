function MagInfo = Spikelet_Op_reduceSpikeByInclusion(MagInfo) 
  Plot_On = true; % Plot_On = false;
  
  FuncName = 'Spikelet_Op_reduceSpikeByInclusion';  % Magnitude window
  OpName = 'reduceSpikeByInclusion';
  %% parameters
  Param = MagInfo.param.operation.(OpName);
   
  %% Input and ouput
  % Data_org = MagInfo.data_org; 
  Data = MagInfo.data;
  Mag = MagInfo.magnitude;
  Left = MagInfo.left;
  Right = MagInfo.right;
   
  % restrict support
  Left_new = Left;
  Right_new = Right;
  M_pos = find(Mag ~= 0);
  for i = 1:length(M_pos)
    t = M_pos(i);
    left_t = Left(t); right_t = Right(t); 

    %fprintf('[%s] debug t = %d\n',FuncName,t); 
    Debug = false;
    if Debug && t == 2; 
      fprintf('[%s] debug t = %d\n',FuncName,t); 
      figure(); 
      plot((left_t:right_t),Data(left_t:right_t)); 
      title([num2str(left_t),' ',num2str(t),' ',num2str(right_t),' len:',num2str(right_t-left_t+1)]);
    end

    IncCand = find(Left <= left_t & right_t <= Right);
    if length(IncCand) >= 2
      Left_new(t) = left_t;
      Right_new(t) = right_t;
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
  MagInfo.param_str = [];
  MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param);

   %% PLOT
  if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
  end
end
