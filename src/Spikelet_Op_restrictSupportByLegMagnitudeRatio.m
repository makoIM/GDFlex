function MagInfo = Spikelet_Op_restrictSupportByLegMagnitudeRatio(MagInfo) 
  Debug_On = true; Debug_On = false;
  FuncName = 'Spikelet_restrictSupportByLegMagnitudeRatio';  % Magnitude Time
  
  OpName = 'restrictSupportByLegMagnitudeRatio';
  %% parameters
  ParamSRM = MagInfo.param.operation.(OpName);
  MRatio = ParamSRM.magnitude_ratio; 
  
  %% Input and ouput
  % Data_org = MagInfo.data_org; 
  Data = MagInfo.data;
  Mag = MagInfo.magnitude;
  Center = find(MagInfo.magnitude ~= 0);
  Left = MagInfo.left;
  Right = MagInfo.right;
  MagLeft = zeros(length(Data),1);
  MagLeft(Center) = Data(Center)-Data(Left(Center));
  MagRight = zeros(length(Data),1);
  MagRight(Center) = Data(Right(Center))-Data(Center); 
  
  %% from here 
  %% support restriction by magnitude ( 8 fields)
  [Mag_new, Left_new, Right_new] = support_restriction_by_legMagnitude_ratio(Data,Mag,MagLeft,MagRight,Left,Right,ParamSRM);
  % revised feilds (magnitude, left, right,leg_magnitude )
  MagInfo.magnitude = Mag_new;
  MagInfo.left = Left_new;
  MagInfo.right = Right_new;
  % unvaryling feilds (center,value,type)
  % post-caluculated fields (data)
  
  %% post processing
  MagInfo.param_str = ['magnitude_ratio_threshold =', num2str(MRatio)];
  MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, ParamSRM);
 
 %% PLOT
 if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
 end
end

%% support_restriction_by_magnitude_ratio
function  [Mag, Left, Right] = support_restriction_by_legMagnitude_ratio(Data_org,Mag_org,MagLeft_org,MagRight_org,Left_org,Right_org,ParamSRM)
  FuncName = 'support_restriction_by_legMagnitude_ratio';
  
  % argument
  MRatio = ParamSRM.magnitude_ratio; 
  %MBand = ParamSRM.magnitude_band;    M_from = MBand(1); M_to = MBand(2); 
  
  Supp = [Left_org, Right_org];
  
  % output
  Left = Left_org; 
  Right = Right_org;
  
  %% main
  M_pos = find(Mag_org ~= 0);
 

  for i = 1:length(M_pos)
    t = M_pos(i); 
    %{
    if t == 140; 
      from = Left_org(t); to = Right_org(t); Len = to-from+1;  I = ceil(Len/4);
      figure(); plot(Data_org(from:to));  xlim([1 to-from+1]); xticks(1:I:Len); xticklabels(from:I:to);
      title(['[from,center,right]=[',num2str(from-from+1),'->',num2str(from)...
                                ,' ',num2str(t-from+1),'->',num2str(t)...
                                ,' ',num2str(to-from+1),'->',num2str(to),']']);
      fprintf('[%s] debug t = %d\n',FuncName,t);
    end
    %}
    % find left terminal
    Mag_org_rel = Mag_org(Supp(t,1):t-1);
    
    left_mag_rel = find(...
        ( Mag_org_rel*sign(MagLeft_org(t)) > 0 ) & sign(Mag_org(t)) == sign(Mag_org_rel) & ( abs(Mag_org_rel) > abs(MagLeft_org(t))*MRatio ) ...
         ,1,'last');
    if ~isempty(left_mag_rel) && ( Supp(t,1) + left_mag_rel - 1 ) ~= Left_org(t)
      left_mag = Supp(t,1) + left_mag_rel - 1;
   
      left_boundary = Right_org(left_mag);
      [left_val, left_rel] = min(Data_org(left_boundary:t)*sign(Mag_org(t)));
      Left(t) = left_boundary + left_rel - 1;
    end
    
    % find right terminal
    Mag_org_rel = Mag_org(t+1:Supp(t,2));
    right_boundary_rel = find(...
        ( Mag_org_rel*sign(MagRight_org(t)) > 0 ) & sign(Mag_org(t)) == sign(Mag_org_rel) & (abs(Mag_org_rel) > abs(MagRight_org(t))*MRatio)...
       ,1, 'first');
    if ~isempty(right_boundary_rel) && (t + 1 + right_boundary_rel - 1) ~= Right_org(t)
      right_mag = t + 1 + right_boundary_rel - 1;  
      right_boundary = Left_org(right_mag);
      [right_val, right_rel] = min(Data_org(t:right_boundary)*sign(Mag_org(t)));

      Right(t) = t + right_rel - 1;
    end
  end 
  %% revise magnitude
  Mag_left = Data_org - Data_org(Left);
  Mag_right = Data_org - Data_org(Right);
  Mag = sign(Mag_left) .* min(abs(Mag_left), abs(Mag_right));
end
