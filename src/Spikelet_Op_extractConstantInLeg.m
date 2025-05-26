function MagInfo = Spikelet_Op_extractConstantInLeg(MagInfo)
  PLOT = true; % PLOT = false;

  %%
  FuncName = 'Spikelet_Op_extractConstant_inLeg';
  OpName = 'extractConstantInLeg';

  %% parameters
  Param = MagInfo.param.operation.(OpName);
  CLT = Param.constant_length_threshold; 
  CMT = Param.constant_magnitude_threshold;
  LimitLen = inf;

  %% argument
  D = MagInfo.data;
  D_org = MagInfo.data_org;
  Mag = MagInfo.magnitude;
  Center = find(MagInfo.magnitude ~= 0);
  Left = MagInfo.left;
  Right = MagInfo.right;

  Mag_new = Mag;
  Left_new = Left;
  Right_new = Right;

  %% maimal constant Region
  [MCreg, Names] = Spikelet_extractConsantSegment(D_org,CLT,CMT,LimitLen);

  %% revise terminal of spikes
  % for each constant segment
  for i = 1:size(MCreg,1)
    left_i = MCreg(i,Names=="left"); right_i = MCreg(i,Names=="right");  
    mean_i = MCreg(i,Names=="mean");

    % restrict the support of a left leg 
    % for each left leg intersecting with constant segment_i
    for j = 1:length(Center)
      center_j = Center(j);
      left_j = Left(center_j);
      sign_j = Mag(center_j);

      %{
      DEBUG = true;
      fprintf('[debug] center=%d (i,j)=(%d,%d)\n',center_j ,i,j);
      % if DEBUG && i == 13 && j == 76
      if DEBUG && ismember(center_j, [4152,4155,4200,4202,4248]) && 
        fprintf('[debug] center=%d (i,j)=(%d,%d)\n',center_j ,i,j);
      end
      %}
      if left_j <= right_i && right_i <= center_j
        LeftLeg_j = (left_j:center_j);  
        if sign_j > 0
          left_rel = find(D_org(LeftLeg_j) <= mean_i,1,'last');
        else
          left_rel = find(D_org(LeftLeg_j) >= mean_i,1,'last');
        end
        
        if ~isempty(left_rel)
          %{  
          if center_j == 10334
            fprintf('[%s] debug\n',FuncName);
          end
          %}
          Left_new(center_j) = left_j + left_rel - 1;   
        end
      end
    end
    % restrict the support of a right leg 
    % for each right leg intersecting with constant segment_i
    for j = 1:length(Center)
      
      center_j = Center(j);
      sign_j = Mag(center_j);
      right_j = Right(center_j);
      %{
      DEBUG = true;
      fprintf('[debug] (i,j)=(%d,%d)\n',i,j);
      if DEBUG && center_j == 10018 && left_i == 10100 % i == 13 && j == 73
        figure();   
        fprintf('[debug] (i,j)=(%d,%d)\n',i,j);
        Title = ['j:',num2str(center_j),' ',num2str(right_j)]; 
        PlotRange = (center_j:max(right_j,right_i)); plot(PlotRange,D_org(PlotRange)); 
        ConstRange = (left_i:right_i);  hold on; plot(ConstRange,mean_i*ones(length(ConstRange),1));
        xline(right_j);
        title(Title);    
      end
      %}
      if center_j <= left_i && left_i <= right_j
        RightLeg_j = (center_j:right_j);     
        if sign_j > 0
          right_rel = find(D_org(RightLeg_j) <= mean_i,1,'first');
        else
          right_rel = find(D_org(RightLeg_j) >= mean_i,1,'first');
        end
        if ~isempty(right_rel)
          Right_new(center_j) = center_j + right_rel - 1;   
        end
      end       
    end
  end

  %% set result
  MagInfo.left = Left_new;
  MagInfo.right = Right_new;
  MagInfo.magnitude = Mag_new;
  
  %% post processing
  MagInfo.param_str = ['CLT:',num2str(CLT), ' CMT:', num2str(CMT)];
  MagInfo.constant = MCreg;      % Operation depended field
  MagInfo.constant_names = Names;

  MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param);

  MagInfo.output.(OpName).constant = MCreg;  % Operation depended field
  MagInfo.output.(OpName).constant_names = Names;

   %% PLOT
  if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
  end
  
end


%{
function [Cmax, Index_AC] = maximal_disjoint(len, Clen,D,Names,Index_AC)
   Index_len = Clen(:,Names=="length") == len;
   Clen_len = Clen(Index_len,:);
   % Clen_l_cand = disjoint_list(Clen_len,Names,Index_AC);
   [Cmax, Index_AC] = disjoint_and_maximal(Clen_len,D,Names,Index_AC);
end
%}

