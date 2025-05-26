function MagInfo = Spikelet_Op_extendLegSupportMaximal(MagInfo) 
  PLOT = true; PLOT = false;
  DEBUG = true; DEBUG = false;
  FuncName = 'Spikelet_Op_extendLegSupportMaximal';  % Magnitude window
  OpName = 'extendLegSupportMaximal';
  %% parameters
  Param = MagInfo.param.operation.(OpName);
  if isfield(Param,'length')
    SuppLen = Param.length;
  else
    SuppLen = inf;  
  end

  %% Input and ouput
  D = MagInfo.data; 
  Mag = MagInfo.magnitude;
  Center = find(MagInfo.magnitude ~= 0);
  Left = MagInfo.left;
  Right = MagInfo.right;
  LeftPoint = Left(Center);
  RightPoint = Right(Center);
  ImpPoint = [LeftPoint; Center; RightPoint];
  ImpPoint = sort(unique(ImpPoint),'ascend');

  %% support restriction by window
  % output
  Left_new = Left; 
  Right_new = Right;
  
  %% kokokara
  % find target times 
  for i = 1:length(Center)
    center_i = Center(i);
    CenterImpId = find(ImpPoint == center_i);
    centerValue_i = D(center_i);
    centerSign_i = sign(Mag(center_i));
    %{
    if center_i == 4152
      fprintf('debug\n');
    end
    %}
    if i ~= 1
      ImpId_rightTerminal_inLeftLeg_i = find(ImpPoint == center_i) - 1;
      % LeftTime_i
      if center_i - 1 + 1 > SuppLen
        BoundaryTime_i = max(center_i - SuppLen + 1, 1); 
        BoundaryImpId_i = find(ImpPoint <= BoundaryTime_i,1,'last');
        if BoundaryTime_i == ImpPoint(BoundaryImpId_i)
          LeftTime_i = (ImpPoint(BoundaryImpId_i:ImpId_rightTerminal_inLeftLeg_i));  
        else
          LeftTime_i = [(BoundaryTime_i:ImpPoint(BoundaryImpId_i+1)-1)'; (ImpPoint(BoundaryImpId_i+1:ImpId_rightTerminal_inLeftLeg_i))]; 
        end
      else
        LeftTime_i = ImpPoint(1:ImpId_rightTerminal_inLeftLeg_i);
      end
      % LeftValue_i
      LeftValue_i = D(LeftTime_i);     
    end
    % RightTime_i
    if i ~= length(Center)
      impId_leftTerminal_inRightLeg_i = find(ImpPoint == center_i) + 1;
      % RightTime_i
      if ImpPoint(end) - center_i + 1 > SuppLen
        BoundaryTime_i = min(center_i + SuppLen -1, length(D));
        BoundaryImpId_i =  find(ImpPoint >= BoundaryTime_i,1,'first');
        if BoundaryTime_i == ImpPoint(BoundaryImpId_i)
          rightTime_i = ImpPoint(impId_leftTerminal_inRightLeg_i:BoundaryImpId_i);  
        else
          rightTime_i = [ImpPoint(impId_leftTerminal_inRightLeg_i:BoundaryImpId_i-1);(ImpPoint(BoundaryImpId_i-1)+1:BoundaryTime_i)']; 
        end
      else
         rightTime_i = ImpPoint(impId_leftTerminal_inRightLeg_i:end);
      end
      % RightValue_i
      RightValue_i = D(rightTime_i);
    end

    if centerSign_i == 1
      PN = "positive";
      if i ~= 1  
        [leftTime_i, leftOverTime_i] = find_left_maximal(LeftValue_i,LeftTime_i,centerValue_i,PN);
      end
      if i ~= length(Center)
        [rightTime_i, rightOverTime_i] = find_right_maximal(RightValue_i,rightTime_i,centerValue_i,PN);
      end
    elseif centerSign_i == -1
      PN = "negative";  
      if i ~= 1
        [leftTime_i, leftOverTime_i] = find_left_maximal(LeftValue_i,LeftTime_i,centerValue_i,PN);
      end
      if i ~= length(Center)
        [rightTime_i, rightOverTime_i] = find_right_maximal(RightValue_i,rightTime_i,centerValue_i,PN);
      end
    end

    if i ~= 1 
      Left_new(center_i) = leftTime_i;
      
      if DEBUG && leftTime_i < Left(center_i)
        LocalPlot = true;  
        fprintf('[%s] extend left leg %d->%d at center(%d)\n',FuncName,Left(center_i),leftTime_i,center_i);  
      elseif  DEBUG && leftTime_i > Left(center_i)
        LocalPlot = true;   
        fprintf('[%s] restrict left %d->%d at center(%d)\n',FuncName,Left(center_i),leftTime_i,center_i);
      end
    end
    if i ~= length(Center)
      Right_new(center_i) = rightTime_i;

      if DEBUG && Right(center_i) < rightTime_i
        LocalPlot = true;   
        fprintf('[%s] extend right leg %d->%d at center(%d)\n',FuncName,Right(center_i),rightTime_i,center_i);  
      elseif  DEBUG && rightTime_i > Right(center_i)
        LocalPlot = true;   
        fprintf('[%s] restrict right %d->%d at center(%d)\n',FuncName,Right(center_i),rightTime_i,center_i);
      end    
    end

    if exist('LocalPlot','var') && LocalPlot
      if exist('leftOverTime_i','var')
          from = leftOverTime_i;  
      else
          from = Left_new(center_i);  
      end
      if exist('rightOverTime_i','var')
          to = rightOverTime_i;
      else
          to = Right_new(center_i);
      end  
      
      figure(); 
      fromMin = min([from,Left(center_i),Left_new(center_i)]);
      toMax = max([to,Right(center_i),Right_new(center_i)]);
      XLIM = lim_margin([fromMin,toMax],0.05);
      YLIM = lim_margin(D(fromMin:toMax),0.05);

      subplot(2,1,1)
        plot((Left(center_i):Right(center_i)),D(Left(center_i):Right(center_i))); xlim(XLIM); ylim(YLIM);
        hold on;
        xline(Left(center_i),'--'); xline(Right(center_i),'--');
        xline(center_i,':');
        title([num2str(Left(center_i)),'-',num2str(Right(center_i))]);
      subplot(2,1,2)
        plot((from:to),D(from:to)); xlim(XLIM); title('first over time'); xlim(XLIM); ylim(YLIM);
        hold on;
        xline(Left_new(center_i),'--'); xline(Right_new(center_i),'--');
        xline(center_i,':');
        title([num2str(Left_new(from)),'-',num2str(Right_new(to))]);

      Title = [char(PN),'(',num2str(center_i),')',' SuppLen=',num2str(SuppLen)];  
      Fig_figtitle(Title);  

      LocalPlot = false;
    end
    
  end
  %% revise magnitude
  Mag_left_new = D - D(Left_new);
  Mag_right_new = D - D(Right_new);
  Mag_new = sign(Mag_left_new) .* min(abs(Mag_left_new), abs(Mag_right_new));
  
  %% set result
  MagInfo.left = Left_new;
  MagInfo.right = Right_new;
  MagInfo.magnitude = Mag_new;
  
  %% post processing
  MagInfo.param_str = OpName;
  MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param);
  
 %% debug
 if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
 end
end

function [leftTime, leftOverTime] = find_left_maximal(LeftValue, LeftTime,centerValue,PN)
  for i = length(LeftValue):-1:1
    pos_i = i;   
    if PN == "positive" &&LeftValue(i) > centerValue
      break;
    end
    if PN == "negative" &&LeftValue(i) < centerValue
      break;
    end
  end

  RangeVal = LeftValue((pos_i:end));
  if PN == "positive" 
    [MinVal, Pos_rel] = min(RangeVal);
  elseif PN == "negative"
    [MaxVal, Pos_rel] = max(RangeVal);
  end
  leftCenterId = pos_i + Pos_rel -1; 
  
  % leftTime
  leftTime = LeftTime(leftCenterId);
  leftOverTime =  LeftTime(pos_i);
end

function [rightTime, rightOverTime] = find_right_maximal(RightValue, RightTime,centerValue, PN)
  for i = 1:length(RightValue)
    pos_i = i;   
    if PN == "positive" && RightValue(i) > centerValue 
      break;
    end
    if PN == "negative" && RightValue(i) < centerValue
      break;
    end
  end

  RangeVal = RightValue((1:pos_i));
  if PN == "positive" 
    [MinVal, Pos] = min(RangeVal);
  elseif PN == "negative"
    [MinVal, Pos] = max(RangeVal);
  end
  rightCenterId = Pos;

  % leftTime
  rightTime = RightTime(rightCenterId);
  rightOverTime = RightTime(pos_i);
end

function LIM = lim_margin(D,Rate)
  Margin = Rate*(max(D) - min(D));
  if Margin ~= 0
    LIM = [min(D)-Margin, max(D)+Margin];
  else
    LIM = [];
  end
end
