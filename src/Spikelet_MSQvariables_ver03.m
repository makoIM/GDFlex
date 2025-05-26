function MsqV = Spikelet_MSQvariables_ver03(MagInfo, PsrQuery)
  FuncName = 'Spikelet_MSQvariables_ver03';
  % PLOT = true;
 %% data
  Data_org = MagInfo.data_org; 
  Mag = MagInfo.magnitude;
  Peak_time = find(Mag ~= 0);
  Left = MagInfo.left;
  Right = MagInfo.right;
  Mag_left = Data_org - Data_org(Left);
  Mag_right = Data_org - Data_org(Right);

  %% dimension reduction by spikelet
  MagS   =  Mag(Peak_time);
  LeftS = Left(Peak_time);
  RightS = Right(Peak_time);
  SuppS  = RightS - LeftS + 1;
  Supp_leftS = Peak_time - LeftS + 1;
  Supp_rightS = RightS - Peak_time + 1;
  Supp_LRratio = Supp_rightS./Supp_leftS;
  Mag_leftS = Mag_left(Peak_time);
  Mag_rightS = Mag_right(Peak_time);
  Mag_lrRatio =  sign(Mag_rightS).*abs(Mag_rightS./Mag_leftS);

  %% Quartile support length
  qpercent = [0.25,0.5,0.75,1.00];   
  
  %% both support
  Quartile_length_both = zeros(length(MagS),length(qpercent)); 
  Quartile_length_both_total = zeros(length(MagS),length(qpercent)); 
  Quartile_from_both =  repmat((1:length(MagS))', 1, length(qpercent)); 
  Quartile_to_both = repmat((1:length(MagS))', 1, length(qpercent));    
  for i = 1:length(Peak_time)
    t = Peak_time(i);
    % if i == 8; fprintf('[%s] debug',FuncName); end
    mag = Mag(t);
    left = Left(t); right = Right(t);
    support_backward = (left:t);
    support_forward = (t:right);
    support_both = (left:right);
    for j = 1:length(qpercent) 
      DEBUG = false;
      if DEBUG
        fprintf('[%s] debug',FuncName); 
        if j == 1;  
          figure(); plot(Data_org); xlim([left right]); title(num2str(t));  
          hold on; 
          if sign(mag) > 0
            TH = Data_org(t) - qpercent(j)*abs(mag);
          else
            TH = Data_org(t) + qpercent(j)*abs(mag);  
          end
          yline(TH);
         end
      end

      back_rel    = find( abs(Data_org(t) - Data_org(support_backward)) <= qpercent(j)*abs(mag), 1, 'first');
      forward_rel = find( abs(Data_org(t) - Data_org(support_forward))  <= qpercent(j)*abs(mag), 1, 'last');
      Index =           (abs(Data_org(t) - Data_org(support_both))      <= qpercent(j)*abs(mag));
      from = left + (back_rel+1) -1;
      to =    t + (forward_rel - 1) - 1;
      Quartile_length_both(i,j) =  to - from + 1;
      Quartile_length_both_total(i,j) =  sum(Index);
      Quartile_from_both(i,j)   =  from;
      Quartile_to_both(i,j)  =  to;
    end
  end
  %{
  %% left support
  Quartile_length_left = zeros(length(MagS),length(qpercent)+1); 
  Quartile_length_left_total = zeros(length(MagS),length(qpercent)+1); 
  Quartile_from_left =  repmat((1:length(MagS))', 1, length(qpercent)+1); 
  Quartile_to_left = repmat((1:length(MagS))', 1, length(qpercent)+1);   
  
   for i = 1:length(Peak_time)
    t = Peak_time(i);
    % if (t == 2068 || t == 2067);fprintf('[%s] debug',FuncName);  end
    mag = Mag_left(t);
    left = Left(t);
    support_backward = (left:t);
    for j = 1:length(qpercent) 
      back_rel    = find( abs(Data_org(t) - Data_org(support_backward)) >= qpercent(j)*abs(mag), 1, 'last');
      Index       =     ( abs(Data_org(t) - Data_org(support_backward)) >= qpercent(j)*abs(mag));
      from = left + (back_rel+1) -1;
      Quartile_length_left(i,j) =  t - from + 1;
      Quartile_length_left_total(i,j) = sum(Index);
      Quartile_from_left(i,j)   =  from;
      Quartile_to_left(i,j)   =  t; 
    end
   end
  %% right support
  Quartile_length_right = zeros(length(MagS),length(qpercent)+1); 
  Quartile_length_right_total = zeros(length(MagS),length(qpercent)+1); 
  Quartile_from_right =  repmat((1:length(MagS))', 1, length(qpercent)+1); 
  Quartile_to_right = repmat((1:length(MagS))', 1, length(qpercent)+1);   
  
   for i = 1:length(Peak_time)
    t = Peak_time(i);
    % if (t == 2068 || t == 2067);fprintf('[%s] debug',FuncName);  end
    mag = Mag_right(t);
    right = Right(t);
    support_forward = (t:right);
    for j = 1:length(qpercent) 
      forward_rel = find( abs(Data_org(t) - Data_org(support_forward)) >= qpercent(j)*abs(mag), 1, 'first');
      Index =           ( abs(Data_org(t) - Data_org(support_forward)) >= qpercent(j)*abs(mag));
      to =    t + (forward_rel - 1) - 1;
      Quartile_length_right(i,j) =  to - t + 1;
      Quartile_length_right_total(i,j) =  sum(Index);
      Quartile_from_right(i,j)   =  t;
      Quartile_to_right(i,j)   = to; 
    end  
   end
  
  %% 100 qpercent
  qpercent_full = [qpercent,100];

  %% both
  Quartile_length_both(:,length(qpercent_full))    = RightS - LeftS + 1;
  Quartile_from_both(:,length(qpercent_full))      = LeftS;
  Quartile_to_both(:,length(qpercent_full))        = RightS;
  
  Quartile_length_left(:,length(qpercent_full))    = Peak_time - LeftS+1;
  Quartile_from_left(:,length(qpercent_full))      = LeftS;
  Quartile_to_left(:,length(qpercent_full))        = Peak_time;
  
  Quartile_length_right(:,length(qpercent_full))   = RightS - Peak_time+1;
  Quartile_from_right(:,length(qpercent_full))     = Peak_time;
  Quartile_to_right(:,length(qpercent_full))       = RightS;
  
  %% quantile suport
  Quartile_support_def = cell(1,           length(qpercent_full)*2);
  Quartile_support_both = NaN(length(MagS), length(qpercent_full)*2);
  Quartile_support_left = NaN(length(MagS), length(qpercent_full)*2);
  Quartile_support_right= NaN(length(MagS), length(qpercent_full)*2);
  
  % both
  for i = 1:length(qpercent_full)
    Quartile_support_def{2*(i-1)+1} = [num2str(qpercent_full(i)),'_left']; 
    Quartile_support_def{2*(i-1)+2}   = [num2str(qpercent_full(i)),'_right']; 
    
    Quartile_support_both(:, 2*(i-1)+1)  = Quartile_from_both(:,i);
    Quartile_support_both(:, 2*(i-1)+2)  = Quartile_to_both(:,i);
    
    Quartile_support_left(:, 2*(i-1)+1)  = Quartile_from_left(:,i);
    Quartile_support_left(:, 2*(i-1)+2)  = Quartile_to_left(:,i);
    
    Quartile_support_right(:, 2*(i-1)+1)  = Quartile_from_right(:,i);
    Quartile_support_right(:, 2*(i-1)+2)  = Quartile_to_right(:,i);
  end
  %}
  %% output
  MsqV = struct;
  MsqV.meta.data_org = Data_org;
  MsqV.meta.qpercent = qpercent;
  MsqV.Matrix = [];
  MsqV.Names = [];
  %  MsqV.magnitude = Mag;

  MsqV.Matrix = [ Peak_time...
                , MagS, Mag_rightS, Mag_leftS, Mag_lrRatio...
                , Supp_leftS, Supp_rightS,Supp_LRratio,SuppS...
                , Data_org(Peak_time), Data_org(LeftS),Data_org(RightS),Left(Peak_time),Right(Peak_time)...
                , Quartile_length_both_total(:,1)./Quartile_length_both_total(:,2)...
                , Quartile_length_both_total(:,1)./Quartile_length_both_total(:,3)...
                , Quartile_length_both_total(:,1)./Quartile_length_both_total(:,4)...
                , Quartile_length_both_total(:,2)./Quartile_length_both_total(:,3)...
                , Quartile_length_both_total(:,2)./Quartile_length_both_total(:,4)...
                , Quartile_length_both_total(:,3)./Quartile_length_both_total(:,4)...
                %{
                , Quartile_length_both_total, Quartile_length_both...
                , Quartile_length_left_total...
                , Quartile_length_left...
                , Quartile_length_right...
                , Quartile_length_right_total...  
                , Quartile_support_both...
                , Quartile_support_left...
                , Quartile_support_right
                %}
                ];

  MsqV.Names =  ["peak_time"...
                ,"spike_magnitude","spike_magnitude_right","spike_magnitude_left","spike_magnitude_lrRatio"...
                ,"spike_support_left","spike_support_right","spike_support_lrRatio","spike_support"...
                ,"peak_val","left_val","right_val","left","right"...
                ,"quantile_length_both_ratio_12"...
                ,"quantile_length_both_ratio_13"... 
                ,"quantile_length_both_ratio_14"...
                ,"quantile_length_both_ratio_23"...
                ,"quantile_length_both_ratio_24"...
                ,"quantile_length_both_ratio_34"...
                %{
                ,"quantile_length_both_total","quantile_length_both"...  
                ,"quantile_length_left_total"...
                ,"quantile_length_left"...
                ,"quantile_length_right"...
                ,"quantile_length_right_total"...  
                ,"quantile_support_both"... 
                ,"quantile_support_left"...
                ,"quantile_support_right"
                %}
                ];

   if ~isempty(PsrQuery) && contains(PsrQuery,'Mag_OpSpike')
    Mag_OpSpike = mag_OpSpikeSqRatio(MagInfo,Peak_time,MagS,LeftS,RightS);   
    MsqV.Matrix = [MsqV.Matrix, Mag_OpSpike];    
    MsqV.Names = [MsqV.Names,"spike_magnitude_OpSpikeSqRatio"];
   end

   %% plot
   if exist('PLOT','var') && PLOT
     D  = MagInfo.data_org; 
     M = MagInfo.magnitude;
     L  = MagInfo.left;
     R = MagInfo.right;
     X_Range =  [1 length(D)];  X = (1:length(D)); LineWidth = 0.75; Sign = 'positive';
     % Rate = 1.05;  Shift = 1.0;
     % Y_LIM = [min(D)*Rate,max(D)*Rate];  Y_Width = Y_LIM(2) - Y_LIM(1); 
     % Y_LIM_P =  [0-Shift,Y_Width-Shift]; 

     figure(); Row = 2;
     subplot(Row,1,1)
       plot(X,D); xlim(X_Range); % ylim(Y_LIM_P);
     subplot(Row,1,2)
       Op_PLOT =MagInfo.param.operation.operation_sequence{end};
       Spikelet_plot_decomposition_ver_04(M,L,R);
       title(Op_PLOT);
   end
end

function Mag_opsq = mag_OpSpikeSqRatio(MagInfo,Peak_time,MagS,LeftS,RightS)
    InitialSpike = MagInfo.output.generateInitialSpikelet;
    PI = find(InitialSpike.magnitude ~= 0);
    MI = InitialSpike.magnitude(PI);
    TI = sign(MI);
    LI = InitialSpike.left(PI);
    RI = InitialSpike.right(PI);
    Mag_opsq = zeros(length(MagS),1);
    for i = 1:length(Peak_time)
      t_i = sign(MagS(i));
      l_i = LeftS(i);
      r_i = RightS(i);
      Index = l_i <= LI &  RI <= r_i & TI == t_i*(-1);
      mlist_i = MI(Index);
      if isempty(mlist_i)
        Mag_opsq(i) = 0;   
      else
        Mag_opsq(i) = sum(mlist_i.^2);
      end
    end
    Mag_OpSpikeSqRatio = sqrt((Mag_opsq.^2)./(MagS.^2));
  end