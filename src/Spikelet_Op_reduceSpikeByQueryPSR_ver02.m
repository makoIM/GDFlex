function MagInfo = Spikelet_Op_reduceSpikeByQueryPSR_ver02(MagInfo) 
  FuncName = 'Spikelet_Op_reduceSpikeByQueryPSR_ver02';
  OpName = 'reduceSpikeByQueryPSR';
  
  %% parameters
  ParamPSR = MagInfo.param.operation.(OpName);
  QueryPSR = ParamPSR.query;
  QueryName = ParamPSR.query_name;
  
  %% find center
  [SpikeCenter, ReducedTime, FeatureMatrix, FeatureMatrix_Names, MsqV]= Spikelet_eval_queryPSR(QueryPSR, MagInfo);
  
  %% revise spikelet decomposition
  
  % magnitude
  if ~isempty(ReducedTime)
    MagInfo.magnitude(ReducedTime) = zeros(length(ReducedTime),1);    
    MagInfo.value(ReducedTime) = nan(length(ReducedTime),1);
    MagInfo.type(ReducedTime) = nan(length(ReducedTime),1);
    MagInfo.left(ReducedTime) = ReducedTime;
    MagInfo.right(ReducedTime) = ReducedTime;
  end
  PSRqueryInfo = struct;
  PSRqueryInfo.Feature = FeatureMatrix;
  PSRqueryInfo.Feature_Names = FeatureMatrix_Names;
  PSRqueryInfo.MsqV = MsqV;

  %% output
  MagInfo.param_str = [QueryName, ' ',QueryPSR];
  MagInfo = Spikelet_MagInfo_post_processing(MagInfo, OpName);
  MagInfo.output.(OpName).PSRqueryInfo = PSRqueryInfo;  % Operation depended field

  %% PLOT
  if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
  end   
end

function  [SpikeCenter, ReducedTime,Matrix_true,Matrix_Names,MsqV]= Spikelet_eval_queryPSR(QueryPSR, MagInfo)
  % PlOT = true;

  MsqV = Spikelet_MSQvariables_ver03(MagInfo,QueryPSR);
  Names = MsqV.Names;  Matrix_Names = Names;
  
  % Data_org = MsqV.meta.data_org;

  %% Matrix
  Peak_time = MsqV.Matrix(:,Names=="peak_time");
  % mag 4/5
  Mag = MsqV.Matrix(:,Names=="spike_magnitude");
  MagR = MsqV.Matrix(:,Names=="spike_magnitude_right");
  MagL = MsqV.Matrix(:,Names=="spike_magnitude_left");
  MagLRratio = MsqV.Matrix(:,Names=="spike_magnitude_lrRatio");
  MagOpRatio = MsqV.Matrix(:,Names=="spike_magnitude_OpSpikeSqRatio");

  % support 4/9
  Supp = MsqV.Matrix(:,Names=="spike_support");
  SuppR = MsqV.Matrix(:,Names=="spike_support_right");
  SuppL = MsqV.Matrix(:,Names=="spike_support_left");
  SuppLRratio = MsqV.Matrix(:,Names=="spike_support_lrRatio");
  % val/LR 5/14
  Peak_val = MsqV.Matrix(:,Names=="peak_val");
  Left_val = MsqV.Matrix(:,Names=="left_val");
  Right_val = MsqV.Matrix(:,Names=="right_val");
  Left =  MsqV.Matrix(:,Names=="left");
  Right =  MsqV.Matrix(:,Names=="right");

  % Qrario 6/20
  Qratio12_len = MsqV.Matrix(:,Names=="quantile_length_both_ratio_12");
  Qratio13_len = MsqV.Matrix(:,Names=="quantile_length_both_ratio_13");
  Qratio14_len = MsqV.Matrix(:,Names=="quantile_length_both_ratio_14");
  Qratio23_len = MsqV.Matrix(:,Names=="quantile_length_both_ratio_23");
  Qratio24_len = MsqV.Matrix(:,Names=="quantile_length_both_ratio_24");
  Qratio34_len = MsqV.Matrix(:,Names=="quantile_length_both_ratio_34");
 
  % QL  = MsqV.Matr23ix(:,Names=="quantile_length_both");
  % QLT  = MsqV.Matrix(:,Names=="quantile_length_both_total");
  % QLR = MsqV.Matrix(:,Names=="quantile_length_right");
  % QLRT = MsqV.Matrix(:,Names=="quantile_length_right_total");
  % QLL = MsqV.Matrix(:,Names=="quantile_length_left");
  % QLLT = MsqV.Matrix(:,Names=="quantile_length_left_total");

  if ~isempty(QueryPSR) && contains(QueryPSR,'Mag_OpSpike')
    Mag_OpSpike = MsqV.Matrix(:,Names=="spike_magnitude_OpSpikeSqRatio");
  end
  %% evaluate
  if isempty(QueryPSR)
   Q_index = true(length(Peak_time),1);
  else
   Q_index = eval(QueryPSR);
  end
  SpikeCenter = Peak_time(Q_index);
  Matrix_true = MsqV.Matrix(Q_index,:);
  ReducedTime= setdiff(MagInfo.center, SpikeCenter);

  %% plot
  if exist('PLOT','var') && PLOT
     D  = MagInfo.data; 
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
