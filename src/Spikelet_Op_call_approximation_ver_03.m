function  MagInfo = Spikelet_Op_call_approximation_ver_03(MagInfo)
  InitialOpName = 'generateInitialSpikelet';
  % PLOT = true; 
  %% call operation sequence
  Ops = MagInfo.param.operation.operation_sequence;
  for i = 1:length(Ops)
    Op = Ops{i};
    MagInfo.current_operation = i;
    tic;  
   %% 1. Spike reduction
   if strcmp(Op,'reduceSpikeByMagnitude')                            % 5trial new
     MagInfo = Spikelet_Op_reduceSpikeByMagnitude_ver_02(MagInfo); 
   elseif strcmp(Op,'reduceSpikeByInclusion')                        % 5trial new
     MagInfo = Spikelet_Op_reduceSpikeByInclusion(MagInfo);  
   elseif strcmp(Op,'reduceSpikeByQueryPSR')                         % 5trial new % PSR: Percentile Support Ratio
     MagInfo = Spikelet_Op_reduceSpikeByQueryPSR_ver02(MagInfo);     
   elseif strcmp(Op,'reduceSpikeByMagnitudeRatio') 
     MagInfo = Spikelet_reduceSpikeByMagnitudeRatio(MagInfo);  
   elseif strcmp(Op,'reduceSpikeByQuery') 
     MagInfo = Spikelet_reduceSpikeByQuery(MagInfo);  
       
   %% 2. Support restriction/extension
   elseif strcmp(Op,'restrictSupportByLegMagnitudeRatio')         % 5trial new
     MagInfo = Spikelet_Op_restrictSupportByLegMagnitudeRatio(MagInfo);
   elseif strcmp(Op,'restrictSupportByLegMagnitudeRatioFirstMinimum')         % 5trial new
     MagInfo = Spikelet_Op_restrictSupportByLegMagnitudeRatioFirstMinimum(MagInfo);  
   elseif strcmp(Op,'restrictSupportByLegSupportLength')          % 5trial new
     MagInfo = Spikelet_Op_restrictSupportByLegSupportLength(MagInfo);
   elseif strcmp(Op,'restrictSupportBySupportLength')              % 5trial new
    MagInfo = Spikelet_Op_restrictSupportBySupportLength(MagInfo);
   elseif strcmp(Op,'restrictSupportByMagnitudeRatio')            % 5trial new
    MagInfo = Spikelet_Op_restrictSupportByMagnitudeRatio(MagInfo);
   elseif strcmp(Op,'extendLegSupportMaximal')                    % 5trial new
    MagInfo = Spikelet_Op_extendLegSupportMaximal(MagInfo);
    
   elseif strcmp(Op,'restrictSupportByWindowLength')
     MagInfo = Spikelet_restrictSupportByWindowLength(MagInfo); 
   elseif strcmp(Op,'restrictSupportByMagnitudeRatioInitial')
     MagInfo = Spikelet_restrictSupportByMagnitudeRatioInitial(MagInfo); 
   elseif strcmp(Op,'restrictSupportBySolidarity')
     MagInfo = Spikelet_Op_restrictSupportBySolidarity(MagInfo);        
   elseif strcmp(Op,'restrictSupportByReferenceValue')
     MagInfo = Spikelet_Op_restrictSupportByReferenceValue(MagInfo);   
    %%   
  
  %% 3. Constant segment extraction
   elseif strcmp(Op,'extractConstantInLeg')
     MagInfo = Spikelet_Op_extractConstantInLeg(MagInfo);   % 5trial new
   elseif strcmp(Op,'extractConstantSegment')
     % MagInfo = Spikelet_Op_extractConstantSegment(MagInfo);   % correct version
     MagInfo = Spikelet_Op_extractConstantSegment_ver02(MagInfo); % latest version is 03
     CLThr = MagInfo.output.(Op).length_threshold;
     fprintf('%s: Constant Length Threshold = %d\n',Op,CLThr);
   elseif strcmp(Op,'findConstantInSpike')
     MagInfo = Spikelet_Op_findConstantInSpike(MagInfo); 
     CLThr = MagInfo.output.(Op).length_threshold;
     fprintf('%s: Constant Length Threshold = %d\n',Op,CLThr);       
   else
     fprintf('unknown operation (%s)\n', Op)
   end
   
   %% time and size
   MagInfo.time.(Op) = toc;
   spike_num = sum(MagInfo.output.(Op).magnitude ~= 0);
   if isfield(MagInfo.output.(Op),'constant')
     constant_num =  size(MagInfo.output.(Op).constant,1);
   else
     constant_num = 0;  
   end
   MagInfo.size.(Op) = spike_num + constant_num;
   
   % message
   fprintf([Op,': time = ',       num2str(MagInfo.time.(Op)),            ', size= ',num2str(MagInfo.size.(Op)),'\n']); 
  end
  
  %% total tile 
  Total_time = MagInfo.time.(InitialOpName);
  for i = 1:length(Ops)
    Op = Ops{i};   
    Total_time = Total_time + MagInfo.time.(Op);
  end
  MagInfo.time.operation_all = Total_time;
  MagInfo.size.operation_all =  MagInfo.size.(Ops{end});

  %% plot
  if exist('PLOT','var') && PLOT
     D  = MagInfo.data_org; 
     X_Range =  [1 length(D)];  X = (1:length(D)); LineWidth = 0.75; Sign = 'positive';
     % Rate = 1.05;  Shift = 1.0;
     % Y_LIM = [min(D)*Rate,max(D)*Rate];  Y_Width = Y_LIM(2) - Y_LIM(1); 
     % Y_LIM_P =  [0-Shift,Y_Width-Shift]; 

     figure(); Row = 3;
     subplot(Row,1,1)
       plot(X,D); xlim(X_Range); % ylim(Y_LIM_P);
     subplot(Row,1,2)
       Op_PLOT =MagInfo.param.operation.operation_sequence{end-1};
       Spikelet_plot_decomposition_ver_03(MagInfo,Op_PLOT,X_Range,LineWidth,Sign); % ylim(Y_LIM_P);
       title(Op_PLOT);
     subplot(Row,1,3)
       Op_PLOT =MagInfo.param.operation.operation_sequence{end};
       Spikelet_plot_decomposition_ver_03(MagInfo,Op_PLOT,X_Range,LineWidth,Sign); % ylim(Y_LIM_P);
       title(Op_PLOT);
  end
end