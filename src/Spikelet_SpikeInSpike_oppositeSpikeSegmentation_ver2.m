function Info = Spikelet_SpikeInSpike_oppositeSpikeSegmentation_ver2(MagInfo,Param,SegDb_i,Names,D)

  %% parameters
  ParamIntra = Param.intra;
  Opp_NoiseThresRate = ParamIntra.Opp_NoiseThresRate; 
  Opp_ConsecuteMagRate = ParamIntra.Opp_ConsecuteMagRate;
  Opp_VibPVnum = ParamIntra.Opp_VibPVnum; 
  Opp_VibLengthRate = ParamIntra.Opp_VibLengthRate;
  
  %% original data (before z-normalization)
  from_org = SegDb_i(Names=="from");
  to_org = SegDb_i(Names=="to");
  D_org = D(from_org:to_org);

  %% SegStatistics (raw data)
  MatPos = Param.SegStatistics.positive2.MaT;
  MatNeg = Param.SegStatistics.negative2.MaT;
  MaT = min([abs(MatPos),abs(MatNeg)]);

  % Opp_Median_Ratio = 2;
  % Names = ["from","to","type","center","oppTerminal","parentSpikeId"];
  Names = ["from","to","type","center","parentSpikeId"];
  
  %% Oppsite spike in spike
  Mag = MagInfo.magnitude;
  [AbsMag_sorted, AbsMag_sorted_pos] = sort(abs(Mag),"descend");
  peak = AbsMag_sorted_pos(1);
  MagAbs_max = AbsMag_sorted(1);
  MagAbs_2nd = AbsMag_sorted(2);
  if MagAbs_2nd == 0
    Info = struct;
    Info.segmentList = [];
    Info.segmentList_names = Names;  
    return;
  end

  Left = MagInfo.left;
  Right = MagInfo.right;
  Supp = Right - Left + 1;
  Dnorm = MagInfo.data_org;
  FluxDnorm = max(Dnorm) - min(Dnorm);  

  ThresMatMin = MagAbs_2nd*Opp_ConsecuteMagRate;
  TimeMag = find(ThresMatMin <= abs(Mag));
  ValueMag = Mag(TimeMag);
  % ThresMatMax = MagAbs_max*Opp_MagRate(2);
  % TimeMag = find(ThresMatMin <= abs(Mag) & abs(Mag) <= ThresMatMax);
  SpikeIntervalThreshold = (length(Dnorm)*Opp_VibLengthRate)/Opp_VibPVnum;
  % FluctuationThreshold = 
  [CS,CS_sequence] = consecutive_spikes_PV(TimeMag,Opp_VibPVnum+1,SpikeIntervalThreshold);
  %{
   figure(); T = (1:length(Dnorm)); plot(T,Dnorm(T)); hold on;
             T1 = (CS(1,1):CS(1,end)); plot(T1,Dnorm(T1),'r','LineWidth',2);
  %}
  [CS,CS_sequence] = selectCS(CS_sequence,Mag,Opp_VibPVnum);

  SegmentList_oppSpike = zeros(size(CS,1),length(Names)); 
  for i = 1:size(CS,1)
      from_i = CS(i,1); to_i = CS(i,2);
      T_i = (from_i:to_i);
      [max_v,max_t] = max(D_org(T_i));
      [min_v,min_t] = min(D_org(T_i));
      mat_abs_i = max_v - min_v;
  
      Dnorm_i = Dnorm(T_i);
      CS_i = CS_sequence{i};
      DnormCs_i = Dnorm(CS_i);
      DnormCsDiff_i = diff(DnormCs_i);
      DnormCsDiffAbs = abs(DnormCsDiff_i);
      LegAbsMax_i = max(DnormCsDiffAbs);
      MeanFlux_i = mean(DnormCsDiffAbs);
      StdFlux_i = std(DnormCsDiffAbs);
      StdFluxRate_i = StdFlux_i/MeanFlux_i;
      DnormAprx_i = interp1(CS_i,Dnorm(CS_i),(T_i(1):T_i(end)));
      error_i = norm(Dnorm_i - DnormAprx_i')/norm(Dnorm_i-mean(Dnorm_i));
      %{
      
      lags = 2;
      [h, pValue, stat, cValue] = my_adftest(Dnorm_i - DnormAprx_i', lags);
      %}
      
      % LegAbsMax = mean(abs(diff(Dnorm(CS_i))));
      % if (LegAbsMax >= MagAbs_max*NoiseThresRate && error_i <= 0.3)
      % if LegAbsMax_i >= MagAbs_max*NoiseThresRate && (StdFluxRate_i <= 0.3 || error_i <= 0.3)
      ERR_Thr = 0.6; StdFlux_Thr = 0.6;
      % if LegAbsMax_i >= MagAbs_max*NoiseThresRate && error_i <= ERR_Thr && StdFluxRate_i <= StdFlux_Thr && (error_i <= ERR_Thr/2 || StdFluxRate_i <= StdFlux_Thr/2)
      if LegAbsMax_i >= MagAbs_max*Opp_NoiseThresRate
        Noise = false;  
      else
        Noise = true;
      end

      %{
      [maxNorm_v,maxNorm_t] = max(Dnorm(T_i));
      [minNorm_v,minNorm_t] = min(Dnorm(T_i));
      fluxNorm_i = maxNorm_v - minNorm_v;
      if fluxNorm_i < FluxDnorm*NoiseThresRate 
        Noise = true;
      else
        Noise = false;  
      end
      %}

      type_i = 0; 
      if min_t < max_t
        if mat_abs_i >= MatPos
          type_i = 1;
        end
      else
        if mat_abs_i >= abs(MatNeg)
          type_i = -1;  
        end   
      end
      center_i = floor((from_i+to_i)/2);
      SegmentList_oppSpike(i,:) =  [from_i,to_i,type_i,center_i,peak];
      if Noise
        SegmentList_oppSpike(i,1) = NaN;   
      else
        SegmentList_oppSpike(i,:) =  [from_i,to_i,type_i,center_i,peak];
      end
  end

  Index_NonNan = ~isnan(SegmentList_oppSpike(:,1));
  SegmentList_oppSpike = SegmentList_oppSpike(Index_NonNan,:);

  Info = struct;
  Info.segmentList = SegmentList_oppSpike;
  Info.segmentList_names = Names;
end

function [CS, CSseq_out] = selectCS(CSseq,Mag,Opp_VibPVnum)
  CS = NaN(length(CSseq),2);CSseq_out = cell(length(CSseq),1);
  for i = 1:length(CSseq)
    CSseq_i = CSseq{i};
    Mag_i = Mag(CSseq_i);
    max_sum_i = 0;
    for j = 1:length(CSseq_i)-Opp_VibPVnum+1
      Conseq_j = (j:j+Opp_VibPVnum-1);
      max_sum_ij = sum(abs(Mag_i(Conseq_j)));
      if max_sum_ij > max_sum_i
        CSseq_selected_i = CSseq_i(Conseq_j);
      end
    end
    CS(i,:) = [CSseq_selected_i(1),CSseq_selected_i(end)];
    CSseq_out{i} = CSseq_selected_i;
  end
end

function [CS,CS_seq] = consecutive_spikes_PV(TimeMag,Opp_VibPVnum,SpikeIntervalThreshold)
  CS = [];  CS_seq = {};
  if length(TimeMag) < Opp_VibPVnum
    return;  
  end

  cur_i = [TimeMag(1)];
  for i = 2:length(TimeMag)
    interval_i = TimeMag(i) - TimeMag(i-1);
    % fluctuation_i = ValueMag(i)- ValueMag(i-1);
    if interval_i > SpikeIntervalThreshold
      if length(cur_i) >=   Opp_VibPVnum
        CS = [CS; cur_i(1), cur_i(end)];  
        CS_seq = [CS_seq,cur_i];
      end
      cur_i = [TimeMag(i)];  
    else
      cur_i = [cur_i,TimeMag(i)];  
    end
    if i == length(TimeMag) && length(cur_i) >= Opp_VibPVnum
      CS = [CS; cur_i(1), cur_i(end)]; 
      CS_seq = [CS_seq,cur_i];
    end
  end
end


%% Option == "default",  
function to_constant =  findBoundaryFrom_rightTerminal(start,D,magSign,boundary_val,Option)
   if ~exist('Option','var')
     Option = "default";
   end

   P2T = (start:length(D));
   if magSign > 0
       to_cross = find( D(P2T) > boundary_val, 1,'last');
   else
       to_cross = find( D(P2T) < boundary_val, 1,'last');
   end
   
   if Option == "default"
     to_cross_pre = to_cross +1;
     [min_val, min_pos] = min([abs(D(to_cross)-boundary_val),abs(D(to_cross_pre)-boundary_val)]);
     if min_pos == 1
       to_constant = to_cross;  
     elseif min_pos == 2
       to_constant = to_cross_pre;
     end

   elseif Option == "monotone_extension"
     to_constant = to_cross;  
     if magSign > 0   
       for i = to_cross+1:length(D)
         if (D(i+1) - D(i)) < 0
           to_constant = i;
         else
           break;  
         end
       end
     elseif magSign < 0   
       for i = to_cross+1:length(D)
         if (D(i+1) - D(i)) > 0
           to_constant = i;
         else
           break;  
         end
       end
     end
   end
end

function [h, pValue, stat, cValue] = my_adftest(y, lags)
    % Augmented Dickey-Fuller (ADF) Test - Custom Implementation
    % 
    % INPUTS:
    %   y    - Time series data (vector)
    %   lags - Number of lags to include in the regression model
    % OUTPUTS:
    %   h      - Test result (0: Unit root exists / 1: Stationary)
    %   pValue - p-value of the test
    %   stat   - ADF test statistic
    %   cValue - Critical value at 5% significance level

    % Length of the data
    T = length(y);
    
    % Compute first-order differences
    dy = diff(y);
    
    % Construct lagged matrix X
    X = [y(1:end-1), ones(T-1,1)]; % Lagged y (including a constant term)
    
    if lags > 0
        for i = 1:lags
            X = [X, [zeros(i,1); dy(1:end-i)]]; % Add lagged differences
        end
    end
    
    % Define dependent variable (differenced series)
    Y = dy; 
    
    % Estimate coefficients using Ordinary Least Squares (OLS)
    beta = X \ Y; % OLS estimation
    residuals = Y - X * beta; % Compute residuals
    
    % Compute residual variance
    sigma2 = sum(residuals.^2) / (T - lags - 1);
    
    % Compute ADF test statistic
    stat = (beta(1) / sqrt(sigma2 / sum(X(:,1).^2))); 
    
    % Approximate 5% critical value for ADF test
    cValue = -2.86; % This value can be adjusted based on sample size
    
    % Approximate p-value (based on normal distribution)
    pValue = 1 - normcdf(stat);
    
    % Determine test result: reject H0 if stat < cValue
    h = (stat < cValue);
end