function Pub_GDFlex_decomposition_chap3(dataIdStr,Range,Command,MaT)
  %{
dataIdStr = '091'; Range = [113601,115000];
MaT = [1500, -2500];
Command = 'ScutMaTr-0p6_LScutMaTrFM-0p9';
Pub_GDFlex_decomposition_chap3(dataIdStr,Range,Command,MaT);
  %}

  %% Parameters for plotting
  Name = 'chapt3_decomp';
  %% directory
  TSAD2021_Dir = '../data';
  
  %%
  RawDataQuery = [TSAD2021_Dir,'/',dataIdStr,'*.mat']; 
  RawDataList = dir(RawDataQuery); 
  RawFile =  [TSAD2021_Dir,'/',RawDataList.name];
  RawData = load(RawFile);
  DataName = RawData.name;
  D = RawData.data;
  Anomaly_region =  RawData.anomaly_region;

  %% argument
  if ~exist('Command','var')
    Command = 'ScutMaTr-0p6_LScutMaTrFM-0p9';
  end
  if ~exist('Range','var')
    Range = [1,length(D)]; 
  end

  %% spike decomposition
  Arg = command2Arg(Command);
  ParamAprox = Arg2SpikeletParam(Arg);
  MagInfo = Spikelet_aproximation_ver_05(D, ParamAprox);
  MagInfo_initial = MagInfo.output.generateInitialSpikelet;

  %% main
  close all;
 
  fig_id = figure();
  Row = 2; Column = 2;
  t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'none');

  %subplot(Row,Column,1)
  nexttile(1);
   D = MagInfo.data_org; 
   PlotTimes = (Range(1):Range(2));
   plot(PlotTimes,D(PlotTimes))
   xlim([Range(1) Range(2)]); 
   axis off;
   
  % decomposition
  % subplot(Row,Column,2)
  nexttile(2);
   plot_decomposition(MagInfo,Range);
   axis off;

  %subplot(Row,Column,2+1)
  nexttile(2+1);
   Thres = MaT(1);
   plot_PN(MagInfo,Range,"positive",Thres);
   axis off;  
   
  % subplot(Row,Column,2+2)
  nexttile(2+2);
   Thres  = MaT(2);
   plot_PN(MagInfo,Range,"negative",Thres);
   axis off;

  % output
  %{
  FigFile = [OutputDir,'/',dataIdStr,'.fig'];
  saveas(fig_id,FigFile);
  FigFile = [OutputDir,'/',dataIdStr,'.svg'];
  saveas(fig_id,FigFile);
  %}
end

function  plot_PN(MagInfo,Range,PN,Thres)
  LineWidth = 1; 
  M = MagInfo.magnitude;  
  L = MagInfo.left;       
  R = MagInfo.right; 
  D = MagInfo.data_org;
  %% main
  if PN == "positive"
    P = find(M >= Thres);  
  else
    P = find(M <= Thres);  
  end
  index = (Range(1) <= P(:) &  P(:) <= Range(2));
  P = P(index);

  for i = 1:length(P)
    if rem(i,2) == 0
      Color = 'k'; 
    elseif rem(i,2) == 1
      Color = 'r';  
    end
    lt = L(P(i)); pt = P(i); rt = R(P(i)); 
    plot((lt:rt),D(lt:rt),Color, 'LineWidth', LineWidth);        
    if i == 1; hold on; end
    if i == length(P); hold off; end
  end
end

function Arg = command2Arg(Command)
  Arg = struct;
  Split = regexp(Command,'_','split');
  for i = 1:length(Split)
    Arg_i = Split{i};
    if length(Arg_i) >= 10 && strcmp(Arg_i(1:9),'ScutMaTr-')  
      Arg.SuppCutMaTr = str2num(strrep(Arg_i(10:end),'p','.'));   
    elseif length(Arg_i) >= 13 && strcmp(Arg_i(1:12),'LScutMaTrFM-')  
      Arg.LegSuppCutMaTrFM = str2num(strrep(Arg_i(13:end),'p','.'));  
    end
  end
end

function ParamAprx = Arg2SpikeletParam(Arg)
  % The order of the fields in if-statements is the same as the order of operations.
  Operation = struct;
  OS = {};  OSname = [];
  if isfield(Arg,"SuppCutMaTr")
    OS =  [OS,'restrictSupportByMagnitudeRatio'];  
    Operation.restrictSupportByMagnitudeRatio.magnitude_ratio = Arg.SuppCutMaTr;
    OSname = [OSname,'_RSMR'];
  end
  if isfield(Arg,"LegSuppCutMaTrFM")
    OS = [OS,'restrictSupportByLegMagnitudeRatioFirstMinimum'];   
    OSname = [OSname,'_RSLMRFM'];  
    Operation.restrictSupportByLegMagnitudeRatioFirstMinimum.magnitude_ratio = Arg.LegSuppCutMaTrFM;
  end
  Operation.operation_sequence = OS;
  Operation.operation_name = OSname;
  ParamAprx = struct;
  ParamAprx.operation = Operation;
end

function Range = plot_range(X, Ratio)
  max_x = max(X);
  min_x = min(X);
  
  Max = max(X) + Ratio*(max_x - min_x);
  Min = min(X) - Ratio*(max_x - min_x);
  Range = [Min, Max];
end

function plot_decomposition(MagInfo,Range)
  %% parameters
  Ratio = 0; 
  ExactOn = false; 
  ColorPos = 'k'; ColorNeg = 'b';
  Sign = 'all';  
  LineWidth = 1; 
  %%
  M = MagInfo.magnitude;  
  L = MagInfo.left;       
  R = MagInfo.right;  
  %% main
  P = find(M ~= 0);
  index = (Range(1) <= P(:) &  P(:) <= Range(2));
  P = P(index);
   
  min_v = inf; max_v = -inf;
  for i = 1:length(P)
     % if Debug_On && (i == 6 || i == 7); fprintf('[%s] debug\n', FuncName); end
     mv = M(P(i));
     lt = L(P(i)); pt = P(i); rt = R(P(i));
     
     if ExactOn
       lv = D(lt); pv = D(pt); rv = D(rt);  % spike values      
       LM = [(pv - lv), (pv - rv)];         % leg magnitude  
       [m_abs, base]  = min(sign(mv)*LM);
       if base == 1
         Spike = [lv, pv, rv] -lv;
       elseif base == 2
         Spike = [lv, pv, rv] - rv;
       else
         fprintf('[%s] unexpected minimum\n',FuncName)  
       end
     else
       Spike = [0, M(pt), 0]; 
     end
     
     Y = interp1([lt, pt, rt],Spike,(lt:rt)) ;
     if M(pt) >= 0
       Col = ColorPos;
     else
       Col = ColorNeg;
     end
     
     if strcmp(Sign,'positive') && M(pt) >= 0
       plot((lt:rt),Y,Col, 'LineWidth', LineWidth);
       min_v = min(min(Y), min_v); max_v = max(max(Y), max_v);
     elseif strcmp(Sign,'negative') && M(pt) <= 0
       plot((lt:rt),Y,Col, 'LineWidth', LineWidth);  
       min_v = min(min(Y), min_v); max_v = max(max(Y), max_v);
     elseif strcmp(Sign,'all') 
       plot((lt:rt),Y,Col, 'LineWidth', LineWidth);        
       min_v = min(min(Y), min_v); max_v = max(max(Y), max_v);
     end  
     if i == 1; hold on; end
     if i == length(P); hold off; end
  end
   
  xlim([Range(1) Range(2)]); 
  ylim(plot_range([min_v, max_v], Ratio));
end