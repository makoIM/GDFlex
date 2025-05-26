function  MagInfo = Spikelet_generateInitialSpikelet_ver_02(X, RightWindowLength)  
%{
 ver.02: Intorucing the data-structure for Spikelet Decomposition.
[Input]
X time series (1 row, N columns)
RightWindowLength: forward window length (the length to the next index is 2)
[Output]
Mag: Magnitude
L:   Support length
C: [sign,left,right,balanced_left,balanced_right] % sign:1 convex, -1:convex
   sign; 1 convex, -1:convex, 2, -2: terminal continuation (trapezoidal maximum/minimum value)C3,-3=terminal 
%}
DEBUG = true; DEBUG = false;
  %% argument
  if nargin == 1
    RightWindowLength = Inf; 
  end
  
  %% Function Declaration
  tic
  OpName = 'generateInitialSpikelet';
  FuncName = 'Spikelet_generateInitialSpikelet';
  
  Debug_time = 3048; Debug_time = 0;
  
  % 0. 1 Global Declaratoin
  global  LEG_WINDOW_LENGTH_ON WINDOW_LENGTH;
  global  LEG_MAXIMAL_CONVEX_AMPLITUDE_2_SUPPORT_ON;
  global  DEBUG_FIND_LR;
  if ~isempty(LEG_MAXIMAL_CONVEX_AMPLITUDE_2_SUPPORT_ON)
     LEG_MAXIMAL_CONVEX_AMPLITUDE_2_SUPPORT_ON = true;
  end
  
  DEBUG_FIND_LR = 0; % DEBUG_FIND_LR = 1; 
  Debug_Input =0; % Debug_Input = 1; 
  Debug_GetNext = 0; % Debug_GetNext = 1;    
  Debug_ForcedReduction = 0; % Debug_ForcedReduction = 1;
  Debug_MiddleReduction = 0; % Debug_MiddleReduction = 1; 
  Debug_LeftReduction = 0; % Debug_LeftReduction = 1; 
  Debug_RightReduction = 0; % Debug_RightReduction = 1;
  
  %% 0.3 Check Input

  % 1st Argument
  if size(X,2) > size(X,1)
    X = X';  
  end
  X = X(:,1);   if Debug_Input == 1; disp(X'); end

  % 0.3 Set local variables
  
  %------------
  % 0. Initialization
  %------------ 
  n = length(X);   % the length of an input time series 
  A = zeros(n,1);  % Amplitude
  L = zeros(n,1);  % Support length
  SubseqId = (1:n)';  
  C = repmat(SubseqId,1,5); % Support of convex [sign,left,right,valanced_left,valanced_right] % sign:1 convex, -1:convex
  S= Stack_create(length(X)); % stack
  m = NaN; % Previous local maximaltime
  s = NaN(1,4); % Local stacks
  %------------
  % 1. Main loop
  %------------  
  Apl = false; 
  i = 1;  S = Stack_push(S, 1); 
  while ~isnan(i) && (i <= n-1 || (i == n && Apl == true && S.p >= 3 )) % When all the value had been read and no reduction rules are appliable
       if Debug_GetNext && i == 12; fprintf('Dubug i = %f\n', i); end
      
    %% (1) get next vertex   
    if i <= n-1 && Apl == false    % Read next maximum   
       % (1.1a) get next Local maximal transformation
       [next,C] = getNextLocalMaximum(X, C, i);  
       
       if next == Debug_time
         fprintf('debug\n');   
       end
       if Debug_GetNext == 1;  if next <= length(X);  fprintf('read: X(%d)=', next);  disp(X(next));  else;  fprintf('read: completed\n'); end; end;

       % (1.2) Window Length Limit check % The bottom time in stack S is earlier than (i - window_length + 1)
       while  (S.p >= 2) && ( S.v(2) < (next - RightWindowLength + 1)) 
          right_end = S.v(2) +  RightWindowLength-1;
          [S, A, C, L] = forced_left_reduction(X, S, A, C, L, right_end);
          if Debug_ForcedReduction == true; forced_left_reduction_debug_write(S,i,next); end;  
       end
       
       % (1.1b) post processing of getNextLocalMaximum
       i = next;
       if ~isnan(next) && next <= length(X)
          S = Stack_push(S, next); % push the next local maximum
       end
    end
    if Debug_GetNext == 1; fprintf("Stack = ");fprintf("%d,",X(S.v(1:S.p))); fprintf("\n");end 
    
    %% (2) Leg reductions
    if S.p >= 4     % (2.1) Middle leg reduction
      p = S.p;   s1=S.v(p - 3);  s2=S.v(p - 2); s3=S.v(p - 1); s4=S.v(p);  % s1,s2,s3,s4 ordered by time (opposite to the order in the paper)
      if abs(X(s4) - X(s3)) >= abs(X(s3) - X(s2)) ...
        && abs(X(s2) - X(s1)) >  abs(X(s3) - X(s2))        
        if Debug_MiddleReduction == true; fprintf("Middle leg reduction (%d,%d,%d,%d)\n",s1,s2,s3,s4); end
        if s2 == Debug_time || s3 == Debug_time
           fprintf('[%s] debug \n',FuncName);  
           [s1,s2,s3,s4]
           figure(); plot([s1:s4],X(s1:s4)); 
           xlim([s1 s4]);
        end;

        [S, A, C, L] = middle_leg_reduction(X,S,A,C,L,s1,s2,s3,s4);
        Apl1 = true;
      else
        Apl1 = false;
      end
    else
      Apl1 = false;
    end
    if S.p >= 3  %(2.2) Left leg reduction
      s1 = S.v(1); s2 = S.v(2); s3 = S.v(3); 
      if abs(X(s3) - X(s2))  >= abs(X(s2) - X(s1))    
        if Debug_LeftReduction == true; fprintf("Left leg reduction (%d,%d,%d)\n",s1,s2,s3); end    
        if s2 == Debug_time
           fprintf('[%s] debug \n',FuncName);  
           [s1,s2,s3]
           figure(); plot([s1:s3],X(s1:s3)); 
           xlim([s1 s3]);
        end
        
        [S, A, C, L] = left_leg_reduction(X,S,A,C,L,s1,s2,s3);
        Apl2 = true;
      else
       Apl2 = false;  
      end
    else
      Apl2 = false;  
    end
    
    Apl = Apl1 + Apl2;
  end
    
    % (3) right leg reduction (Reversed terminal leg reduction)
    while S.p >= 3
      p = S.p;   s1= S.v(p-2); s2= S.v(p-1); s3= S.v(p);
      
      if Debug_RightReduction == true; fprintf("Right leg reduction (%d,%d,%d)\n",s1,s2,s3); end;     
      if s2 == Debug_time
           fprintf('[%s] debug \n',FuncName);  
           [s1,s2,s3]
           figure();plot([s1:s3],X(s1:s3)); 
           xlim([s1 s3]);
      end
      
      [S, A, C, L] = right_leg_reduction(X,S,A,C,L,s1,s2,s3);
      
      %X(S.v(1:S.p))
    end
    
    % time series representation of initial Spikelet decomposition
    initial_mag_time = find(A ~= 0);
    initial_mag_trend = zeros(length(A),1); 
    initial_mag_trend(initial_mag_time) = A(initial_mag_time);
    if ~isempty(initial_mag_time)
      if initial_mag_time(1) ~= 1
        initial_mag_time = [1;initial_mag_time];  
      end
      if initial_mag_time(end) ~= length(X)
        initial_mag_time = [initial_mag_time;length(X)];    
      end
      Data_initial = interp1(initial_mag_time,X(initial_mag_time),[1:length(X)]);
      Data_initial = Data_initial';
    else
      Data_initial = initial_mag_trend;  
    end
  
   %% set Spikelet Decomposition
    % set original
    MagInfo = struct;
    MagInfo.fields = {'center','type','value','left','right'...
                     ,'data','magnitude','leg_magnitude'}; 
                 
    MagInfo.data_org = X;
    MagInfo.leg_magnitude = zeros(length(X),1);
    
    % initial spikelet
    Len = length(MagInfo.data_org);
    MagInfo.data = Data_initial;
    MagInfo.magnitude = A;
    MagInfo.center = find(A ~= 0);
    
    MagInfo.type = NaN(Len,1);
    MagInfo.value = NaN(size(MagInfo.data));
    if ~isempty(MagInfo.center)
      MagInfo.type(MagInfo.center) = 2 * sign(A(MagInfo.center));
      MagInfo.value(MagInfo.center) = X(MagInfo.center);
    else
      MagInfo.type = NaN(Len,1);
      MagInfo.value= NaN(Len,1);
    end
    MagInfo.left  = C(:,2);
    MagInfo.right = C(:,3);
    
    %% save_initial spikelet
    MagInfo.output =[];
    MagInfo = Spikelet_MagInfo_post_processing(MagInfo, OpName);
  if DEBUG
    [SpikeDb, SpNames] = Spikelet_MagInfo_extractSpikeDB(MagInfo);   
  end
end

%% A: Find left and right terminal of a convex
%==========================================
% search from terminals to the peak to a spike
%==========================================
% [Backward] Serach the first time-index whose value exceeds to_short for the first time from 'from' to 'to'
% If there are not such a time,return the nearest value of 'to'

function b_left = Leg_find_left(X,from,to,to_short)
  global DEBUG_FIND_LR;
  FuncName = 'Leg_find_left';
  if DEBUG_FIND_LR == 1
    fprintf('[%s] (%d,%d,%d)\n', FuncName, from, to, to_short); 
  end
  
  % 0. 1 Global Declaratoin
  global  LEG_MAXIMAL_CONVEX_AMPLITUDE_2__SUPPORT_ON; 

  if LEG_MAXIMAL_CONVEX_AMPLITUDE_2__SUPPORT_ON == true
   diff_ft = X(to_short) - X(from);    
   if diff_ft > 0      % increasing towawrd left
     cand_rel = find(X(to:from-1) >= X(to_short));
   elseif diff_ft < 0  % decreasing toward left
     cand_rel = find(X(to:from-1) <= X(to_short)); 
   else
     fprintf('[%s] diff_ft == 0 (from, to, to_short)=(%d, %d, %d)', FuncName, from,to,to_short);  
   end
   if isempty(cand_rel)
       if diff_ft > 0
          [val_max, cand_rel] =  max(X(to:from-1));
       elseif diff_ft <0
          [val_min, cand_rel] = min(X(to:from-1));
       end    
   end
     from_to_right_rel = max(cand_rel);   
     b_left = to + from_to_right_rel -1; 
  else 
    b_left = to;
  end
end

% [Forward] Serach the first time-index whose value exceeds to_short for the first time from 'from' to 'to' 
function b_right = Leg_find_right(X,from,to,to_short) 
  global DEBUG_FIND_LR;
  FuncName = 'Leg_find_right';
  if DEBUG_FIND_LR == 1
     fprintf('[%s] (%d,%d,%d)\n', FuncName, from, to, to_short); 
  end 
  
  % 0. 1 Global Declaratoin
  global  LEG_MAXIMAL_CONVEX_AMPLITUDE_2_SUPPORT_ON; 

  if LEG_MAXIMAL_CONVEX_AMPLITUDE_2_SUPPORT_ON == true
   diff_ft = X(from) - X(to_short);    
   if diff_ft > 0      % decreasing toward right
     cand_rel = find(X(from+1:to) <= X(to_short));
   elseif diff_ft < 0  % increasing toward right
     cand_rel = find(X(from+1:to) >= X(to_short)); 
   else
     fprintf('[%s] diff_ft == 0 (from, to, to_short)=(%d, %d, %d)', FuncName, from,to,to_short);  
   end
   if isempty(cand_rel)
       if diff_ft > 0
          [val_min, cand_rel] =  min(X(from+1:to));
       elseif diff_ft <0
          [val_max, cand_rel] = max(X(from+1:to));
       end
   end
        
   from_to_right_rel = min(cand_rel);
   b_right = from + from_to_right_rel; 
   
  else 
    b_right = to;
  end
end

%===================
% B: Reduction
%===================
% (1) Two-side reduction s1,s2,s3,s4 is in time-order [s1,s2,s3],[s2,s3,s4] are spikes.
% Each center is s2 and s3, respectively.
function [S, A, C, L] = middle_leg_reduction(X,S,A,C,L,s1,s2,s3,s4)
       global  LEG_WINDOW_LENGTH_ON WINDOW_LENGTH;
       FuncName = 'middle_leg_reduction';

       %% (1) spike s3
        % (1.1) find left and right in the convex s3
        s3_left = s2; 
        if LEG_WINDOW_LENGTH_ON 
           s2_limit = s3-WINDOW_LENGTH+1;
           if  s2< s2_limit
               s3_left = Leg_find_left(X,s3,s2_limit,s4);
           end
        end
        s3_right = Leg_find_right(X,s3,s4,s2);

       % (1.2) calculate the convex s3
       s3_sign =  sign((X(s3) - X(s2)));      
       if s3_sign ~= sign(C(s3,1)); fprintf('1:The sign of local maximum (%d) is mismached.\n',s3); end; 
        C(s3,:) = [s3_sign, s2, s4, s3_left, s3_right]; 
        if LEG_WINDOW_LENGTH_ON
          A(s3) = s3_sign*min(abs(X(s3) - X(s3_left)), abs(X(s3) - X(s3_right))); 
        else
          A(s3) = X(s3) - X(s2);  
        end    
        L(s3) = s3_right - s3_left + 1;
  
       %% (2) spike s2
        %if s2 == 508; fprintf('[%s]\n',FuncName); end % debug
        
        % (2.1) find left and right in the convex s2
        s2_right = s3; s1_restricted = s1; 
        if LEG_WINDOW_LENGTH_ON 
           s1_limit = s2-WINDOW_LENGTH+1;
           if  s1< s1_limit
             s1_restricted = Leg_find_left(X,s2,s1_limit,s3);
           end
        end
        s2_left = Leg_find_left(X,s2,s1_restricted,s3);
        
        % (2.2) calculate the convex s2
        s2_sign =  sign(X(s2) - X(s1)); 
        if s2_sign ~= sign(-1*s3_sign); fprintf('2:The sign of local maximum (%d) is mismached to that of %d.\n',s2, s3);  end;
        C(s2,:) = [s2_sign, s1, s3, s2_left, s2_right];  % 
        if LEG_WINDOW_LENGTH_ON
           A(s2) = s2_sign*min(abs(X(s2) - X(s2_left)), abs(X(s2) - X(s2_right))); 
        else
          A(s2) = X(s2) - X(s3);  
        end
        L(s2) = s2_right - s2_left + 1;
        
       %% (3) pop stack for s2 and s3
        S = Stack_pop_2(S);     
        S = Stack_pop_2(S);
end

%% (2.a) left reduction s1,s2,s3 is in time-order. (s1 is left terminal)@ 
%         [s1,s2,s3] is a spikeCeach center is s2
function [S, A, C, L] = left_leg_reduction(X,S,A,C,L,s1,s2,s3)
        global  LEG_WINDOW_LENGTH_ON WINDOW_LENGTH;     
        %if s2 == 508; fprintf('[%s]\n',FuncName); end % debug
        
        % (1) find left anb right in the convex
        s2_left = s1;
        s2_right = s3;
        if LEG_WINDOW_LENGTH_ON 
           s1_limit = s2 - WINDOW_LENGTH+1;
           if  s1< s1_limit
             s1_restricted = Leg_find_left(X,s2,s1_limit,s3);
             s2_left = Leg_find_left(X,s2,s1_restricted,s3);
           end
        end
        
        % (2) calculate convex
        s2_sign =  sign((X(s2) - X(s1)));   %s2_sign =  sign((X(s2) - X(s2_left)));
        if s2_sign ~= sign(C(s2,1));fprintf('3:The sign of local maximum (%d) is mismached.\n',s2); end    
        C(s2,:) = [s2_sign, s1, s3, s2_left, s2_right]; 
        if LEG_WINDOW_LENGTH_ON == true
           A(s2) = s2_sign*min(abs(X(s2) - X(s2_left)), abs(X(s2) - X(s2_right))); 
        else
           A(s2) = X(s2) - X(s1);  % the magnitude of a left leg is smaller than that of a right leg
        end
        L(s2) = s2_right - s2_left + 1;
       % if LEG_WINDOW_LENGTH_ON && L(s2) > 200; fprintf('[%s]\d',FuncName); end % debug
       
       % (3) pop stack
       S = Stack_pop_bottom(S); 
end

%% (2.b) Forcefully reduce left reduction because the window length of a left leg is over RightWindowLength
% This algoritm is not symmetric because it is online-algorithm. 
function [S, A, C, L] = forced_left_reduction(X, S, A, C, L, i)
      global  LEG_WINDOW_LENGTH_ON WINDOW_LENGTH;
      FuncName = 'forced_left_reduction';
      % if s2 == 8; fprintf('[%s] s2 = %d\n',FuncName, s2); end % for debug      
      if LEG_WINDOW_LENGTH_ON ~= true
       fprint('[%s] Unexpected LEG_WINDOW_LENGTH_ON = false. i = %d \n', i);   
      end
      
      % (1) set local stack
      if S.p >= 2
          s1 = S.v(1);
          s2 = S.v(2);

         % (2) find left and right of a convex
         if S.p >= 3
             s3 = S.v(3);   % S.v(3) is the right most of reduction sequence.
          else
             s3 = Leg_find_right(X, s2, i, s1); 
         end
         s2_right = s3;
         s2_left = s1;
         if LEG_WINDOW_LENGTH_ON 
             s1_restricted = s2-WINDOW_LENGTH+1;
          if  s1< s1_restricted
             s2_left = Leg_find_left(X,s2,s1_restricted,s3);
          end
         end

          % (3) calculate convex
          c_sign =  sign((X(s2) - X(s3)));
          if c_sign ~= sign(C(s2,1)); fprintf('3:The sign of local maximum (%d) is mismached.\n',s2); end  
          A(s2) =c_sign*min(abs(X(s2) - X(s2_left)), abs(X(s2) - X(s2_right)));
          C(s2,:) = [c_sign, s1, s3, s2_left, s2_right];   
          L(s2) = s2_right - s2_left + 1;
          
          % (3) pop stack:  remove the left most time s(1) of the local stack s(1:3)    
          [S,val] = Stack_pop_bottom(S);   
      else
          S = Stack_replace_bottom(S, S.v(1)+1);    
      end
end

function forced_left_reduction_debug_write(S,i,next)
  s1 = S.v(1); s2 = S.v(2);  
  if S.p >= 3
    s3 = S.v(3); 
    fprintf('Forced left leg reduction. S.p>=3. (%d,%d,%d). next=%d ',s1,s2,s3);
  elseif S.p == 2
    fprintf('Forced left leg reduction. S.p=2. (%d,%d,%d) ',s1,s2,i);
  end
  fprintf('next=%d \n',next); 
end


%% (3) right reduction. s1,s2,s3 are in time-order. s3 is the right terminal. [s1,s2,s3] is a spike. its center is s2.
function [S, A, C, L] = right_leg_reduction(X,S,A,C,L,s1,s2,s3)
      global  LEG_WINDOW_LENGTH_ON WINDOW_LENGTH;  
      FuncName = 'right_leg_reduction';
      
      % (1) find left anb right in the convex
      s2_left = s1;
      s2_right = s3;
      if LEG_WINDOW_LENGTH_ON 
           s1_limit = s1-WINDOW_LENGTH+1;
           if  s1< s1_limit
             s1_restricted = Leg_find_left(X,s2,s1_limit,s3);
             s2_left = Leg_find_left(X,s2,s1_restricted,s3);
           end
      end
      
      % (2) calculate convex
      s2_sign =  sign((X(s2) - X(s1))); % s2_sign =  sign((X(s2) - X(s2_left)));
      if s2_sign ~= sign(C(s2,1));
          fprintf('4:The sign of local maximum (%d) is mismached.\n',s2);
      end  
      
      if LEG_WINDOW_LENGTH_ON
         A(s2) = s2_sign*min(abs(X(s2) - X(s2_left)), abs(X(s2) - X(s2_right))); 
      else
         A(s2) = (X(s2) - X(s3));   % the magnitude of a left leg is smaller tham that of a right one.
      end
      C(s2,:) = [s2_sign, s1, s3, s2_left, s2_right];
      L(s2) = s2_right - s2_left + 1;
      %{
      if LEG_WINDOW_LENGTH_ON &&  L(s2) > 200
         fprintf('[%s]\d', FuncName); 
      end
      %}
      % pop
      S= Stack_pop_1(S);
end

%% C. Get Next vertex
% ===================
%  getNextLocalMaximum: find next local minimam/maximam.
% ===================
function [next,C] = getNextLocalMaximum(X, C, present)
 if present >= length(X)   
    next = NaN;
    return;
 end
 
 % (1) find next time-index whose value is different. = next_diff ( while i has a minimum/maximum value.
 
 i = present + 1;
 if present == 0
   present_val = X(i);  
   present_sign = 0;  % unknown at present
 else
   present_val = X(present);  
   present_sign = sign(C(present,1)); 
 end
 
 while X(i) == present_val
    C(i,1) = 2*present_sign;
    if i == length(X)
      next = NaN;  
      return;  
    else    
      i = i+1;  % next value that has a different value from that of a present time.
    end  
 end
 next_diff= i;  % next value that is not a local minimum/maximum.
 
 % (2) When next_diff is end-index, exit. 
 if next_diff >= length(X)
    C(i,1) = -3*C(present,1);   
    next = length(X);  
     return; 
 end
 
 % (3) When next_diff is not an end-incex, find an local minimum/maximum.
 current_slope = X(next_diff) - present_val; 
 while current_slope*(X(i+1)-X(i)) >= 0 
    C(i,1) = 0; 
    if i== length(X)
       C(i,1) = -3*present_sign;   
       return;
    else
       i = i+1;   
       % (3-1) When i is end-index
       if i >= length(X)
          C(i,1) = -3*present_sign;   
          next = length(X);  
          return; 
       end
    end  
 end
    
 % (3-2) When i is not end-index
   next = i;  % The next to the present one is local maximal or minimum.
   C(i,1) = sign((X(i)-X(i+1)));
   if present == 0
      C(1,1) =  -3*C(i,1);
      C(2:i-1,1) = 2*sign(C(1,1));
   end
   return;
 end

%% D. Stack Management
% ===================
%  The operations on a stack
% ===================
function Stack_plot(S)
  S.v(1:S.p)
end

function S = Stack_create(size)
  S = struct;
  S.v = NaN(size, 1);
  S.p = 0;
  S.size = size;
end

% =================== 
% Push dt to stac S
% ===================
function S = Stack_push(S, val)
   if( S.p >= S.size ) 
     S.p = Inf;
     fprintf('Push failed because of Stack overflow.\n');
     return;
   end
   
   S.p = S.p + 1;
   S.v(S.p) = val;
end

% ================================= 
% Pop a value at the top of S
% =================================
function [S,val] = Stack_pop_1(S)

    if( S.p <= 0 ) 
      fprintf('Pop failed because of Empty Stack.\n');
    end
    
    val = S.v(S.p); 
    S.p = S.p - 1;
end

% ================================= 
% Pop the second value from the top
% =================================
function [S,val] = Stack_pop_2(S)

    if( S.p <= 0 ) 
      fprintf('Pop failed because of Empty Stack.\n');
    end
    
    val = S.v(S.p - 1); 
    S.v(S.p - 1) = S.v(S.p);
    S.p = S.p - 1;
end

% ================================= 
% Pop the third value from the top
% =================================
%{
function [S,val] = Stack_pop_3(S)

    if( S.p <= 0 ) 
      fprintf('Pop failed because of Empty Stack.\n');
    end
    
    val = S.v(S.p - 2); 
    S.v(S.p - 2) = S.v(S.p -1);
    S.v(S.p - 1) = S.v(S.p);
    S.p = S.p - 1;
end
%}

% ================================= 
% Pop the third value from the top
% =================================
function [S,val] = Stack_pop_bottom(S)

    if( S.p <= 0 ) 
      fprintf('Pop failed because of Empty Stack.\n');
    end
    
    val = S.v(1);
    S.v(1:S.p-1) = S.v(2:S.p);
    S.p = S.p - 1;
end

% ================================= 
% Replace  the third value from the top
% =================================
function [S,val] = Stack_replace_bottom(S, s_val)

    if( S.p <= 0 ) 
      fprintf('Pop failed because of Empty Stack.\n');
    end
    
    val = S.v(1);
    S.v(1) = s_val;
end

% ================================= 
% Pop abd return a data dt
% =================================
%{
function [S,val] = Stack_pop(S)

    if( S.p <= 0 ) 
      fprintf('Pop failed because of Empty Stack.\n');
    end
    
    val = S.v(S.p);
    S.pos = S.pos - 1;
end
%}