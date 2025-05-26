function  DistanceProfile = Spikelet_STOMP_DistanceProfile(D, subseq, Param)
 FuncName = 'Spikelet_STOMP_DistanceProfile';
%%Anonymous-Author information blinded for review
%% This is the revised version of source code for the ICDM 2017 paper 
%% "Time Series Chains: A New Primitive for Time Series Data Mining"

%%input:
%%A: Time Series
%%SubsequenceLength: Subsequence Length

SubseqNorm = Param.subseq_norm;

%% Argument
if size(D,1) < size(D,2)
   D = D'; 
end
SubsequenceLength = length(subseq);
if SubsequenceLength > length(D)/2
    error('Error: Time series is too short relative to desired subsequence length');
end
if SubsequenceLength < 4
    error('Error: Subsequence length must be at least 4');
end

%% main  
[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(D, SubsequenceLength);
[DistanceProfile lastz dropval lastsumy lastsumy2]= fastfindNN(X, subseq, n, SubsequenceLength, sumx2, sumx, meanx, sigmax2, sigmax, SubseqNorm);

end

%figure;
%cur=ChainStart;
    
%for i=1:count-1
%    next=MPindexRight(cur);
%    subplot(ceil(count/3),3,i);
%    plot(cur:(cur+SubsequenceLength-1),-zscore(A(cur:(cur+SubsequenceLength-1)),1)+zscore(A(next:(next+SubsequenceLength-1)),1));
%    xlim([cur,cur+SubsequenceLength-1]);
%    ylim([curmindiff,curmaxdiff]);
%    if i==1
%        title('evolution trend');
%    end
%    cur=next;
%end

% m is winSize
function [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
     
  n = length(x);
  
  x(n+1:2*n) = 0;
  X = fft(x);
  
  cum_sumx = cumsum(x);
  cum_sumx2 =  cumsum(x.^2);
  sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
  sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
  
  meanx = sumx./m;
  
  sigmax2 = (sumx2./m)-(meanx.^2);
  sigmax = sqrt(sigmax2);
end

% m is winSieze
function [dist lastz dropval sumy sumy2] = fastfindNN(X, y, n, m, sumx2, sumx, meanx, sigmax2, sigmax, DistType)
%x is the data, y is the query
%y = (y-mean(y))./std(y,1);                      %Normalize the query
dropval=y(1);
y = y(end:-1:1);                                %Reverse the query
y(m+1:2*n) = 0;

%The main trick of getting dot products in O(n log n) time
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

%compute y stats -- O(n)
sumy = sum(y);
sumy2 = sum(y.^2);
meany=sumy/m;
sigmay2 = sumy2/m-meany^2;
sigmay = sqrt(sigmay2);

 if DistType == "Znorm"
    dist2 = 2*(m-(z(m:n)-m*meanx*meany)./(sigmax*sigmay));  % z-normalilzed Euclidean distance     
 else
    dist2 = sumx2 + sumy2 - 2*z(m:n);  
 end
dist = sqrt(dist2);
lastz=real(z(m:n));

end