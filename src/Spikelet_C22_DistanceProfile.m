function  DistanceProfile = Spikelet_C22_DistanceProfile(D, subseq, Param,ExclusionZone)
 FuncName = 'Spikelet_C22_DistanceProfile';

 SubseqNorm = Param.subseq_norm;
 %% Argument
 if size(D,1) < size(D,2)
   D = D'; 
 end
 SubsequenceLength = length(subseq);
 if SubsequenceLength > length(D)/2
    error('Error: Time series is too short relative to desired subsequence length');
 end
 if SubsequenceLength < 3
    error('Error: Subsequence length must be at least 3');
 end

 %% main  
 DistanceProfile_length = length(D)-length(subseq)+1;
 DistanceProfile = NaN(DistanceProfile_length,1);
 for i = 1:length(D)-length(subseq)+1
  subseq_i = D(i:i+length(subseq)-1);     
  if ExclusionZone(i) == true
   DistanceProfile(i) = Inf;   
  else
   DistanceProfile(i) = c22_norm(subseq,subseq_i,SubseqNorm);
  end
 end

end

function dist_ij = c22_norm(subseq_i,subseq_j,SubseqNrom)
  if SubseqNrom == "Znorm"
    subseq_i = normalize(subseq_i);
    subseq_j = normalize(subseq_j);
  end

  c22_i = C22_call(subseq_i);
  c22_j = C22_call(subseq_j);
  %　inner product
  dot_product = dot(c22_i, c22_j);

  % inner produ
  dist_ij = 1 - dot_product;
  % other options: arccos(dot_product), sqrt(1-dot_product^2), (1-dot_product)/2,Mahalanobis distance)

end

function dist_ij = euclidean_norm(subseq_i,subseq_j,DistType)
  if DistType == "Znorm"
    subseq_i = normalize(subseq_i);
    subseq_j = normalize(subseq_j);
  end
    dist_ij = norm(subseq_i-subseq_j);
end