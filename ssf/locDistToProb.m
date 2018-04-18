function probDist = locDistToProb(locDist,var)
% converts distances between superpixels to a discrete distribution 
% so that samples can be drawn using discretesample

  nS = size(locDist,1);  
  prob = exp(-(locDist./var).^2);
  prob(logical(eye(nS))) = 0;
  probDist = prob./(sum(prob,2)*ones(1,nS));

end