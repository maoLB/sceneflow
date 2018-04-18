function [objectParticles,objectParticles_struct] = osf_sampleObjectParticles(object,stdObject,objectParticles,stdBackground)
%SAMPLESHAPEPARTICLES generates particles starting from current motion 
% parameters

    nO = size(object,1);                % number of objects
    [~,~,nP] = size(objectParticles);   % number of particles

    objectParticles_struct = {};
    
    %% The first object is treated seperately since it is the background motion
    if ~exist('stdBackground','var') || isempty(stdBackground)
      stdBackground = stdObject;
    end
    
    iO = 1;
    
    % first particle is current MAP solution
    objectParticles(iO,:,1) = object(iO,:,1);

    % sample particle 2:nP (=end) of object i0 (for each dimension)
    for d = 1:6
        objectParticles(iO,d,2:end) = normrnd(object(iO,d,1), stdBackground(d), nP-1, 1);
    end
    
    %% Remaining objects
    if nO > 1
      for iO = 2:nO % all objects

          % first particle is current MAP solution
          objectParticles(iO,:,1) = object(iO,:,1);

          % sample particle 2:nP (=end) of object i0 (for each dimension)
          for d = 1:6
              objectParticles(iO,d,2:end) = normrnd(object(iO,d,1), stdObject(d), nP-1, 1);
          end
      end   
    end
    
end