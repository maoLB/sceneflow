1
change

data = dlmread(filename,' ',0,1);

to

datatemp  = data = dlmread(filename,' ',25,1);
data=[datatemp(1,:);datatemp(9,:)];

in loadOSFInput.m
to allow the datasheet of KITTI

2
in flow_write.m
change 

I(:,:,1) = uint16(max(min(shiftdim(F(:,:,1))*64+2^15,2^16-1),0));
I(:,:,2) = uint16(max(min(shiftdim(F(:,:,2))*64+2^15,2^16-1),0));
I(:,:,3) = uint16(max(min(shiftdim(F(:,:,3)),1),0));

to

I(:,:,1) = uint16(max(min(shiftdim(F(:,:,1))*512+2^15,2^16-1),0));
I(:,:,2) = uint16(max(min(shiftdim(F(:,:,2))*512+2^15,2^16-1),0));
I(:,:,3) = uint16(max(min(shiftdim(F(:,:,3)),1),0));
