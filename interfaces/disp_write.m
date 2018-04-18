function disp_write (D,filename)
% saves disparity map D to png file

D = double(D);

I = D*256;
I(D==0) = 1;
I(I<0) = 0;
I(I>65535) = 0;
I = uint16(I);
imwrite(I,filename);

