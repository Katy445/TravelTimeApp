function [biv]=bilinear(nv,dsx,dsz,dnx,dnz)
% 
% 
% ! nv = four node vertex values
% ! dsx,dsz = distance between internal point and top left node
% ! dnx,dnz = width and height of node rectangle
% ! biv = value at internal point calculated by bilinear interpolation
% ! produ = product variable

biv=0;
for i=1:2
   for j=1:2
      produ=(1-abs(((i-1)*dnx-dsx)/dnx))*(1-abs(((j-1)*dnz-dsz)/dnz));
      biv=biv+nv(i,j)*produ;
   end
end

end