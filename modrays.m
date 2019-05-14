
function [travel_time_field]=modrays(scx,scy,map)
format long
global nsts btg veln dnx dnz goz gox nnx nnz ttn

veln=map;
nnz=size(veln,1);
nnx=size(veln,2);
dnx=1;
dnz=1;
gox=0;
goz=0;

% !
% ! Allocate memory for node status and binary trees
% !

snb=0.5;
nsts=zeros(nnz,nnx);
maxbt=round(snb*nnx*nnz);
btg=zeros(maxbt,2);
% !
% ! Loop through all sources and find traveltime fields
% !
x=scx;
z=scy;
% !     Call a subroutine that works out the first-arrival traveltime
% !     field.
travel(x,z)
travel_time_field=ttn;
end