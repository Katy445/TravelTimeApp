function travel(scx,scz)

global nsts btg ntr ttn veln gox goz dnx dnz nnx nnz

% ! isx,isz = grid cell indices (i,j,k) which contains source
% ! scx,scz = (r,x,y) location of source
% ! sw = a switch (0=off,1=on)
% ! ix,iz = j,k position of "close" point with minimum traveltime
% ! maxbt = maximum size of narrow band binary tree
% ! rd2,rd3 = substitution variables
% ! vsrc = velocity at source
% ! vss = velocity at nodes surrounding source
% ! dsx, dsz = distance from source to cell boundary in x and z
% ! ds = distance from source to nearby node
% ! urg = use refined grid (0=no,1=yes,2=previously used)
% ! swrg = switch to end refined source grid computation
% !
% ! The first step is to find out where the source resides
% ! in the grid of nodes. The cell in which it resides is
% ! identified by the "north-west" node of the cell. If the
% ! source lies on the edge or corner (a node) of the cell, then
% ! this scheme still applies.


isx=round((scx-gox)/dnx)+1;
isz=round((scz-goz)/dnz)+1;

if isx==nnx
    isx=isx-1;
end
if isz==nnz
    isz=isz-1;
end

% Set all values of nsts to -1 if beginning from a source point.

nsts=nsts-1;

% set initial size of binary tree to zero

ntr=0;

% !  In general, the source point need not lie on a grid point.
% !  Bi-linear interpolation is used to find velocity at the
% !  source point.

vss=zeros(2,2);

for i=1:2
    for j=1:2
        vss(i,j)=veln(isz-1+j,isx-1+i);

    end
end

dsx=(scx-gox)-(isx-1)*dnx;
dsz=(scz-goz)-(isz-1)*dnz;
%    
vsrc=bilinear(vss,dsx,dsz,dnx,dnz);


% !  Now find the traveltime at the four surrounding grid points. This
% !  is calculated approximately by assuming the traveltime from the
% !  source point to each node is equal to the the distance between
% !  the two points divided by the average velocity of the points

for i=1:2
    for j=1:2
        ds=sqrt((dsx-(i-1)*dnx)^2+(dsz-(j-1)*dnz)^2);
         vss_temp=vss(i,j);
         vsrc_temp=round(vsrc);
        ttn(isz-1+j,isx-1+i)=2.0*ds/(vss_temp+vsrc_temp);
        addtree(isz-1+j,isx-1+i)
    end
end

% ! Now calculate the first-arrival traveltimes at the
% ! remaining grid points. This is done via a loop which
% ! repeats the procedure of finding the first-arrival
% ! of all "close" points, adding it to the set of "alive"
% ! points and updating the points surrounding the new "alive"
% ! point. The process ceases when the binary tree is empty,
% ! in which case all grid points are "alive".

while ntr>0

% ! Set the "close" point with minimum traveltime
% ! to "alive"
% !
   ix=btg(1,2);
   iz=btg(1,1);
   nsts(iz,ix)=0;

% ! Update the binary tree by removing the root and
% ! sweeping down the tree.
% !
   downtree

% ! Now update or find values of up to four grid points
% ! that surround the new "alive" point.
% !
% ! Test points that vary in x
% !
for i=ix-1:2:ix+1
    if i>=1 && i<=nnx
        if nsts(iz,i)==-1

% ! This option occurs when a far point is added to the list
% ! of "close" points
% !
         fouds14(iz,i)

         addtree(iz,i)
        elseif nsts(iz,i)>0
% !
% ! This happens when a "close" point is updated
% !
          fouds14(iz,i)

          updtree(iz,i)
        end
    end
end

% !
% ! Test points that vary in z
% !
for i=iz-1:2:iz+1
    if i>=1 && i<=nnz
        if nsts(i,ix)==-1

% ! This option occurs when a far point is added to the list
% ! of "close" points
% !

            fouds14(i,ix)

            addtree(i,ix)
        elseif nsts(i,ix)>0
% !
% ! This happens when a "close" point is updated
% !
            fouds14(i,ix)

            updtree(i,ix)
        end
    end
end
end
end
