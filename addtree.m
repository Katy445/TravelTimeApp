%%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ! This subroutine adds a value to the binary tree by
% ! placing a value at the bottom and pushing it up
% ! to its correct position.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function addtree(iz,ix)

global nsts btg ntr ttn
% !
% ! ix,iz = grid position of new addition to tree
% ! tpp = tree position of parent
% ! tpc = tree position of child
% ! exch = dummy to exchange btg values
% !
% ! First, increase the size of the tree by one.
% !
ntr=ntr+1;
% !
% ! Put new value at base of tree
% !
nsts(iz,ix)=ntr;
btg(ntr,2)=ix;
btg(ntr,1)=iz;
% !
% ! Now filter the new value up to its correct position
% !
tpc=ntr;
tpp=round(tpc/2);   %%%%%%%%%% KT
while tpp>0
    aa=btg(tpp,1);
    bb=btg(tpp,2);
   if ttn(iz,ix)<ttn(aa,bb)
      nsts(iz,ix)=tpp;
      nsts(btg(tpp,1),btg(tpp,2))=tpc;
      exch=btg(tpc,:);
      btg(tpc,:)=btg(tpp,:);
      btg(tpp,:)=exch;
      tpc=tpp;
      tpp=round(tpc/2);
   else
      tpp=0;
   end
end
% END SUBROUTINE addtree
end
