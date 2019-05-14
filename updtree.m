%%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ! This subroutine updates a value on the binary tree. The FMM
% ! should only produce updated values that are less than their
% ! prior values.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function updtree(iz,ix)

global nsts btg ttn
% 
% ! ix,iz = grid position of new addition to tree
% ! tpp = tree position of parent
% ! tpc = tree position of child
% ! exch = dummy to exchange btg values
% !
% ! Filter the updated value to its correct position
% !
tpc=nsts(iz,ix);
tpp=round(tpc/2);       %%%%%%%% kt
while tpp>0
    if ttn(iz,ix)<ttn(btg(tpp,1),btg(tpp,2))
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

end


