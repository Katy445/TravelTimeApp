%%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ! This subroutine updates the binary tree after the root
% ! value has been used. The root is replaced by the value
% ! at the bottom of the tree, which is then filtered down
% ! to its correct position. This ensures that the tree remains
% ! balanced.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function downtree
global nsts btg ntr ttn
% !
% ! tpp = tree position of parent
% ! tpc = tree position of child
% ! exch = dummy to exchange btg values
% ! rd1,rd2 = substitution variables
% !
% ! Replace root of tree with its last value

if ntr==1
    ntr=ntr-1;
    return
end
nsts(btg(ntr,1),btg(ntr,2))=1;
btg(1,:)=btg(ntr,:);
% !
% ! Reduce size of tree by one
% !
ntr=ntr-1;
% !
% ! Now filter new root down to its correct position
% !
tpp=1;
tpc=2*tpp;
while tpc<ntr

% ! Check which of the two children is smallest - use the smallest
% !

   rd1=ttn(btg(tpc,1),btg(tpc,2));
   rd2=ttn(btg(tpc+1,1),btg(tpc+1,2));
   if rd1>rd2
       tpc=tpc+1;
   end
% !
% !  Check whether the child is smaller than the parent; if so, then swap,
% !  if not, then we are done
% !
   rd1=ttn(btg(tpc,1),btg(tpc,2));
   rd2=ttn(btg(tpp,1),btg(tpp,2));
   if rd1<rd2
      nsts(btg(tpp,1),btg(tpp,2))=tpc;
      nsts(btg(tpc,1),btg(tpc,2))=tpp;
      exch=btg(tpc,:);
      btg(tpc,:)=btg(tpp,:);
      btg(tpp,:)=exch;
      tpp=tpc;
      tpc=2*tpp;
   else
      tpc=ntr+1;
   end
end

% !
% ! If ntr is an even number, then we still have one more test to do
% !

if tpc==ntr
   rd1=ttn(btg(tpc,1),btg(tpc,2));
   rd2=ttn(btg(tpp,1),btg(tpp,2));
   if rd1<rd2
      nsts(btg(tpp,1),btg(tpp,2))=tpc;
      nsts(btg(tpc,1),btg(tpc,2))=tpp;
      exch=btg(tpc,:);
      btg(tpc,:)=btg(tpp,:);
      btg(tpp,:)=exch;
   end
end
end
