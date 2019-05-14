%%
function fouds14(iz,ix)
global nsts ttn dnx dnz nnx nnz veln 

% ! ix = NS position of node coordinate for determination
% ! iz = EW vertical position of node coordinate for determination
% ! trav = traveltime calculated for trial node
% ! travm = minimum traveltime calculated for trial node
% ! slown = slowness at (iz,ix)
% ! tsw1 = traveltime switch (0=first time,1=previously)
% ! a,b,c,u,v,em = Convenience variables for solving quadratic
% ! tdsh = local traveltime from neighbouring node
% ! tref = reference traveltime at neighbouring node
% ! ri = Radial distance
% ! risti = ri*sin(theta) at point (iz,ix)
% ! rd1 = dummy variable
% ! swsol = switch for solution (0=no solution, 1=solution)
% !
% ! Inspect each of the four quadrants for the minimum time
% ! solution.

% CALCULATE 2ND ORDER TRAVEL TIMES ON 0DEG STENCIL

tsw1=0;
travm=0;
slown=1.0/veln(iz,ix);

for j=ix-1:2:ix+1
    if j>=1 && j<=nnx
        swj=-1;
        if j==ix-1
            j2=j-1;
            if j2>=1
                if nsts(iz,j2)==0
                    swj=0;
                end
            end
        else
            j2=j+1;
            if j2<=nnx
               if nsts(iz,j2)==0
                        swj=0;
               end
            end
        end
        if nsts(iz,j)==0 && swj==0
            swj=-1;
            if ttn(iz,j)>=ttn(iz,j2)
            swj=0;
            end
        else
            swj=-1;
        end
   for k=iz-1:2:iz+1
       if k>=1 && k<=nnz
           swk=-1;
           if k==iz-1
               k2=k-1;
               if k2>=1
                   if nsts(k2,ix)==0
                       swk=0;
                   end
               end
           else
               k2=k+1;
               if k2<=nnz
                   if nsts(k2,ix)==0
                       swk=0;
                   end
               end
           end
           if nsts(k,ix)==0 && swk==0
               swk=-1;
               if ttn(k,ix)>=ttn(k2,ix)
                   swk=0;
               end
           else
               swk=-1;
           end
            
           swsol=0;
           
           if swj==0
               swsol=1;
               if swk==0
                  u=2.0*dnx;
                  a=18*3.0;
                  b=-18*(4.0*ttn(iz,j)-ttn(iz,j2)+4.0*ttn(k,ix)-ttn(k2,ix));
                  c=9*(4.0*ttn(iz,j)-ttn(iz,j2))^2+9*(4.0*ttn(k,ix)-ttn(k2,ix))^2-3*4*u^2*slown^2;
                  tref=4.0*ttn(iz,j)-ttn(iz,j2);
                  tdiv=3.0;           
               elseif nsts(k,ix)==0
                  u=dnz;
                  v=2.0*dnx;
                  em=3.0*ttn(k,ix)+4.0*ttn(iz,j)+ttn(iz,j2);
                  a=18;
                  b=-6.0*em;
                  c=(3.0*ttn(k,ix))^2+(4.0*ttn(iz,j)+ttn(iz,j2))^2-3*4*v^2*slown^2;
                  tref=ttn(k,ix);
                  tdiv=1.0;
               else
                  u=2.0*dnx;
                  a=1.0;
                  b=0.0;
                  c=-u^2*slown^2;
                  tref=4.0*ttn(iz,j)-ttn(iz,j2);
                  tdiv=3.0;
               end
           elseif nsts(iz,j)==0
               swsol=1;
               if swk==0              
                  u=dnx;
                  v=2.0*dnz;
                  em=3.0*ttn(iz,j)+4.0*ttn(k,ix)+ttn(k2,ix);
                  a=18;
                  b=-6.0*em;
                  c=(3.0*ttn(iz,j))^2+(4.0*ttn(k,ix)+ttn(k2,ix))^2-3*4*u^2*slown^2;
                  tref=ttn(iz,j);
                  tdiv=1.0;       
               elseif nsts(k,ix)==0
                  u=dnx;
                  v=dnz;
                  a=2;
                  b=-2*(ttn(k,ix)+ttn(iz,j));
                  c=ttn(k,ix)^2+ttn(iz,j)^2-(u*slown)^2;
                  tref=ttn(iz,j);
                  tdiv=1.0; %%               
               else
                  a=1.0;
                  b=0.0;
                  c=-slown^2*dnx^2;
                  tref=ttn(iz,j);
                  tdiv=1.0;
               end
           else
               if swk==0
                  swsol=1;
                  u=2.0*dnz;
                  a=1.0;
                  b=0.0;
                  c=-u^2*slown^2;
                  tref=4.0*ttn(k,ix)-ttn(k2,ix);
                  tdiv=3.0;   
               elseif nsts(k,ix)==0
                  swsol=1;
                  a=1.0;
                  b=0.0;
                  c=-slown^2*dnz^2;
                  tref=ttn(k,ix);
                  tdiv=1.0;
               end
           end
          %Now find the solution of the quadratic equation  
            if swsol==1
               rd1=b^2-4.0*a*c;
               if rd1<0
                   rd1=0;
               end
               tdsh=(-b+sqrt(rd1))/(2.0*a);
               trav=(tref+tdsh)/tdiv;
               if tsw1==1
                   travm=min(trav,travm);
               else
                   travm=trav;
                   tsw1=1;
               end
            end
       end
   end
    end
end

% if travm~=0 
%     ttn(iz,ix)=min(travm,trav);
% else
%     ttn(iz,ix)=trav;
% end

tsw2=0;
travmd=0;
slown=1.0/veln(iz,ix);

mf2=sqrt(2);  % distance between two diagonally connected points  

for j=ix-1:2:ix+1
    for k=iz-1:2:iz+1
    if j>=1 && j<=nnx && k>=1 && k<=nnz
        swjj=-1;
        if j==ix-1 && k==iz-1
            j2=j-1;
            k2=k-1;
            if j2>=1 && k2>=1
                if nsts(k2,j2)==0
                    swjj=0;
                end
            end
        elseif j==ix-1 && k==iz+1
            j2=j-1;
            k2=k+1;
            if j2>=1 && k2<=nnz
               if nsts(k2,j2)==0
                  swjj=0;
               end
            end
        elseif j==ix+1 && k==iz-1
            j2=j+1;
            k2=k-1;
            if j2<=nnx && k2>=1
                if nsts(k2,j2)==0
                    swjj=0;
                end
            end
        else
            j2=j+1;
            k2=k+1;
            if j2<=nnx && k2<=nnz
                if nsts(k2,j2)==0
                    swjj=0;
                end
            end           
        end
        if nsts(k,j)==0 && swjj==0
            swjj=-1;
            if ttn(k,j)>=ttn(k2,j2)
            swjj=0;
            end
        else
            swjj=-1;
        end  
        
      
	 for jj=ix-1:2:ix+1		
	   for kk=iz-1:2:iz+1
        if jj>=1 && jj<=nnx && kk>=1 && kk<=nnz
           swkk=-1;
           if kk==iz-1 && jj==ix-1
               kk2=kk-1;
               jj2=jj-1;
               if kk2>=1 && jj2>=1
                   if nsts(kk2,jj2)==0
                       swkk=0;
                   end
               end
           elseif kk==iz-1 && jj==ix+1
               kk2=kk-1;
               jj2=jj+1;
               if kk2>=1 && jj2<=nnz
                   if nsts(kk2,jj2)==0
                       swkk=0;
                   end
               end
           elseif kk==iz+1 && jj==ix-1
               kk2=kk+1;
               jj2=jj-1;
               if kk2<=nnz && jj2>=1
                   if nsts(kk2,jj2)==0
                       swkk=0;
                   end
               end
           else 
               kk2=kk+1;
               jj2=jj+1;
               if kk2<=nnz && jj2<=nnx
                   if nsts(kk2,jj2)==0
                       swkk=0;
                   end
               end
           end

           if nsts(kk,jj)==0 && swkk==0
               swkk=-1;
               if ttn(kk,jj)>=ttn(kk2,jj2)
                   swkk=0;
               end
           else
               swkk=-1;
           end
           
           swsol=0;
           
           if swjj==0 
               swsol=1;
               if swkk==0
                  u=mf2*2.0*dnx;
                  a=18*3.0;
                  b=-18*(4.0*ttn(k,j)-ttn(k2,j2)+4.0*ttn(kk,jj)-ttn(kk2,jj2));
                  c=9*(4.0*ttn(k,j)-ttn(k2,j2))^2+9*(4.0*ttn(kk,jj)-ttn(kk2,jj2))^2-3*4*u^2*slown^2;
                  tref=4.0*ttn(k,j)-ttn(k2,j2);
                  tdiv=3.0;
               elseif nsts(kk,jj)==0
                  u=mf2*dnz;
                  v=mf2*2.0*dnx;
                  em=3.0*ttn(kk,jj)+4.0*ttn(k,j)+ttn(k2,j2);
                  a=18;
                  b=-6.0*em;
                  c=(3.0*ttn(kk,jj))^2+(4.0*ttn(k,j)+ttn(k2,j2))^2-3*4*v^2*slown^2;
                  tref=ttn(k,j);
                  tdiv=1.0;
               else
                  u=mf2*2.0*dnx;
                  a=1.0;
                  b=0.0;
                  c=-u^2*slown^2;
                  tref=4.0*ttn(k,j)-ttn(k2,j2);
                  tdiv=3.0; %%
               end
           elseif nsts(k,j)==0
               swsol=1;
               if swkk==0
                  u=mf2*dnx;
                  v=mf2*2.0*dnz;
                  em=3.0*ttn(k,j)+4.0*ttn(kk,jj)+ttn(kk2,jj2);
                  a=18;
                  b=-6.0*em;
                  c=(3.0*ttn(k,j))^2+(4.0*ttn(kk,jj)+ttn(kk2,jj2))^2-3*4*u^2*slown^2;
                  tref=ttn(k,j);
                  tdiv=1.0;
               elseif nsts(kk,jj)==0
                  u=mf2*dnx;
                  v=mf2*dnz;
                  a=2;
                  b=-2*(ttn(k,ix)+ttn(iz,j));
                  c=ttn(k,ix)^2+ttn(iz,j)^2-(u*slown)^2;
                  tref=ttn(iz,j);
                  tdiv=1.0; %%
               else
                  u=mf2*dnx; 
                  a=1.0;
                  b=0.0;
                  c=-slown^2*u^2;
                  tref=ttn(k,j);
                  tdiv=1.0; %%
               end
           else
               if swkk==0
                  swsol=1;
                  u=2.0*mf2*dnz;
                  a=1.0;
                  b=0.0;
                  c=-u^2*slown^2;
                  tref=4.0*ttn(kk,jj)-ttn(kk2,jj2);
                  tdiv=3.0; %%
               elseif nsts(kk,jj)==0
                  swsol=1;
                  u=mf2*dnx;
                  a=1.0;
                  b=0.0;
                  c=-slown^2*u^2;
                  tref=ttn(kk,jj);
                  tdiv=1.0; %%
               end
           end
          %Now find the solution of the quadratic equation  
            if swsol==1
               rd1=b^2-4.0*a*c;
               if rd1<0
                   rd1=0;
               end
               tdsh=(-b+sqrt(rd1))/(2.0*a);
               trav=(tref+tdsh)/tdiv;
               if tsw2==1
                   travmd=min(trav,travmd);
               else
                   travmd=trav;
                   tsw2=1;
               end
            end
       end
   end
    end
    end
    end
end


if travmd~=0 
    travmd=min(travm,travmd);
else
    travmd=travm;
end

ttn(iz,ix)=travmd;
end
