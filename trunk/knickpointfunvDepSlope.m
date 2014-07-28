%%This function uses a constant rate (R) of knickpoint erosion to move the
%%headcell and headedge upstream.
%Edited May 23, 2012. Took out Si

function[headedge,headcell, xcell,zcell, xedge,cumero]=Knickpoint(xcellin,zcellin, t,dT, R,knickthresh,xedge,headedge,headcell, xstar,Stop,dx,xedgeOri,cumero, Sbtm, hcthresh,stillretreat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Constants
xcell0=xcellin;%m There are more cells than edges
xcell=xcell0;% use if not in if loop.
zcell0=zcellin;%m

dx0=diff(xcell0);
%Now set up variables needed to go into the loop
xedge0=xedge;%I need the original xcell, called xcell0, to be preserved. So I initialize this variable and make changes to it, preserving xcell0.
zcell=zcell0;

zheadO=zcell(1,headcell);
xedgeO=xedge(1,headedge);
xcellO=xcell(1,headcell);
smth=xcell(1,headcell)-xcell(1,headcell-1);
xedge(1,headedge)=xedge(1,headedge)-(R*dT);
xcell(1,headcell)=mean([xedge(1,headedge) xedge(1,headedge-1)]);
xcell(1,headcell+1)=mean([xedge(1,headedge) xedge(1,headedge+1)]);
dx4=xedgeO-xedge(1,headedge);%R*dt; %xcell(1,headcell)-xcell(1,headcell-1);%dx for slope at headedge-1
%The calc. below should be right. However, if I want to lock in the slope
%instead of just making the slope equal to the local upstream slope, I
%could use a constant slope in the equation instead of dz/dx which I am
%currently using and calculating manually.
zcell(1,headcell)=zheadO-((zheadO-zcell(1,headcell-1))/(smth))*(xcellO-xcell(1,headcell));
%Now deposit the sediment downstream
erodevol=R*dT*(zcell(1,headcell)-zcell(1,headcell+1));
lendwnstr=length(xcell(1,headcell+1:end-1));   
dumx=(1:lendwnstr);%Can't just use x value b/c that would lead to really small xprop values
xprop=(1/xstar)*exp(-dumx./xstar);%proportion of sed going into each cell determined by the exponential distribution
sedpercell=(erodevol).*xprop;%erodevol.*xprop;%
zcell(1,headcell+1:end-1)=zcell(1,headcell+1:end-1)+sedpercell;


 if xedge(1,headedge) <= (xedge(1,headedge-1)+dx0(1)/2);%xedgeOri(1,headedge)-xedge(1,headedge)== 0.5; %If you erode the half of the headcell
     xedge(1,headedge-1)=xedge(1,headedge);%Move the cell center to the left over to the right
     xedge(1,headedge)=xedgeOri(1,headedge);% and make the current xcell center equal to where it started
     xcell(1,headcell)=mean([xedge(1,headedge) xedge(1,headedge-1)]);
     xcell(1,headcell-1)=mean([xedge(1,headedge-1) xedge(1,headedge-2)]);
     dx3=(xcell(1,headcell-1))-(xcell(1,headcell-2));%dx for slope at headedge-1 Remember haven't changed the headcell yet
     %what happens if I change from end-2 to end-1
     %Not much. I'll leave it like that for now
     zcell(1,headcell-1)=zcell(1,headcell-2)-Stop(1,end-2)*dx3;%Stop(1,end-1)*dx3;
     erodevol=R*dT*(zcell(1,headcell)-zcell(1,headcell+1));
     zcell(1,headcell)=zcell(1,headcell+1);
     lendwnstr=length(xcell(1,headcell:end-1));
     dumx=(1:lendwnstr);%Can't just use x value b/c that would lead to really small xprop values
     xprop=(1/xstar)*exp(-dumx./xstar);%proportion of sed going into each cell determined by the exponential distribution
     sedpercell=erodevol.*xprop;%(erodevol./dx(1,headcell:end-1)).*xprop;
     zcell(1,headcell:end-1)=zcell(1,headcell:end-1)+sedpercell;
     headedge=headedge-1;%and change the headcell location
     headcell=headedge-1;
     
 end
 end
     
     
 
 