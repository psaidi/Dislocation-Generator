function atomsInNew_Del_Shifted=StretchedShiftbackNoComp(atomsInNew_Del,Wmin,hmin,b,SPH)
%relaxing the system at the distance bigger than h min
%Total displacement as a function of Z is determined
%this part works for Z>=0
%Slip plane height=SPH
atomsInNew_Del(:,4)=abs(b./2.*(1-abs(atomsInNew_Del(:,3))./hmin));


M = atomsInNew_Del(:,4);
M(abs(atomsInNew_Del(:,3))>=hmin)=0;       % displacement is zero for atoms further than hmin in z direction
atomsInNew_Del(:,4)=M(:,1);

atomsInNew_Del(:,5)=atomsInNew_Del(:,4)./(atomsInNew_Del(:,4)-Wmin).*(abs(atomsInNew_Del(:,1))-Wmin);
MM=atomsInNew_Del(:,5);
MM(abs(atomsInNew_Del(:,3))>=hmin)=0;
MM(abs(atomsInNew_Del(:,1))>=Wmin)=0;
MM(atomsInNew_Del(:,3)<SPH)=0;
atomsInNew_Del(:,5)=MM(:,1);

% if x<=0 then we add the abs(dx) and if it is positinve we subtract abs(dx)
atomsInNew_Del(:,1)=atomsInNew_Del(:,1)-sign(atomsInNew_Del(:,1)).*abs(atomsInNew_Del(:,5));
%atomsInNew_Del_Shifted=atomsInNew_Del(:,1:3);



% %We displace the atoms in the lower part to compensate for the half of
% %distortion 
% %dl
% atomsInNew_Del(:,6)=abs(b./4.*(1-abs(atomsInNew_Del(:,3))./hmin));          %This is displcement as a function of Z
% 
% 
% N = atomsInNew_Del(:,6);
% N(atomsInNew_Del(:,3)<=-hmin | atomsInNew_Del(:,3)>SPH)=0;       % displacement is zero for atoms further than hmin in z direction or at the upper half
% atomsInNew_Del(:,6)=N(:,1);
% [R,C]=size(atomsInNew_Del);
% %dx
% %abs(b./4.*(1-abs(atomsInNew_Del(:,3))./hmin))=dl
% %dx=dl/Wmin*x  if x<wmin
% NN=zeros(R,1);
% NN(abs(atomsInNew_Del(:,1))<Wmin)=atomsInNew_Del(abs(atomsInNew_Del(:,1))<Wmin,6).*(abs(atomsInNew_Del(abs(atomsInNew_Del(:,1))<Wmin,1)))./Wmin; 
% NN(abs(atomsInNew_Del(:,1))>=Wmin)=atomsInNew_Del(abs(atomsInNew_Del(:,1))>=Wmin,6);
% NN(atomsInNew_Del(:,3)>=SPH | abs(atomsInNew_Del(:,3))>=hmin)=0;
% atomsInNew_Del(:,7)=NN(:,1);
% 
% % if x<=0 then we add the abs(dx) and if it is positinve we subtract abs(dx)
% atomsInNew_Del(:,1)=atomsInNew_Del(:,1)-sign(atomsInNew_Del(:,1)).*abs(atomsInNew_Del(:,7));
% atomsInNew_Del_Shifted=atomsInNew_Del(:,1:3);
% 

%Shifting back the atoms
atomsInNew_Del_Shifted=atomsInNew_Del(:,1:3);
dx=b/2;
atomsInNew_Del_Shifted(:,1)=atomsInNew_Del_Shifted(:,1)+dx;












end