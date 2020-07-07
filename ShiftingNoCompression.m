function atomsInNew_Del_Shifted=ShiftingNoCompression(atomsInNew_Del,hmin,b,safty)
%Shifting all atoms so the centre of burgers vector is at x=0  
dx=b/2;
atomsInNew_Del(:,1)=atomsInNew_Del(:,1)-dx;

%Shifting the atoms to make up for the parts out of the domain of deformation burgers vector.
atomsInNew_Del(:,3+1)=b./2.*atomsInNew_Del(:,3)./hmin;
Nx=atomsInNew_Del(:,1);
Nz=atomsInNew_Del(:,3);
dx2=b/2;
N=atomsInNew_Del(:,4);
%if Z>0 atoms will shift otherwise no. 
N(Nz>=hmin)=dx2;
N(Nz<safty)=0;
% if x<=0 then we add the abs(dx) and if it is positinve we subtract abs(dx)
atomsInNew_Del(:,1)=atomsInNew_Del(:,1)-sign(atomsInNew_Del(:,1)).*abs(N(:,1));
atomsInNew_Del_Shifted=atomsInNew_Del(:,1:3);

end