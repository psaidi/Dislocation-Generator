clear
% Terms that we determine 
%P1: on the slip plane beginiing of burgers vector
%P2: on the slip plane, end of burgers vector
%P3: A thirs point on the slip plane
%N_ucX= Size of the initital system in x direction
%N_ucY= Size of the initital system in y direction
%N_ucZ= Size of the initital system in z direction
%safty Safty for making sure that the layer of atoms at the slip plane remains onchanged
%hmin the height that beyond that the distortion of the system is zero;
%Wmin the width that beyond that the distortion of the system is zero;
%Sizex is the size of the system after trimming in the direction of the
%burgers vector
%Sizey is the size of the system after trimming in the direction of the
%dislocation line 
%SizeZ is the size of the system after trimming in the direction of the
%normal to the slip plane vector

% Dumpfile 
savedir0= pwd; savedir=[savedir0 '/'];
%position=importdata([savedir FileName]);

%Atom_Pos=position.data;

LP=3.234;
a=LP*[1,sqrt(3),1.598];
NewDir=[a(1),0,0
        0,a(2),0
        0,0,a(3)];
    
coeff=1;
%Three atoms on the slip system. P2-P1 is the burgers vector and P3-P1 is
%the other vector in the slip system uses for finding the normal
P1=[0,0,0];
P2=[a(1),0,-a(3)];
P3=[a(1),a(2),0];

b_x=P2-P1;     %equivalent to X direction (burgers vector)
L_b_x=norm(b_x);
A=P3-P1;
N_z=cross(b_x,A);     %equivalent to Z direction (Normal to slip plane)
L_N_z=norm(N_z);
Y_y=cross(N_z,b_x);     %equivalent to Y direction (third vector) [Dislocation line]
L_Y_y=norm(Y_y);
a_new=[L_b_x,L_Y_y,L_N_z];


%Mapping the original orientation to the Global orientation
OrgDir0=[
        b_x/norm(b_x);
        Y_y/norm(Y_y);
        N_z/norm(N_z)];

%normalizing the original crystal orientation    
    OrgDir=zeros(3,3);
    for i=1:3
        OrgDir1(i,:)=OrgDir0(i,:)/norm(OrgDir0(i,:));
        NewDir1(i,:)=NewDir(i,:)/norm(NewDir(i,:));
    end
%finding the mapping matrix   
B=inv(OrgDir1)*NewDir1;

% Basis atoms for parents
 BasisAtoms=[0.000000,0.000000,0.000000;
             0.500000,0.166667,0.500000;
             0.500000,0.500000,0.000000;
             0.000000,0.666667,0.500000];

         
     
%Expanding the strucutre in the original coordination
[NBasis,c]=size(BasisAtoms);
N_ucX=100;
N_ucY=60;
N_ucZ=60;
LowX=-N_ucX; HighX=N_ucX;  
LowY=-N_ucY; HighY=N_ucY;  
LowZ=-N_ucZ; HighZ=N_ucZ;  


for i=1:NBasis
Xvec(i,:)=LowX:1:HighX;
Xvec(i,:)=Xvec(i,:)+BasisAtoms(i,1);
Yvec(i,:)=LowY:1:HighY;
Yvec(i,:)=Yvec(i,:)+BasisAtoms(i,2);
Zvec(i,:)=LowZ:1:HighZ;
Zvec(i,:)=Zvec(i,:)+BasisAtoms(i,3);
atomsInOriginal0(i,:,:)=combvec(Xvec(i,:),Yvec(i,:),Zvec(i,:))';
end

[Q,QQ,QQQ]=size(atomsInOriginal0);
p=0;
q=QQ;
for i=1:NBasis
  atomsInOriginal(p+1:q,1:3)= atomsInOriginal0(i,:,:);
  p=q;
  q=(i+1)*QQ;
end


atomsInOriginal(:,1)=atomsInOriginal(:,1).*a(1,1);
atomsInOriginal(:,2)=atomsInOriginal(:,2).*a(1,2);
atomsInOriginal(:,3)=atomsInOriginal(:,3).*a(1,3);

[Natoms,c]=size(atomsInOriginal);

%Original atoms

%FileName0='Original.dump';
%DumpMaker(savedir, atomsInOriginal,FileName0)


%Rotated atoms
atoms_Rot=atomsInOriginal*B;

%FileName0='Rot.dump';
%DumpMaker(savedir, atoms_Rot,FileName0)


%Triming the area around the box in the new direction
Sizex=L_b_x*38; %burgers vector
Sizey=L_b_x*38; %normal to slip plane
Sizez=L_b_x*21; %Dislocation line
atoms_Rot_trimmed1=atoms_Rot(atoms_Rot(:,1)>=-Sizex & atoms_Rot(:,2)>=-Sizey & atoms_Rot(:,3)>=-Sizez,:);
atoms_Rot_trimmed=atoms_Rot_trimmed1(atoms_Rot_trimmed1(:,1)<Sizex & atoms_Rot_trimmed1(:,2)<Sizey & atoms_Rot_trimmed1(:,3)<Sizez,:);


%FileName0='Rot__trimmed.dump';
%DumpMaker(savedir,atoms_Rot_trimmed,FileName0)



b_x_norm=b_x./norm(b_x);
N_z_norm=N_z./norm(N_z);
LMT_1=b_x_norm(1)*P1(1)+b_x_norm(2)*P1(2)+b_x_norm(3)*P1(3);
LMT_2=b_x_norm(1)*P2(1)+b_x_norm(2)*P2(2)+b_x_norm(3)*P2(3);
LMT_SLIP=N_z_norm(1)*P1(1)+N_z_norm(2)*P1(2)+N_z_norm(3)*P1(3);
X_0=atoms_Rot_trimmed(:,1);
Y_0=atoms_Rot_trimmed(:,2);
Z_0=atoms_Rot_trimmed(:,3);
atoms_Rot_trimmed_Del=atoms_Rot_trimmed;

safty=-0.01; %Safty for making sure that the layer of atoms at the slip plane remains onchanged
% safty is zero for forming stacking fault
%this command deleted the atoms in the half plane to form the dislocation 
% atoms_Rot_trimmed_Del(N_z_norm(1)*X_0+N_z_norm(2)*Y_0+N_z_norm(3)*Z_0>LMT_SLIP+safty & ...
%           b_x_norm(1)*X_0+b_x_norm(2)*Y_0+b_x_norm(3)*Z_0>LMT_1 & ...
%           b_x_norm(1)*X_0+b_x_norm(2)*Y_0+b_x_norm(3)*Z_0<=LMT_2,:)=[];
      
      
atoms_Rot_trimmed_Del(Z_0>LMT_SLIP+safty & ...
          X_0>LMT_1 & ...
          X_0<=LMT_2,:)=[];
    

%FileName0='Rot_trimmed_Del.dump';
%DumpMaker(savedir, atoms_Rot_trimmed_Del,FileName0)



hmin=40;
atoms_Rot_trimmed_Del_Shifted=ShiftingNoCompression(atoms_Rot_trimmed_Del,hmin,L_b_x,safty);

%FileName0='Rot_trimmed_Del_Shifted.dump';
%DumpMaker(savedir, atoms_Rot_trimmed_Del_Shifted,FileName0)

Wmin=10;
SPH=safty;
atoms_Rot_trimmed_Del_Shifted_Stretched=StretchedShiftbackNoComp(atoms_Rot_trimmed_Del_Shifted,Wmin,hmin,L_b_x,SPH);

%FileName0='Rot_trimmed_Del_Shifted_Stretched.dump';
%DumpMaker(savedir, atoms_Rot_trimmed_Del_Shifted_Stretched,FileName0)

% Produce the files for Lammps
%We use the results after trimming without rotation back
FileName1=['Disl_ReadableW' num2str(Wmin) '.dump'];
DumpreadableMaker(savedir,atoms_Rot_trimmed_Del_Shifted_Stretched,FileName1)

% %Triming the area around the box in the new direction
% Sizex=L_b_x*26; %burgers vector
% Sizey=L_b_x*20; %normal to slip plane
% Sizez=L_b_x*20; %Dislocation line
% atoms_Del_Rot_Shifted_Stretched_trimmed1=atoms_Del_Rot_Shifted_Stretched(atoms_Del_Rot_Shifted_Stretched(:,1)>=-Sizex & atoms_Del_Rot_Shifted_Stretched(:,2)>=-Sizey & atoms_Del_Rot_Shifted_Stretched(:,3)>=-Sizez,:);
% atoms_Del_Rot_Shifted_Stretched_trimmed=atoms_Del_Rot_Shifted_Stretched_trimmed1(atoms_Del_Rot_Shifted_Stretched_trimmed1(:,1)<Sizex & atoms_Del_Rot_Shifted_Stretched_trimmed1(:,2)<Sizey & atoms_Del_Rot_Shifted_Stretched_trimmed1(:,3)<Sizez,:);
% 
% 
% 
% FileName0='Del_Rot_Sifted_Stretched_trimmed.dump';
% DumpMaker(savedir,atoms_Del_Rot_Shifted_Stretched_trimmed,FileName0)
% 
% 
% %Rotate atoms back
% atoms_Final=atoms_Del_Rot_Shifted_Stretched*inv(B);
% %FileName0='Del_Rot_Sifted_Stretched_RotB.dump';
% %DumpMaker(savedir,atoms_Final,FileName0)














