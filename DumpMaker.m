function DumpMaker(savedir, AtomsTotal,FileName)
[natoms,col] = size(AtomsTotal);

AtomsTotal0(:,1)=1:natoms;
AtomsTotal0(:,2)=1;
AtomsTotal0(:,3:5)=AtomsTotal(:,1:3);
%Dumpfile maker
fid = fopen([savedir FileName], 'w');



Extention=20;
Xlo=min(AtomsTotal(:,1))-Extention;
Xhi=max(AtomsTotal(:,1))+Extention;
Ylo=min(AtomsTotal(:,2))-Extention;
Yhi=max(AtomsTotal(:,2))+Extention;
Zlo=min(AtomsTotal(:,3))-Extention;
Zhi=max(AtomsTotal(:,3))+Extention;


line1='ITEM: TIMESTEP';
line2='1';
line3='ITEM: NUMBER OF ATOMS';
line4=num2str(natoms);
line5='ITEM: BOX BOUNDS pp pp pp';
line6=[num2str(Xlo) ' ' num2str(Xhi)];
line7=[num2str(Ylo) ' ' num2str(Yhi)];
line8=[num2str(Zlo) ' ' num2str(Zhi)];
line9='ITEM: ATOMS id type x y z';


fprintf(fid,'%s\n', line1);
fprintf(fid,'%s\n', line2);
fprintf(fid,'%s\n', line3);
fprintf(fid,'%s\n', line4);
fprintf(fid,'%s\n', line5);
fprintf(fid,'%s\n', line6);
fprintf(fid,'%s\n', line7);
fprintf(fid,'%s\n', line8);
fprintf(fid,'%s \n', line9);

fprintf(fid,'%d %d %f %f %f \n',AtomsTotal0')


fclose(fid);




        
end
