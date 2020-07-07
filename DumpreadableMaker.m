function DumpreadableMaker(savedir, AtomsTotal,FileName)
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





line1=[' ' num2str(natoms) 'atoms'];
line2='1 atom types';
line3=[num2str(Xlo) ' ' num2str(Xhi) ' xlo xhi'];
line4=[num2str(Ylo) ' ' num2str(Yhi) ' ylo yhi'];
line5=[num2str(Zlo) ' ' num2str(Zhi) ' zlo zhi'];
line6='Atoms';


fprintf(fid,'\n\n%s\n\n', line1);
fprintf(fid,'%s\n\n', line2);
fprintf(fid,'%s\n', line3);
fprintf(fid,'%s\n', line4);
fprintf(fid,'%s\n\n', line5);
fprintf(fid,'%s\n\n', line6);


fprintf(fid,'%d %d %f %f %f \n',AtomsTotal0')

fclose(fid);




        
end
