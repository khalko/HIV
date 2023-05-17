function matlabFunction_new(F,filename,vars)
%
% Latest revision 09.10.2021 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V. Sklyarova      (MIPT)
%
% Transformation of symbolic function to matlab function written to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%        F - symbolic function;
%        filename - name of the file where matlab function is written;
%        vars - symbolic variables of symbolic function.
%
% OUTPUT:   none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
matlabFunction(F,'File',filename,'Vars',vars);
fid=fopen(filename,'r');
C=fscanf(fid,'%s');
i=1;
while (C(i)~=')')
    i=i+1;
    j=i;
end
for k=1:size(C,2)
    if C(k)=='%'
        q=k;
    end
end
D=[C(1:8),' ',C(9:j-1),',varargin',C(j:end)];
fid_new=fopen(filename,'w');
fprintf(fid_new,'%s\n %s\n',D(1:j+10), D(j+11:29+q));
p=29+q;
for i=30+q:size(D,2)
    if D(i)==';'
        fprintf(fid_new,'%s\n', D(p+1:i));
        p=i;
    end    
end
fclose(fid);
fclose(fid_new);