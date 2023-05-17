function Ploter(t,U,mark,iclean,h,v,varname,logind) 
%
% Latest revision 11.10.2019 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V. Sklyarova      (MIPT)
%
% Plotting model solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%        t - time grid (1-by-[T/delta] matrix);
%        U - solution in the grid nodes (n-by-[T/delta] matrix);
%        mark - color and style of the line;
%        iclean - cleaning of figures;
%        v,h - number of figures in vertical and horisontal in subplots 
%              matrix;
%        varname - cell 1-by-n structure consisting of variable names;
%        logind - choosing between log scale (1) and ordinary scale (0).
%
% OUTPUT: none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
varname=struct2cell(varname);
for i=1:size(U,1)
subplot(v,h,i)
   if(iclean) 
       hold off
   end    
   if (logind)
       semilogy(t,U(i,:),mark);
       symlog2('y');
   else
       plot(t,U(i,:),mark)
   end
   xlabel('$t$','Interpreter','latex');           
   title(varname{i, 1},'Interpreter','latex'); 
   hold on
end

