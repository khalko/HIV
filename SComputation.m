function [SS,Resid]=SComputation(coef, RHS, SSFinder, tol)
%SComputation Computes non negative stationary solutions to a system of DDEs
%
% Latest revision 16.11.2021 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V. Sklyarova      (MIPT)
%
% External functions: Residual, matlabFunction_new
%
% Internal functions: none
%
% Algorithm: ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
%        coef - model parameters excluding delays;
%
%        RHS  - a callback function specifying right hand side of the
%               system;
%        
%        SSFinder - a callback function for computation of the solutions of
%                   nonlinear system with Wolfram Mathematica;
%
%        tol - positive double scalar value specifying required 2-norm of
%              residual of computation of stationary solutions.
%
%
% OUTPUT:    
%         SS - N-by-S double matrix of computed stationary solutions
%              where N is the total number of the model variables and S
%              is the total number of the computed stationary
%              solutions;
%           
%         Resid - 1-by-S double row vector of 2-norm of residuals where
%                 S is the total number of the computed stationary
%                 solutions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating working path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(filepath));
%
format long
%
g=pwd;
D=dir('D:\13.1\');%%%%%%%%%%C:\Program Files\Wolfram Research\Mathematica\
D=struct2cell(D);
for j=1:size(D,2)
    addpath(['D:\13.1\',D{1,j}]);%%%%%%%%%%C:\Program Files\Wolfram Research\Mathematica\
end
mapa=which('math.exe');
mapa=mapa(1:end-8);
cd(mapa);
addpath(g);
SS=SSFinder(coef);
cd(g);
SS=SS';
Resid=zeros(1,size(SS,2));
for i=1:size(SS,2)
   if Residual(SS(:,i),coef, RHS)>tol
        l=0;
        [UZ,UD,F]=RHS(coef);
        UD=struct2cell(UD);
        UZn=cell(1,1);
        UZn{1,1}=UZ;
        Var=[UZn;UD];
        matlabFunction_new(F,'rhsm.m',Var);
        Y=rhsm(SS(:,i),SS(:,i),SS(:,i),SS(:,i),SS(:,i),SS(:,i))';
        normY=norm(Y,2);
        for j=2:size(Var,1)
            for k=1:size(Var{1,1},1)
                F=subs(F,Var{j,1}(k),Var{1,1}(k));
            end
        end
        Jacob=jacobian(F,Var{1,1});
        while(normY>tol)
            jacob1=subs(Jacob,Var{1,1},SS(:,i));
            jacob1=double(jacob1);
            SS(:,i)=SS(:,i)-jacob1\Y;
            l=l+1;
            Y=rhsm(SS(:,i),SS(:,i),SS(:,i),SS(:,i),SS(:,i),SS(:,i))';
            normY=norm(Y,2);
        end
   end
   Resid(i)=Residual(SS(:,i),coef, RHS);
end

  
