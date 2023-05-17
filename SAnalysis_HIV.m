function SAnalysis_HIV
% SAnalysis_HBV Computes all stationary solutions of HBV model and analyzes their stability
%
% Latest revision 08.01.2022 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V.  Sklyarova     (MIPT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating working path
[filepath,~,~] = fileparts(mfilename('fullpath'));
cd(filepath);
cd ..
cd ..
addpath(genpath(pwd));
cd(filepath);
%
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING PARAMETERS:
%
% I.  Parameters for stationary solutions computation:
%
tol=1e-6; % required accuracy for computing the stationary solution
%
% II. Parameters for stationary solutions stability analysis:
%
p=7; % required number of eigenvalues
%
delta=1.0e-2; % grid step for the scheme to use for the approximation
%               of nonlinear eigenvalue problem 
%
scheme=2; % integration scheme to use for approximation
%           of nonlinear eigenvalue problem:
%               1 - implicit Euler method,
%               2 - BDF2. 
%
% Parameters of the succesive linear problems method:
R=10.0; % radius of a circle |z-z0|<=R on the complex plane
%         beyond which the calculated eigenvalue z must not escape. 
%
lmaxlp=10; % maximum number of iterations of the method
%
resmin=1.0e-12; % accuracy of computaion of eigenvalues
%
deflation=0; % deflation regime (-1,0,1,2 or 3) 
%              0 - no deflation, 
%              1, 2, 3 - deflation regimes, described in Demyanko K.V., 
%              Nechepurenko Yu.M. and Sadkane M. A Newton-type method for
%              non-linear eigenproblems // Russian Journal of Numerical
%              Analysis and Mathematical Modelling, vol. 32, no. 4, 2017,
%              pp. 237-244. This regimes differ in the nonlinear
%              eigenvalue problem which is solved;
%
addinf=0; % indicates whether to write step-by-step information on how
%           the partial eigenproblem solver works (0/1)
%
ind=0; % indicates whether to draw eigenvalues on the complex plane (0/1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters:
[coefs,delays]=GetPar_HIV1;
% File name for saving information:
resFile='SS_HIV_file';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. COMPUTATION OF ALL STATIONARY SOLUTIONS:
SSFinder = @SSFinder_HIV; % function to compute steady states with Wolfram 
                          % Mathematica (was generated for this model by
                          % NewSSFinderGeneration.m)
[SS,~]=SComputation(coefs, @RHS_HIV, SSFinder, tol);
%Adding trivial stationary solution
% U0= [0; 0; 0; 0; coefs.H_E_0; coefs.E_0; coefs.H_B_0; coefs.B_0; coefs.P_0; ...
%     coefs.rho_F*coefs.P_0/coefs.alpha_F];
% SS=[SS,U0];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % II. STABILITY ANALYSIS. COMPUTATIONS OF EIGENVALUES:
% %
% EiVals=[];
% EiVects=[];
for k=1:size(SS,2)
    fprintf(1,'For stationary solution %d: \n',k);
    Ui=SS(:,k);
    disp('Stationary solution:');
    disp(Ui);
    tic
end
    [Ei,Xi]=SStability(Ui,coefs,delays,@RHS_HIV,p,delta,...
                                  scheme,R,lmaxlp,resmin,deflation,addinf);
    fprintf(1,'%d leading eigenvalues: \n',p);
    disp(Ei);
    toc
%     disp(' ');
%     EiVals=[EiVals,Ei];
%     EiVects=[EiVects,Xi];
%     % Drawing eigenvalues:
%     if ind==1
%         figure(ifig)
%             plot(real(Ei),imag(Ei),'bx')
%             title('$\Lambda$','Interpreter','latex'); 
%             hold on
%     end
%     save (resFile,'SS','EiVals','EiVects', 'coefs', 'delays');
end
%save (resFile,'SS','EiVals','EiVects', 'coefs', 'delays');
