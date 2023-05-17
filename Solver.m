function [tnew,Unew, desp_array]=Solver(coef,delays,RHS,t,U,Tadd,delta,err,...
                            interp1_method,iresult,nfig,iclean,h,v,marks,...
                            inc,scheme,varnames,logind)
%Solver Numerical integrates of the delay differential equation over time 
%
% Latest revision 20.01.2022 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V. Sklyarova      (MIPT)
%
% External functions: Ploter, matlabFunction_new
%
% Internal functions: Step, Res, Jacob
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%        coef - the system parameter values excluding delays;
%
%        delays - values of delays;
%
%        RHS  - a callback function specifying right hand side of the
%               system;
%
%        t - time grid;
%
%        U - the solution in the nodes of the grid t;
%
%        Tadd -  time of integration;
%
%        delta - grid step;
%
%        err - accuracy of Newton's method;
%
%        interp1_method - interpolation method: 'linear', 'spline',
%                         'pchip';
%
%        iresult=1 - supplement the old solution with a new computation,    
%                2 - to keep only the new computation,       
%                3 - keep only the "tail" of the new computation, which is 
%                    necessary to continue the integration, 
%                4 - save the new computation and the initial data for it;
%
%        nfig - number of figure (if nfig<1, there will be no figure);
%
%        iclean - images pre-cleaning;
%
%        marks - type of the plot lines;
% 
%        inc - screen output increment;
%
%        integration scheme: 1 - implicit Euler method,
%                            2 - BDF2.
%
%  OUTPUT:
%          tnew - new time grid;
%
%          Unew - solution in the nodes of the new grid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating working path


% displaying parameter values on the screen:
%
% fprintf(1,'T=%g\n',Tadd);
% fprintf(1,'delta=%g\n',delta);

% Plotting initial value:
%
% if (nfig>0)
%    figure(nfig+100)
% Ploter(t,U,marks{1},iclean,h,v,varnames,logind); 
% end 

%[filepath,~,~] = fileparts(mfilename('fullpath'));
format long
%addpath(genpath(filepath));
%
vtau = struct2cell(delays); %% 
vtau=cell2mat(vtau);
[vtau,ord] = sort(vtau);
[UZ,UD,F]=RHS(coef);
delnam=fieldnames(delays);
UD=orderfields(UD,delnam);
UD=struct2cell(UD);
UD=UD(ord);

UZn=cell(1,1);
UZn{1,1}=UZ;
Var=[UZn;UD];
%Var{6,1} = T_now;
matlabFunction_new(F,'rhsm.m',Var);


%
% Forming a MATLAB function that computes the Jacobi matrix
%
dGdU=jacobian(F,Var{1,1});
matlabFunction_new(dGdU,'jac.m',Var);
%
% Time grid for new computations:
%
vm=ceil(vtau/delta);
vm1=zeros(7,1);
vm1(1:size(vm,1))=vm;
mp=vm(end);
K=ceil(Tadd/delta);
desp_array = zeros(1, K);
%
if(K==0)
   disp('Integration time is less than one step')
   STOP
end    
tadd=t(end)+((1-mp):1:K)*delta; 
%
% Formation of the matrix of initial values
%
Uadd=zeros(size(Var{1,1},1),mp+K);
k=sum(t-t(end)<=-vtau(end));
if(k<1)
   disp('t(1)>t(end)-tau_max')
end
%
Uadd(:,1:mp)=interp1(t(k:end),U(:,k:end)',tadd(1:mp),interp1_method)';
%
% Nonlinear time integration:
%   
  
   nns=0; 
   for k=mp+1:mp+K
   %ab = mp+1;
   [Uadd(:,k),nnsk, desp]=Step(Uadd(:,k-1),Uadd(:,k-2),Uadd(:,k-vm1(1)),...
                                    Uadd(:,k-vm1(2)),Uadd(:,k-vm1(3)),...
                                    Uadd(:,k-vm1(4)),Uadd(:,k-vm1(5)),...
                                    Uadd(:,k-vm1(6)),Uadd(:,k-vm1(7)),...
                                    delta,err,scheme, mp, Uadd, coef, tadd, delays);
      nns=nns+nnsk;
      %desp_array(1, k-mp) = desp;
% % % % % % % % %       fprintf("%d\n", result); %%%%%%%%%%%%%%
%       fprintf('%g\n', nnsk);
      if (mod(k-mp,5e3)==0) % method information
         fprintf("%d \n", (k-mp)/K*100);
         
      end   
   end    
 Ploter(tadd,Uadd,marks{1},iclean,h,v,varnames,logind) 
% Drawing a new computation and its initial data 
%
% x = 1:K;
% scatter(x, desp_array)
if (nfig>0)
   figure(nfig)
 Ploter(tadd,Uadd,marks{1},iclean,h,v,varnames,logind) 
end

if(iresult==1)   
   tnew=[t,tadd(mp+1:mp+K)];
   Unew=[U,Uadd(:,mp+1:mp+K)]; 
elseif(iresult==2)       
   tnew=tadd(mp+1:mp+K); 
   Unew=Uadd(:,mp+1:mp+K);
elseif(iresult==3)     
   tnew=tadd(K:mp+K);
   Unew=Uadd(:,K:mp+K);
elseif(iresult==4)      
   tnew=tadd; 
   Unew=Uadd;   
end 
%clear t U tadd Uadd
delete jac.m
delete rhsm.m

function [U,l, normG]=Step(Udelta,Udelta2,U1,U2,U3,U4,U5,U6,U7,delta,err,scheme, k, Uadd, coef, tadd, delays)
%
% Latest revision 03.01.2017 
%
% The step of solving the Cauchy problem.
%
U=Udelta;
G=Res(U,Udelta,Udelta2,U1,U2,U3,U4,U5,U6,U7,delta,scheme, k, Uadd, coef, tadd, delays);
normG=norm(G,2);
err=err*norm(U,2);
l=0;
%desp_n = zeros(1, 15);
 
%fprintf("number of step %d and initial disp is %d\n", k, norm(G,2));
%min_discrepancy = zeros(1, 15);
while(normG>err &&l<15)%ограничить число шагов
   U=U-Jacob(U,U1,U2,U3,U4,U5,U6,U7,delta,scheme)\G;
   U(U<0)=0;
   G=Res(U,Udelta,Udelta2,U1,U2,U3,U4,U5,U6,U7,delta,scheme, k, Uadd, coef, tadd, delays);
   normG=norm(G,2);
   %desp_n(:, l+1) = normG;
   %min_discrepancy = normG; 
   %fprintf("discrepancy = %d\n", norm(G, 2));
   l=l+1;
   if ((k > 33345 && k < 33365))
       fprintf("k = %d, l = %d, desp = %d\n",k ,l,  normG);
       
   end
end

    

function G=Res(U,Udelta,Udelta2,U1,U2,U3,U4,U5,U6,U7,delta,scheme, k, Uadd, coef, tadd, delays)
%
% Latest revision 03.01.2017 
%
if(scheme==1)
   G=(U-Udelta)/delta; 
elseif(scheme==2)    
   G=(1.5*(U-Udelta)+0.5*(Udelta2-Udelta))/delta;
end
g=rhsm(U,U1,U2,U3,U4,U5,U6,U7);
b = g' + integral(k, Uadd, coef, tadd, delays, delta);
%fprintf('%g\n', g);
%asldajsl = integral(k, Uadd, coef, tadd, delays, delta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=G-b; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0;
% if(k == 5001)
%     fprintf('%g\n', g'+integral(k, Uadd, coef, tadd, delays, delta)); %%%%%%%% here we are 
% end
end
  
function dGdU=Jacob(U,U1,U2,U3,U4,U5,U6,U7,delta,scheme)
%
% Latest revision 12.10.2019 
%
if(scheme==1)
   dGdU=(1.0/delta)*eye(size(U,1)); 
elseif(scheme==2)
   dGdU=(1.5/delta)*eye(size(U,1));
end    
dGdU=dGdU-jac(U,U1,U2,U3,U4,U5,U6,U7);
end


function integral_result=integral(k, Uadd, coef, tadd, delays, delta)
delay4int = ceil(delays.w_C/delta);
A_int = Uadd(1,k-delay4int+1:k);
L_int = Uadd(2,k-delay4int+1:k);
I_int = Uadd(3,k-delay4int+1:k);
V_int = Uadd(4,k-delay4int+1:k);

A_int = [A_int Uadd(1,k)];
L_int = [L_int Uadd(2,k)];
I_int = [I_int Uadd(3,k)];
V_int = [V_int Uadd(4,k)];
%t_s = zeros(1, 401);
%t_s = tadd(k-delay4int:k);
t_s = delays.w_C:-delta:0;
%integral_result = zeros(7, 1);
%                                                                                                         %ytyj            
integral_r = trapz(delta,exp(-coef.mu_C*t_s).*((1-coef.p_L)*coef.g_AV*A_int.*V_int+coef.g_AI*A_int.*I_int+coef.g_LI*L_int.*I_int+coef.a_L*L_int)) ; %)); %%%% на 401 пусто а не должно
%trapz(delta,exp(-coef.mu_C*t_s).*((1-coef.p_L)*coef.g_AV*A_int.*V_int+coef.g_AI*A_int.*I_int+
%fprintf("integral equal %d\n", integral_r);
integral_r = -coef.g_CV * integral_r * V_int(end);
integral_result = zeros(7, 1);
integral_result(4, 1) = integral_r;


end    
end
end


% function dGdU=Jacob1(U,U1,U2,U3,U4,U5,U6,U7,delta,scheme)
% %
% % Latest revision 12.10.2019 
% %
% if(scheme==1)
%    dGdU=(1.0/delta)*eye(size(U,1)); 
% elseif(scheme==2)
%    dGdU=(1.5/delta)*eye(size(U,1));
% end    
% dGdU=dGdU-jac(U,U1,U2,U3,U4,U5,U6,U7);
    
    