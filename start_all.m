[tnew,Unew, desp_array]=Solver(coef,delays,@RHS_HIV,t,U,Tadd,delta,err,...
                            interp1_method,iresult,nfig,iclean,h,v,marks,...
                            inc,scheme,varnames,logind);
                        

%
% Forming a MATLAB function that computes the right-hand side of the system
%

