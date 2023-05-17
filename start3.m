Tadd = 150;
delta = 0.5e-3;
[coef, delays] = GetPar_HIV3;
% names = fieldnames(delays);
% [val, idx] = max(struct2array(delays));
% max_fieldname = names(idx);
abc = struct2cell(delays); %% 
abc=cell2mat(abc);
[abc,ord] = sort(abc);
m_q = abc(end)/delta;
%m_q = 5000;
M=ceil(Tadd/delta);
t_V = 0.01;
m_V = ceil(t_V/delta);
t =(-m_q:1:0)*delta;
U = zeros(7,size(t,2));
err = 1e-7;
interp1_method = 'spline';
%%%%22222
% U(1,1:m_q+1)=exp(12.2058)-1;
% U(2,1:m_q+1)=exp(2)-1;
% U(3,1:m_q+1)=exp(0.5872)-1;
% U(4,1:m_q-m_V)=exp(3.3452)-1;
% U(5,1:m_q+1)=exp(9.1195)-1;
% U(6,1:m_q+1)=exp(9.5289)-1;
% 
% U(4,m_q-m_V+1:m_q+1)=exp(3.3452)-1+15000*(t(m_q-m_V+1:m_q+1)+0.01);
% U(7,m_q:m_q+1)=exp(9.3488)-1;


%%%%1111
U(1,1:m_q+1)=coef.r_A/coef.mu_A;
U(5,1:m_q+1)=coef.r_E0/coef.mu_E0;
U(4,m_q-m_V+1:m_q+1)=5000*(t(m_q-m_V+1:m_q+1)+0.01);

%%%%3333
% U(1,1:m_q+1)=exp(12.2058)-1;
% U(2,1:m_q+1)=exp(2)-1;
% U(3,1:m_q+1)=exp(0.5872)-1;
% U(4,1:m_q-m_V)=exp(3.3452)-1;
% U(5,1:m_q+1)=exp(9.1195)-1;
% U(6,1:m_q+1)=exp(9.5289)-1;
% 
% U(4,m_q-m_V+1:m_q+1)=exp(3.3452)-1+15000*(t(m_q-m_V+1:m_q+1)+0.01);
% U(7, :)=1:size(U, 2);

iresult = 1;
nfig = 0;
logind = false;
[varnames, ~] = LatexNames_HIV;
scheme = 1;
inc = 500000;
v = 4;
h = 3;
marks={'r-','b-'}; 
iclean = true;
                 