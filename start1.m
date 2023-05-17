Tadd = 150;
delta = 0.5e-3;
[coef, delays] = GetPar_HIV1;
% names = fieldnames(delays);
% [val, idx] = max(struct2array(delays));
% max_fieldname = names(idx);
abc = struct2cell(delays); %% 
abc=cell2mat(abc);
[abc,ord] = sort(abc);
m_q = abc(end)/delta;

M=ceil(Tadd/delta);
t_V = 0.01;
m_V = ceil(t_V/delta);
t =(-m_q:1:0)*delta;
U = zeros(7,size(t,2));
err = 1e-12;
interp1_method = 'linear';

U(1,1:m_q+1)=coef.r_A/coef.mu_A;
U(5,1:m_q+1)=coef.r_E0/coef.mu_E0;
U(4,m_q-m_V+1:m_q+1)=5000*(t(m_q-m_V+1:m_q+1)+0.01);

iresult = 1;
nfig = 0;
logind = false;
[varnames, ~] = LatexNames_HIV;
scheme = 1;
inc = 500000;
v = 3;
h = 3;
marks={'r-','b-'}; 
iclean = true;



