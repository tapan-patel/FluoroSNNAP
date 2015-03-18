function [sig,pval,F] = FC_granger(dF_cell)
load('params.mat');
X = dF_cell;
momax = params.FC.GC.morder; % Maximum model order
icregmode = 'LWR';
regmode = 'OLS';
mhtc = 'FDR';
alpha = params.FC.GC.alpha;
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
morder = 1;
[A,SIG] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed');

%%
[G,info] = var_to_autocov(A,SIG,params.FC.GC.iter);
var_info(info,true); 
F = autocov_to_pwcgc(G);

%%
assert(~isbad(F,false),'GC calculation failed');
pval = mvgc_pval(F,morder,size(dF_cell,2),1,1,1,size(dF_cell,1)-2,'F');
sig  = significance(pval,alpha,mhtc);
sig(isnan(sig))=0;
sig = sig';
pval = pval';
F = F';
