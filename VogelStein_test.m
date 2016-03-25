
load('/Users/Yves/Documents/MATLAB/2p_General/Template_fitting/TestGRABinfo')

%%

%701,639,689 are all great trace
F = GRABinfo.Traces(639,:);

% Set Parameters
% V.Ncells= 1;
% V.Npixels= 1;
% V.T= 13500;
% V.dt= 0.0333;
V.fast_iter_max= 10;
% V.fast_nonlin = 0;
% V.fast_poiss = 0;
% V.est_gam = 1;


P.gam = 0.9;
P.a = 0.8;
P.b = 0;
%P.sig = 0.3;


[n_best P_best V C]=fast_oopsi(F,V,P);

plot(F)
hold on
plot(round(n_best),'LineWidth',2)