%% Example - Rossi et al.'s Criminal Recidivism Data
% see: https://vincentarelbundock.github.io/Rdatasets/doc/carData/Rossi.html
clear;

load data/rossi;

%% Variables
vars_cat = [grp2idx(Rossi.fin), grp2idx(Rossi.wexp), grp2idx(Rossi.mar), grp2idx(Rossi.paro)];
vars_count = [Rossi.prio];
vars_continuous = [Rossi.age];

target = [Rossi.week Rossi.arrest];

data = [vars_cat, vars_count, vars_continuous, target];

%% Fit mixture model
VarNames = {'Fin. Aid','Work','Married','Parole','#Conv','Age','Week','Arrested'};
mm = snob(data, {'multi',1:4, 'negb',5,'norm',6,'cfixexp',7:8}, 'k', 1, 'varnames', VarNames);
mm_Summary(mm);

% Print KL divergence matrices
mm_KLstats(mm, data);
