function dcdt = TconQD(t,c)
% Function for the differential equations that describe the tumor-inclusive case.
% This function takes in inputs of time span (t) and concentrations (c).
% This function returns as output the differential equations 
% for all of the compartments, to be solved by an ode23s function.

% Concentrations
cp = c(1); %Concentration of drug in the blood/plasma (nM)
ck = c(2); %Concentration of drug in the kidney (nM)
cli = c(3); %Concentration of drug in the liver(nM)
cs = c(4); %Concentration of drug in the spleen(nM)
cl = c(5); %Concentration of drug in the lung(nM)
co = c(6); %Concentration of drug in the other compartments (nM)
ct = c(7); %Concentration of drug in the tumor (nM)


%Volume parameters, V (mL)--based on research DON'T CHANGE
vp = 1.078; % blood
vk = 0.3674; % kidney
vli = 1.2078; % liver
vs = 0.077; % spleen
vl = 0.1606; % lungs
vo = 55.7572; % other
vt = 0.08; % tumor


% Volumetric Flow parameters, Q (mL/h)--based on research
qk = 85.7714; % kidney
qli = 132.8985; % liver
qs = 18.8509; % spleen
ql = 4.7127; % lungs
qo = 700.3091; % other
qt = 1.7664; % tumor (from online research)

qtotal = qk+qli+qs+ql+qo+qt;

% Partition Coefficients, R (unitless)--based on exp data
rk = 43; % kidney
rli = 53; % liver
rs = 95; % spleen
rl = 30; % lungs
ro = 1; % other (unknown)
rt = 55; % tumor (change this only to start)

% Drug removal
k1 = 0.000001; 



% Derivative equations
cblp = (1/vp)*((((qk*ck)/rk)+((qli*cli)/rli)+((qs*cs)/rs)+((ql*cl)/rl)+((qo*co)/ro)+((qt*ct)/rt)) - (qtotal)*(cp)); % blood
ckp = (qk/vk)*(cp - (ck/rk)); % kidney
clip = (qli/vli)*(cp - (cli/rli) - k1*cli); % liver 
csp = (qs/vs)*(cp - (cs/rs)); % spleen
clp = (ql/vl)*(cp - (cl/rl)); % lung
cop = (qo/vo)*(cp - (co/ro)); % other
ctp = (qt/vt)*(cp - (ct/rt)); % tumor


% Cumulative equation
dcdt = [cblp; ckp; clip; csp; clp; cop; ctp];  % Blood, Kidney, Liver, Spleen, Lung, Other, Tumor