function dcdt = conQDnoTumor(t,c)
% Function for the differential equations that describe the non-tumor case,
% untargeted ConQD case.
% This function takes in inputs of time span (t) and concentrations (c). 
% This function returns as output the differential equations for all of the 
% organs, to be solved by an ode23s function.

% Concentrations
cp = c(1); %Concentration of drug in the blood/plasma (nM)
ck = c(2); %Concentration of drug in the kidney (nM)
cli = c(3); %Concentration of drug in the liver(nM)
cs = c(4); %Concentration of drug in the spleen(nM)
cl = c(5); %Concentration of drug in the lung(nM)
co = c(6); %Concentration of drug in the other compartments (nM)


%Volume parameters, V (mL)--based on research
vp = 1.078; % blood
vk = 0.3674; % kidney
vli = 1.2078; % liver
vs = 0.077; % spleen
vl = 0.1606; % lungs
vo = 55.7572; % other


% Volumetric Flow parameters, Q (mL/h)--based on research 
qk = 85.7714; % kidney
qli = 132.8985; % liver
qs = 18.8509; % spleen
ql = 4.7127; % lungs
qo = 700.3091; % other

qtotal = qk+qli+qs+ql+qo;

% Partition Coefficients, R (unitless)
rk = 34; % kidney
rli = 55; % liver
rs = 40; % spleen
rl = 15; % lungs
ro = 1; % other (unknown)

% Drug removal
k1 = 0.000001; % Degradation rate from liver 

% Derivative equations
cblp = (1/vp)*((((qk*ck)/rk)+(((qli+qs)*cli)/rli)+((ql*cl)/rl)+((qo*co)/ro)) - (qtotal)*(cp));
ckp = (qk/vk)*(cp - (ck/rk));
clip = (1/vli)*(qli*cp + ((qs*cs)/rs) - (qli+qs)*(cli/rli) - k1*cli); 
csp = (qs/vs)*(cp - (cs/rs));
clp = (ql/vl)*(cp - (cl/rl));
cop = (qo/vo)*(cp - (co/ro));


% Cumulative equation
dcdt = [cblp; ckp; clip; csp; clp; cop];  % Blood, Kidney, Liver, Spleen, Lung, Other