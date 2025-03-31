% MATLAB script to launch a fiber structure generation and the corresponding LBM simulation
%
%INPUT VARIABLES:
%
% SEED: integer representing the seed for initializing the random
% generator. If seed=0, automatic seed generation. If you want to reproduce
% the same fiber structure, use the same seed (fibers will be located at the same place). 
%
% MEAN_D: contains the mean fiber to be used
%
% STD_D: contains the standard deviation of the fiber diameters
%
% PORO: estimated porosity of the fiber structure to be generated
% 
% NX: domain lateral size in grid cell

% Etude de convergence u_num
L = 0.0002;
dx= 2e-6 ; % grid size in m
n_raf = 2;

mean_poro = 0.9;
std_poro = 7.5e-3;
poro_pdf=normrnd(mean_poro,std_poro,[1,1000]);

results = table();

seed=0;
for i = 1:length(poro_pdf)
    iter = i
    deltaP= 0.1 ; % pressure drop in Pa
    dx_i = dx/n_raf;
    NX= 100*n_raf;
    poro= poro_pdf(i) ;
    mean_fiber_d= 12.5 ; % in microns
    std_d= 2.85 ; % in microns
    
    filename= 'fiber_mat.tiff' ;
    
    % generation of the fiber structure
    [d_equivalent]=Generate_sample(seed,filename,mean_fiber_d,std_d,poro,NX,dx_i);
    
    % calculation of the flow field and the permeability from Darcy Law
    [poro_eff, Re, k_in_micron2] = LBM(filename,NX,deltaP,dx_i,d_equivalent,seed);
    results = [results; table(d_equivalent, dx_i, poro_eff, Re, k_in_micron2)];
    
    results_output = strcat('results_seed_', num2str(seed), '.csv');
    writetable(results, results_output);
end








