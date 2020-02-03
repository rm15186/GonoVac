This repository contains the following MATLAB files, which define and simulate an individual-based
dynamic transmission model of gonorrhoea infection and treatment in MSM. The model incorporates
ciprofloxacin-sensitive and resistant strains, a time-dependent sexual contact network, and allows
a flexible choice of complex (individual-based) treatment, testing and tracing options, including
strain-sensitive diagnostic tests with variable time delays.

MATLAB (https://uk.mathworks.com/products/matlab.html) is required to run this code. It has been tested
in version R2019a, but should work in other versions too.


Model files:
AMR_IBM.m          : main model code (class)
base_params.mat    : data file with preset model parameters (in struct called 'params')
run_simple.m       : a 'simple' script to initialise, run and modify parameters of a single model
subaxis.m          : a function to organise subplots better on a single figure
                     see https://uk.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
                     and subaxis_license.txt
parseArgs.m        : used by subaxis (see above)
                         
