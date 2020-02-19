# GonoVac

Adaptation of code written by Adam Zienkiewicz. Original: AMR_IBM

VacAMR_IBM3 adds 3 different vaccination strategies, the vac parameter given to the class constructor tells it which strategy to use:
[1,0,0] - childhood vaccination
[0,1,0] - vaccination at screening
[0,0,1] - vaccination on diagnosis

run_vac runs one round of the simulation with the selected vaccination strategy 
run_avg runs a set number of simulations and plots the average prevalence, doses of the vaccine, number of people vaccinated at any one time.

