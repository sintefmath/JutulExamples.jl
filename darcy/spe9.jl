using JutulDarcy, Jutul, JutulViz
f = "data/spe9/spe9.mat"
states, reports, pth, setup = simulate_mrst_case(f);
##
wd = full_well_outputs(setup.sim.model, states, setup.forces)
time = report_times(reports)
plot_well_results(wd, time)
