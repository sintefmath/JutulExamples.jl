using JutulDarcy
states, reports, pth, setup = simulate_mrst_case("data/spe1/spe1.mat", extra_outputs = [:Rs, :Saturations]);
mrst_data = setup.mrst
## Plot the gas injection
using JutulViz
g = MRSTWrapMesh(mrst_data["G"])
W = mrst_data["schedule"]["control"][1]["W"]
plot_reservoir(g, states, wells = W)
## Plot well results
wd = full_well_outputs(setup.sim.model, setup.parameters, states)
time = report_times(reports)
plot_well_results(wd, time)