using JutulDarcy, Jutul
# Example demonstrating the SPE1 benchmark case.
#
# Odeh, A.S. 1981. Comparison of Solutions to a Three-Dimensional Black-Oil
# Reservoir Simulation Problem. J Pet Technol 33 (1): 13–25. SPE-9723-PA.
# http://dx.doi.org/10.2118/9723-PA
#
# For comparison against other simulators, see the equivialent JutulSPE1 example 
# in the Jutul module for MRST (www.mrst.no)
states, reports, pth, setup = simulate_mrst_case("data/spe1/spe1.mat", extra_outputs = [:Rs, :Saturations]);
mrst_data = setup.mrst
## Plot the gas injection
using JutulViz
g = MRSTWrapMesh(mrst_data["G"])
W = mrst_data["schedule"]["control"][1]["W"]
plot_reservoir(g, states, wells = W)
## Plot well results
wd = full_well_outputs(setup.sim.model, states, setup.case.forces)
time = report_times(reports)
plot_well_results(wd, time)
