using JutulDarcy, Jutul, JutulViz
# Example of the SPE9 model exported from MRST running in JutulDarcy.
#
#   Killough, J. E. 1995. Ninth SPE comparative solution project: A
#   reexamination of black-oil simulation. In SPE Reservoir Simulation
#   Symposium,  12-15 February 1995, San Antonio, Texas. SPE 29110-MS, doi:
#   10.2118/29110-MS
#
# For comparison against other simulators, see the equivialent JutulSPE9 example 
# in the Jutul module for MRST (www.mrst.no)
f = "data/spe9/spe9.mat"
# Add some more output properties for plotting
extras = [:Saturations, :Rs]
states, reports, pth, setup = simulate_mrst_case(f, extra_outputs = extras, precond = :ilu0);
## Plot reservoir
using JutulViz
mrst_data = setup[:mrst]
g = MRSTWrapMesh(mrst_data["G"])
W = mrst_data["schedule"]["control"][1]["W"]
plot_reservoir(g, states, wells = W)
## Plot the wells
wd = full_well_outputs(setup.sim.model, states, setup.case.forces)
time = report_times(reports)
plot_well_results(wd, time)
