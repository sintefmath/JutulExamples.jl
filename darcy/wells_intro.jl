using JutulDarcy, JutulViz, Jutul
## Define and plot the mesh
nx = 20
ny = 10
nz = 4
# Some useful constants
day = 3600*24
bar = 1e5
# Create and plot the mesh
dims = (nx, ny, nz)
g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
plot_mesh(g)
## Create a layered permeability field
Darcy = 9.869232667160130e-13
nlayer = nx*ny
K = vcat(repeat([0.65], nlayer), repeat([0.3], nlayer), repeat([0.5], nlayer), repeat([0.2], nlayer))*Darcy
## Set up a vertical well in the first corner, perforated in all layers
P = setup_vertical_well(g, K, 1, 1, name = :Producer);
## Set up an injector in the upper left corner
I = setup_well(g, K, [(nx, ny, 1)], name = :Injector);
## Plot the permeability (scaled to Darcy) and the wells
fig, ax, p = plot_cell_data(g, K/Darcy)
plot_well!(ax, g, I, textscale = 0.1)
plot_well!(ax, g, P, color = :darkblue, textscale = 0.1)
## Set up a two-phase immiscible system and define a density secondary variable
phases = (LiquidPhase(), VaporPhase())
sys = ImmiscibleSystem(phases)
rhoLS = 1000.0
rhoGS = 100.0
rhoS = [rhoLS, rhoGS]
c = [1e-6/bar, 1e-4/bar]
ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
model, parameters = setup_reservoir_model(g, sys, wells = [I, P], reference_densities = rhoS)
display(model)
## Replace the density function with our custom version
replace_variables!(model, PhaseMassDensities = ρ)
## Set up initial state
state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
## Set up time-steps
dt = repeat([30.0]*day, 12*5)
## Inject a full pore-volume (at reference conditions) of gas
# We first define an injection rate
reservoir = reservoir_model(model);
pv = pore_volume(model)
inj_rate = sum(pv)/sum(dt)
# We then set up a total rate target (positive value for injection)
# together with a corresponding injection control that specifies the
# mass fractions of the two components/phases for pure gas injection,
# with surface density given by the known gas density.
rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
# The producer operates at a fixed bottom hole pressure
bhp_target = BottomHolePressureTarget(50*bar)
P_ctrl = ProducerControl(bhp_target)
# Set up the controls. One control per well in the Facility.
controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl
# Set up forces for the whole model. For this example, all forces are defaulted
# (amounting to no-flow for the reservoir).
forces = setup_reservoir_forces(model, control = controls)
## Finally simulate!
sim, config = setup_reservoir_simulator(model, state0, parameters, info_level = 0)
states, reports = simulate!(sim, dt, forces = forces, config = config);

## Once the simulation is done, we can plot the states
f, = plot_interactive(g, map(x -> x[:Reservoir], states))
display(f)
## Plot the wells
wd = full_well_outputs(sim.model, parameters, states)
time = report_times(reports)
plot_well_results(wd, time)