using JutulDarcy, JutulViz, Jutul
nx = 20
ny = 10
nz = 4
day = 3600*24
bar = 1e5

dims = (nx, ny, nz)
g = CartesianMesh(dims, (20.0, 15.0, 5.0))
plot_mesh(g)
##
Darcy = 9.869232667160130e-13
nlayer = nx*ny
K = vcat(repeat([0.65], nlayer), repeat([0.3], nlayer), repeat([0.5], nlayer), repeat([0.2], nlayer))*Darcy
##
I = setup_well(g, K, [(1, 1, 1)], name = :Injector);
##
P = setup_vertical_well(g, K, nx, ny, name = :Producer);
##
fig, ax, p = plot_cell_data(g, K/Darcy)
plot_well!(ax, g, I, textscale = 0.1)
plot_well!(ax, g, P, color = :darkblue, textscale = 0.1)
##
phases = (LiquidPhase(), VaporPhase())
sys = ImmiscibleSystem(phases)
rhoLS = 1000.0
rhoGS = 100.0
model, parameters = setup_reservoir_model(g, sys, wells = [I, P], reference_densities = [rhoLS, rhoGS]);
reservoir = reservoir_model(model)
##
state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
##
dt = repeat([30.0]*day, 12*5)
## Inject
pv = reservoir.domain.grid.pore_volumes
inj_rate = sum(pv)/sum(dt)
rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)

bhp_target = BottomHolePressureTarget(50*bar)
P_ctrl = ProducerControl(bhp_target)

controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl

facility = model.models.Facility
surface_forces = setup_forces(facility, control = controls)
# Set up forces for the whole model. For this example, all forces are defaulted
# (amounting to no-flow for the reservoir).
forces = setup_forces(model, Facility = surface_forces)
##
push!(model.models.Reservoir.output_variables, :PhaseMassDensities)
sim, config = setup_reservoir_simulator(model, state0, parameters)
states, reports = simulate(sim, dt, forces = forces, config = config);
##
plot_interactive(g, map(x -> x[:Reservoir], states))