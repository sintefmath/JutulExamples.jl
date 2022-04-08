using JutulDarcy, JutulViz, Jutul
nx = 20
ny = 10
nz = 4
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
model, parameters = setup_reservoir_model(g, sys, wells = [I, P], reference_densities = [1000.0, 100.0]);
##
setup_reservoir_state(model, Pressure = 3, Saturations = [1, 0])
##
