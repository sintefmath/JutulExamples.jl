using Jutul
nx = ny = 100
sys = SimpleHeatSystem()
# Unit square
g = CartesianMesh((nx, ny), (1.0, 1.0))
# Set up a model with the grid and system
D = DiscretizedDomain(g)
model = SimulationModel(D, sys)
# Initial condition is random values
nc = number_of_cells(g)
T0 = rand(nc)
T0 = zeros(nx, ny)

T0[nx÷4:3nx÷4, ny÷4:3ny÷4] .= 1

state0 = setup_state(model, Dict(:T=>vec(T0)))
Δt = repeat([1e-3], 10);
states, = simulate(state0, model, Δt)
pushfirst!(states, state0);

using JutulViz
plot_interactive(g, states)
