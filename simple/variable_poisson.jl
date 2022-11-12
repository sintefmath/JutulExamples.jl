using Jutul
function solve_poisson(nx, ny)
    sys = VariablePoissonSystem()
    # Unit square
    g = CartesianMesh((nx, ny), (1.0, 1.0))
    # Set up a model with the grid and system
    discretization = (poisson = Jutul.PoissonDiscretization(g), )
    D = DiscretizedDomain(g, discretization)
    model = SimulationModel(D, sys)
    # Initial condition doesn't matter
    state0 = setup_state(model, Dict(:U=>1.0))
    param = setup_parameters(model, K = compute_face_trans(g, 1.0))
    nc = number_of_cells(g)
    pos_src = PoissonSource(1, 1.0)
    neg_src = PoissonSource(nc, -1.0)
    forces = setup_forces(model, sources = [pos_src, neg_src])
    states, = simulate(state0, model, [1.0], forces = forces, 
                                             parameters = param)
    return (states, g)
end
nx = ny = 100
states, g = solve_poisson(nx, ny)
using JutulViz
plot_interactive(g, states)
