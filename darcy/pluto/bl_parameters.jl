### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e25f84e0-b7e4-11ec-215a-15d61473a706
begin
	import Pkg
    # careful: this is _not_ a reproducible environment
    # activate local environment where dev-ed packages are added
    Pkg.activate(".")

	using Jutul
	using JutulDarcy
	using JutulViz
end

# ╔═╡ 6ac645bd-fcff-4923-ae5d-a51ef19e8817
Pkg.add("PlutoUI"); using PlutoUI

# ╔═╡ 5f179657-8ec1-43f0-8f16-2f08d87f7ad2
Pkg.add("Plots"); using Plots

# ╔═╡ 931359ca-7d00-4d4e-b9f1-fce440c3d3b9
md"""
Set up BL scenario
"""

# ╔═╡ 963b3405-b02c-4190-8078-2edf15a0a03e
md"""
##### Buckley-Leverett problem
Conservation of mass:

$$\frac{\partial \rho_\alpha S_\alpha \phi}{\partial t} + \nabla \cdot v_\alpha = q_\alpha, \alpha \in \{w, o\}$$

Fluxes are given by Darcy's law:
$$v_t = -K (\lambda_w + \lambda_o) \nabla p, \quad v_\alpha = \frac{\lambda_\alpha}{\lambda_w + \lambda_o} v_t$$

We solve the Buckley-Levrett problem with $\lambda_o = S^N/\mu$ and $\lambda_w = S^N$. The displacement can be stable or unstable depending on the viscosity ratio.
"""

# ╔═╡ d6a3bac0-b5dd-4fcf-bf36-afbd57bc7be7
@bind values PlutoUI.combine() do Child
md"""
$(Child("nx", Slider(10:1000, default = 1000))) Number of grid cells

$(Child("nkr", Slider(1.0:0.1:4, default=2))) Corey exponent

$(Child("sor", Slider(0:0.02:0.9, default=0.0))) Residual oil saturation

$(Child("mu", Slider(0.1:0.1:5, default = 1))) Viscosity ratio

$(Child("pvi", Slider(0:0.01:2, default = 0.5))) Pore volumes injected 

"""
end

# ╔═╡ 7ce78aae-5ba4-41f0-be79-6eef9ef7da58
nx, nkr, sor, mu, pvi = values

# ╔═╡ b53d1c1f-9217-4658-8915-b9e2235e8870
g = CartesianMesh((nx, 1, 1), (1000.0, 1.0, 1.0))

# ╔═╡ 3d207ba8-23fa-4c65-8da0-c17756018143
begin
	K = repeat([9.869232667160130e-14], number_of_cells(g))
	P = setup_well(g, K, [(nx, 1, 1)], name = :Producer)
	I = setup_well(g, K, [(1, 1, 1)], name = :Injector)
	sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
	model, parameters = setup_reservoir_model(g, sys, wells = [I, P],
		                                              block_backend = false)
	kr = BrooksCoreyRelPerm(2, nkr, [0.0, sor])
	μ = ConstantVariables([1, mu].*1e-3)
	model = replace_variables!(model, RelativePermeabilities = kr,
									  PhaseViscosities = μ)
	nothing
end

# ╔═╡ 17d3061a-bd02-4379-b198-ebc9ac040496
begin
	state0 = setup_reservoir_state(model, Pressure = 150*1e5, 
		                                  Saturations = [0.0, 1.0])
	day = 3600*24.0
	T = 1000*day
	nt = 50
	dt = repeat([T/nt], nt)
	reservoir = reservoir_model(model);
	pv = reservoir.domain.grid.pore_volumes
	inj_rate = pvi*sum(pv)/sum(dt)
	
	rate_target = TotalRateTarget(inj_rate)
	I_ctrl = InjectorControl(rate_target, [1.0, 0.0], density = 1000.0)
	# The producer operates at a fixed bottom hole pressure
	bhp_target = BottomHolePressureTarget(50*1e5)
	P_ctrl = ProducerControl(bhp_target)
	# Set up the controls. One control per well in the Facility.
	controls = Dict()
	controls[:Injector] = I_ctrl
	controls[:Producer] = P_ctrl
	facility = model.models.Facility
	surface_forces = setup_forces(facility, control = controls)
	# Set up forces for the whole model. For this example, all forces are defaulted
	# (amounting to no-flow for the reservoir).
	forces = setup_forces(model, Facility = surface_forces)
end

# ╔═╡ 7f3d7a62-65e9-4d64-add7-f0fb81669327
begin
		sim, config = setup_reservoir_simulator(model, state0, parameters,
			                                    info_level = -1)
		states, reports = simulate(sim, dt, forces = forces,
			                                config = config);
		s_w = hcat(map(x -> x[:Reservoir][:Saturations][1, :], states)...);
end;

# ╔═╡ 265591ef-549c-4c69-a585-690ab0965203
plot_title = "$nx cells, Corey exp $nkr, viscosity ratio $mu, pvi $pvi"

# ╔═╡ d667ebe9-7b1e-4066-bc5c-8f92b570a5b4
md"""
plot type: $(@bind plot_type Select(["lines", "contour"]))
"""

# ╔═╡ 21b7dd5c-ef28-4729-bbdf-62fc67a9b91c
if plot_type == "lines"
	plot(s_w, color = :black, title = plot_title, label = nothing, ylim=(0,1))
else
	contourf(s_w, title = plot_title)
end

# ╔═╡ Cell order:
# ╟─e25f84e0-b7e4-11ec-215a-15d61473a706
# ╟─6ac645bd-fcff-4923-ae5d-a51ef19e8817
# ╟─5f179657-8ec1-43f0-8f16-2f08d87f7ad2
# ╠═b53d1c1f-9217-4658-8915-b9e2235e8870
# ╠═3d207ba8-23fa-4c65-8da0-c17756018143
# ╟─931359ca-7d00-4d4e-b9f1-fce440c3d3b9
# ╠═17d3061a-bd02-4379-b198-ebc9ac040496
# ╠═7f3d7a62-65e9-4d64-add7-f0fb81669327
# ╠═265591ef-549c-4c69-a585-690ab0965203
# ╟─963b3405-b02c-4190-8078-2edf15a0a03e
# ╟─d6a3bac0-b5dd-4fcf-bf36-afbd57bc7be7
# ╟─7ce78aae-5ba4-41f0-be79-6eef9ef7da58
# ╟─d667ebe9-7b1e-4066-bc5c-8f92b570a5b4
# ╟─21b7dd5c-ef28-4729-bbdf-62fc67a9b91c
