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

# ╔═╡ 963b3405-b02c-4190-8078-2edf15a0a03e
md"""
##### LaTeX
We solve the Buckley-Levrett problem with $\lambda_o = S^N/\mu$ and $\lambda_w = S^N$
"""

# ╔═╡ 9ef7e5f3-e6f4-43a0-baeb-6a0a1d667e67
@bind nx Slider(10:1000, default = 1000)

# ╔═╡ b53d1c1f-9217-4658-8915-b9e2235e8870
g = CartesianMesh((nx, 1, 1), (1000.0, 1.0, 1.0))

# ╔═╡ 958b0b3e-032b-4e42-af0a-60d1a7521fe3
@bind nkr Slider(1:0.1:4, default=2)

# ╔═╡ 459c6cac-02f6-4dc5-a60f-62ae4e120fe4
@bind mu Slider(0.1:0.1:5, default = 1)

# ╔═╡ 3d207ba8-23fa-4c65-8da0-c17756018143
begin
	K = repeat([9.869232667160130e-14], number_of_cells(g))
	P = setup_well(g, K, [(nx, 1, 1)], name = :Producer)
	I = setup_well(g, K, [(1, 1, 1)], name = :Injector)
	sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
	model, parameters = setup_reservoir_model(g, sys, wells = [I, P], block_backend = false)
	kr = BrooksCoreyRelPerm(2, nkr)
	μ = ConstantVariables([1, mu].*1e-3)
	replace_variables!(model, RelativePermeabilities = kr, PhaseViscosities = μ)

end

# ╔═╡ 265591ef-549c-4c69-a585-690ab0965203
plot_title = "$nx cells, Corey exp $nkr, viscosity ratio $mu"

# ╔═╡ d4b5ade0-ae35-451f-999d-28105306ee3c
@bind pvi Slider(0:0.01:2, default = 0.5)

# ╔═╡ 17d3061a-bd02-4379-b198-ebc9ac040496
begin
	state0 = setup_reservoir_state(model, Pressure = 150*1e5, Saturations = [0.0, 1.0])
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
		sim, config = setup_reservoir_simulator(model, state0, parameters, info_level = -1)
		states, reports = simulate(sim, dt, forces = forces, config = config);
		s_w = hcat(map(x -> x[:Reservoir][:Saturations][1, :], states)...)
end

# ╔═╡ 21b7dd5c-ef28-4729-bbdf-62fc67a9b91c
plot(s_w, color = :black, title = plot_title, label = nothing)

# ╔═╡ 940cfa2a-868c-42da-9039-b0ceef35ed73
contourf(s_w, title = plot_title)

# ╔═╡ Cell order:
# ╠═e25f84e0-b7e4-11ec-215a-15d61473a706
# ╠═6ac645bd-fcff-4923-ae5d-a51ef19e8817
# ╠═5f179657-8ec1-43f0-8f16-2f08d87f7ad2
# ╠═b53d1c1f-9217-4658-8915-b9e2235e8870
# ╠═3d207ba8-23fa-4c65-8da0-c17756018143
# ╠═17d3061a-bd02-4379-b198-ebc9ac040496
# ╠═7f3d7a62-65e9-4d64-add7-f0fb81669327
# ╠═265591ef-549c-4c69-a585-690ab0965203
# ╠═963b3405-b02c-4190-8078-2edf15a0a03e
# ╠═21b7dd5c-ef28-4729-bbdf-62fc67a9b91c
# ╠═940cfa2a-868c-42da-9039-b0ceef35ed73
# ╠═9ef7e5f3-e6f4-43a0-baeb-6a0a1d667e67
# ╠═958b0b3e-032b-4e42-af0a-60d1a7521fe3
# ╠═459c6cac-02f6-4dc5-a60f-62ae4e120fe4
# ╠═d4b5ade0-ae35-451f-999d-28105306ee3c
