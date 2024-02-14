# JutulExamples
> [!CAUTION]
> This repository is deprecated. Examples can now be found in the examples folder of the respective repositories mentioned here.


 
This repositories contained examples for [Jutul.jl](https://github.com/sintefmath/Jutul.jl/) and the porous media simulator [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl/).

To get started, you will have to add the necessary packages. Start by going to the `darcy` subfolder in your favorite terminal, launch a recent version of julia and run the following:
```julia
]          # Enter package mode
activate . # Activate the darcy subfolder environment
add Jutul JutulDarcy
dev https://github.com/sintefmath/JutulViz.jl.git # Not registered, manual URL
instantiate # This will take some time
# hit backspace, then
include("two_phase_buckley_leverett.jl")
```
Congratulations, you have run your first Jutul example. There are additional examples in the folder, as well as notebooks in both [Jupyter](https://jupyter.org/) and [Pluto.jl](https://github.com/fonsp/Pluto.jl) formats.
