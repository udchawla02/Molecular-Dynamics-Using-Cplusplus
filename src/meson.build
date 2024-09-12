lib_sources = [  # All source files (excluding headers)
    'hello.cpp', 'lj_direct_summation.cpp', 'verlet.cpp', 'xyz.cpp', 'energy.cpp', 'neighbors.cpp', 'simulation.cpp', 'simulation_data.cpp', 'ih.cpp', 'ducastelle.cpp', 'domain.cpp', 'export.cpp','th.cpp'
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen, mpi]
)

lib_incdirs = include_directories('.')
