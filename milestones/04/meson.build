executable(
    'milestone04',
    'main.cpp',
    include_directories : [lib_incdirs],
    link_with : [lib],
    dependencies : [eigen, mpi]
)

configure_file(
    input: meson.current_source_dir() / 'lj54.xyz',
    output: 'lj54.xyz',
    copy: true
)