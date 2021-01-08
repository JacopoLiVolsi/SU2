cpp_args = '-Wall'

extract_mesh_src = files(['main.cpp','Point.cpp','Mesh_Extractor.cpp'])

extract_mesh_incl = include_directories('../include','../../Common/include/')

su2_dep = include_directories('../../Common/include/')

extract_mesh_dep = declare_dependency(include_directories: su2_dep)

if get_option('enable-normal')
extract_mesh = executable('mesh_extractor',
                      extract_mesh_src,
                      install : true,
                      include_directories : extract_mesh_incl,
                      dependencies : [extract_mesh_dep, su2_deps, common_dep],
                      cpp_args:  ['-fPIC'] + [cpp_args])
endif


