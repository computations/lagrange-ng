project('lagrange-ng', 'cpp',
  version: '0.1')

source_dir = 'src'
test_source_dir = join_paths('tests', 'src')
test_data_dir = join_paths('tests', 'data')

mkl_dynamic = dependency('mkl-dynamic-lp64-iomp', required:false)
mkl_static = dependency('mkl-static-lp64-iomp', required:false)

test_conf_data = configuration_data()
test_conf_data.set('RANDOM_TREE_FILE', join_paths(test_data_dir, 'random_test_trees'))
test_conf_data.set('PATHOLOGICAL_TREE_FILE', join_paths(test_data_dir,
'pathological_trees'))

configure_file(input : join_paths(test_source_dir, 'config.h.in'),
  output : 'test_data_paths.h',
  configuration : test_conf_data)

configuration_inc = include_directories('.')

if mkl_static.found()
  mkldep = mkl_static
else
  mkldep = mkl_dynamic
endif

cmake = import('cmake')

nloptdep = dependency('NLopt', fallback: ['nlopt', 'nlopt'])

jsondep = dependency('nlohmann_json', version: '>= 3.9.0',
  fallback : ['nlohmann_json'])

if get_option('buildtype') == 'debug'
  cpp_args = ['-Og', '-fsanitize=address']
  link_args = cpp_args
else
  cpp_args = []
  link_args = cpp_args
endif

lagrange_sources = files(
    join_paths(source_dir, 'AncSplit.cpp'),
    join_paths(source_dir, 'Context.cpp'),
    join_paths(source_dir, 'InputReader.cpp'),
    join_paths(source_dir, 'Node.cpp'),
    join_paths(source_dir, 'Operation.cpp'),
    join_paths(source_dir, 'Tree.cpp'),
    join_paths(source_dir, 'TreeReader.cpp'),
    join_paths(source_dir, 'Utils.cpp'),
    join_paths(source_dir, 'WorkerState.cpp'),
    join_paths(source_dir, 'Workspace.cpp'),
    )


lagrange = executable('lagrange-ng',
  [
    join_paths(source_dir, 'main.cpp'),
    lagrange_sources
  ],
  dependencies: [mkldep, nloptdep, jsondep],
  cpp_args: cpp_args,
  link_args: link_args,
  install: true,
  install_dir : join_paths(meson.project_source_root(),'bin'),
  )

gtest_dep = dependency('gtest', required: false, disabler:true)

lagrange_tests = executable('lagrange-test',
  [
    join_paths(test_source_dir, 'tree_reader.cpp'),
    join_paths(test_source_dir, 'workspace.cpp'),
    join_paths(test_source_dir, 'operation.cpp'),
    join_paths(test_source_dir, 'context.cpp'),
    join_paths(test_source_dir, 'tree.cpp'),
    join_paths(test_source_dir, 'main.cpp'),
    lagrange_sources
  ],
  dependencies: [gtest_dep, mkldep, nloptdep, jsondep],
  include_directories: [configuration_inc, source_dir],
  install: true,
  install_dir: join_paths(meson.project_source_root(),'bin'),
  )
