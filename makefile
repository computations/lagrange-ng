deps_dir = $(shell pwd)/deps
blis = $(deps_dir)/blis
blas_libs_prefix = $(deps_dir)/blas
extra_flags = -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -DBLAS_LIBS_INSTALL_PREFIX=$(blas_libs_prefix)

.PHONY: all clean tests debug-flags release-flags install-prefix-flags build

all: release

build:  $(blas_libs_prefix)
	@cmake -Bbuild -H. $(extra_flags)

tests: build
	@cd build && make lagrange-test
	./bin/lagrange-test

release: extra_flags += -DCMAKE_BUILD_TYPE=Release
release: build
	@cd build && $(MAKE)

debug: extra_flags += -DCMAKE_BUILD_TYPE=Debug
debug: build
	@cd build && $(MAKE)

clean:
	rm -rf build bin $(blas_libs_prefix)


$(blas_libs_prefix):
	cd $(deps_dir)/blis && ./configure --prefix=$(blas_libs_prefix) --disable-mixed-dt --enable-threading=openmp auto && $(MAKE) && $(MAKE) install
