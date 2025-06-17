deps_dir = $(shell pwd)/deps
extra_flags = -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -G Ninja

.PHONY: all clean tests debug-flags release-flags install-prefix-flags build

all: release

build:
	@cmake -Bbuild -H. $(extra_flags)

tests: build
	@cd build && make lagrange-test
	./bin/lagrange-test

release: extra_flags += -DCMAKE_BUILD_TYPE=Release
release: build
	@cd build && ninja

debug: extra_flags += -DCMAKE_BUILD_TYPE=Debug
debug: build
	@cd build && ninja

clean:
	rm -rf build bin
