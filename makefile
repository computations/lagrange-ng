deps_dir = $(shell pwd)/deps
extra_flags = -DCMAKE_EXPORT_COMPILE_COMMANDS=YES

.PHONY: all clean tests debug-flags release-flags install-prefix-flags build

all: release

build:
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
	rm -rf build bin
