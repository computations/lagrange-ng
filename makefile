extra_flags = -DCMAKE_EXPORT_COMPILE_COMMANDS=YES

.PHONY: all clean tests

all: release

build:
	@cmake -Bbuild -H. $(extra_flags)

tests: build
	@cd build && make lagrange-test
	./bin/lagrange-test

release: build
	@cd build && make

clean:
	rm -rf build bin
