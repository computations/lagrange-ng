extra_flags = -DCMAKE_EXPORT_COMPILE_COMMANDS=YES

all: release

release:
	@cmake -Bbuild -H. $(extra_flags) && cd build && make

clean:
	rm -rf build bin