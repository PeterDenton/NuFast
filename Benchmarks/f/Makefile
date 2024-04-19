MAKEFLAGS=--jobs=1

Compiler=gfortran

FFlags=-pedantic -Wall
FFlags+=-O5 -ffast-math # can improve speed by 10s of percent over the default of -O0
FFlags+=-fdefault-real-8 # set the precision of reals. commenting this out sets it to 4 (single precision: fastest and enough precision), 8 is double precision (best), 16 is quadruple precision (much slower, but fun to check)
Sources=src/Benchmark.f90 src/NuFast.f90 src/Exact_Cubic.f90
AllObjects=$(addprefix obj/,$(notdir $(Sources:.f90=.o)))
Executables=Benchmark
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	$(Compiler) $(Objects) $(addprefix obj/,$@.o) -o $@

obj/%.o: src/%.f90
	@mkdir -p $(@D)
	$(Compiler) $< -c $(FFlags) -Jobj -o $@

clean:
	rm -rf obj/*.o obj/*.mod $(Executables)

