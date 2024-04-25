CXX=gcc

OMPFlags=-fopenmp

CFlags=-c -Wall -MMD -std=c99#-Werror -pedantic
#CFlags+=-O3 -ffast-math
CFlags+=$(OMPFlags)

GSLFlags=-L/home/pdenton/include/gsl/lib -lgsl -lgslcblas -lm

Sources=$(wildcard src/*.c)
IncludeDir=-I./include
IncludeDir+=-I/home/pdenton/include
IncludeDir+=-I/home/pdenton/include/gsl/include
AllObjects=$(addprefix obj/,$(notdir $(Sources:.c=.o)))
Executables=main
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data obj
	$(CXX) $(Objects) $(addprefix obj/,$@.o) $(GSLFlags) $(OMPFlags) -o $@

obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)

