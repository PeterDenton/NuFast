CXX=gcc

OMPFlags=-fopenmp

CFlags=-c -Wall -MMD -std=c99
CFlags+=-O5 -ffast-math # can improve speed
CFlags+=$(OMPFlags)

GSLFlags=-lgsl -lgslcblas -lm

Sources=$(wildcard src/*.c)
IncludeDir=-I./include
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

clean:
	rm -rf obj/*.o obj/*.d $(Executables)

