CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -fopenmp
CXXOPT = -O0
INCLUDE = -I./$(SRC)
CXXCMD = $(CXX) $(CXXFLAGS) $(CXXOPT) $(INCLUDE)

headers = $(wildcard ./$(SRC)/*.hpp)
sources = $(wildcard ./$(SRC)/*.cpp)
objects_names = $(sources:./$(SRC)/%.cpp=%)
objects = $(objects_names:%=./$(SRC)/%.o)

mains = graphics/burned_probabilities_data graphics/fire_animation_data

all: $(mains)

%.o: %.cpp $(headers)
	$(CXXCMD) -c $< -o $@

$(mains): %: %.cpp $(objects) $(headers)
	$(CXXCMD) $< $(objects) -o $@

data.zip:
	wget https://cs.famaf.unc.edu.ar/~nicolasw/data.zip

data: data.zip
	unzip data.zip

clean:
	rm -f $(objects) $(mains)

.PHONY: all clean
