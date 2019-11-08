EXE = flocking

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CXX = h5c++
CPPFLAGS += -Iinclude -I/home/simongh/anaconda3/include
CXXFLAGS += -Wall
LDFLAGS += -Llib -L/home/simongh/anaconda3/lib
LDLIBS += -lm -lhdf5_cpp

# .PHONY: all clean

all: $(EXE)

debug: CPPFLAGS += -DDEBUG
debug: CXXFLAGS += -g
debug: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ)

print-% : ; @echo $* = $($*)
