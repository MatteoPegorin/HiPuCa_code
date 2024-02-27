# Compiler
CC = g++

# Compiler flags
CFLAGS = -O3 -std=c++17

# Linker flags
LDFLAGS =

# Object directory
OBJDIR = obj

# Source files
LIB_SOURCES = Waveform.cpp MCMC.cpp main.cpp
LIB_HEADERS = Waveform.h MCMC.h
LIB_OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(LIB_SOURCES))

# Executable
EXECUTABLE = HiPuCa_MCMC

# Default target
all: $(EXECUTABLE)

# Rule to compile with OpenMP support
parallel_openMP: CFLAGS += -fopenmp
parallel_openMP: LDFLAGS += -fopenmp
parallel_openMP: $(EXECUTABLE)

# Linking the executable
$(EXECUTABLE): $(LIB_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

# Compiling object files
$(OBJDIR)/%.o: %.cpp $(LIB_HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create object directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean target
clean:
	rm -rf $(OBJDIR) $(EXECUTABLE)
