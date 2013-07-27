BASE   = $(shell pwd)/..

OUTPUT = libsoc.so

OPTIM = fast
#OPTIM  = fast_debug
#OPTIM  = debug

CXXFLAGS = -I../src

### 1st option: compile everything in a single obj file...
#CXXFLAGS += -DMT_IMPLEMENTATION #causes libsoc.cpp to include all *.cpp sources
#OBJS   = libsoc.o


### 2nd option: compile all objs separately
CXXFLAGS += -DMT_IMPLEMENT_TEMPLATES #causes all *.cpp to implement the used templates
OBJS  = $(shell find ../src/MT \( -name '*.cpp' -and -not -name '*_t.cpp' -and -not -name '*_instantiate.cpp' \) )
OBJS := libsoc.o $(OBJS:%.cpp=%.o)


include $(BASE)/make-generic
