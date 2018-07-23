CC=g++
TARGET=802.11
ROOT=./
CFLAGS=-D IEEE80211_MSSE2 -O3

SOURCE = main.cpp \
		Fourier.cpp \
		FFT_fixed.cpp \
		line_fitting.cpp \
		viterbi_decoder_generic.cpp \
		viterbi_decoder_x86.cpp \
		mapper.cpp \


## End sources definition
INCLUDE = 
## end more includes

LIBS = 

OBJ=$(join $(addsuffix ../obj/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cpp=.o))) 

## Fix dependency destination to be ../.dep relative to the src dir
DEPENDS=$(join $(addsuffix ../.dep/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cpp=.d)))

## Default rule executed
all: $(TARGET)
		@true

## Clean Rule
clean:
		@-rm -f $(TARGET) $(OBJ) $(DEPENDS)

## Rule for making the actual target
$(TARGET): $(OBJ)
		@echo "============="
		@echo "Linking the target $@"
		@echo "============="
		@$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
		@strip $@
		@echo -- Link finished --

## Generic compilation rule
%.o : %.cpp
		@mkdir -p $(dir $@)
		@echo "Compiling $<"
		$(CC) $(CFLAGS) -cccccc $< -o $@ $(INCLUDE)


## Rules for object files from cpp files
## Object file for each file is put in obj directory
## one level up from the actual source directory.
../obj/%.o : %.cpp
		@mkdir -p $(dir $@)
		@echo "Compiling $<"
		@$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)

# Rule for "other directory"  You will need one per "other" dir
$(HAL3288)/../obj/%.o : %.cpp
		@mkdir -p $(dir $@)
		@echo "Compiling $<"
		@$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)

# Rule for "other directory"  You will need one per "other" dir
$(FEC)/../obj/%.o : %.cpp
		@mkdir -p $(dir $@)
		@echo "Compiling $<"
		@$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)

## Make dependancy rules
../.dep/%.d: %.cpp
		@mkdir -p $(dir $@)
		@echo "============="
		@echo Building dependencies file for $*.o
		@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^../obj/$*.o^" > $@'

## Dependency rule for "other" directory
$(FEC)/../.dep/%.d: %.cpp
		@mkdir -p $(dir $@)
		@echo "============="
		@echo Building dependencies file for $*.o
		@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OTHERDIR)/../obj/$*.o^" > $@'

## Dependency rule for "other" directory
$(HAL3288)/../.dep/%.d: %.cpp
		@mkdir -p $(dir $@)
		@echo "============="
		@echo Building dependencies file for $*.o
		@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OTHERDIR)/../obj/$*.o^" > $@'

## Include the dependency files
-include $(DEPENDS)