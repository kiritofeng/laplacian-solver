CXX = g++
CXXFLAGS = -march=native


FALSY := 0 off OFF false FALSE
TRUTHY := 1 on ON true TRUE

ifdef DEBUG
ifneq (,$(DEBUG $(VIEW),$(FALSY)))
CXXFLAGS += -O2
else
CXXFLAGS += -g -O0 -fsanitize=address -fsanitize=undefined
endif
else
CXXFLAGS += -O2
endif


%: %_driver.cc
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o bin/$@ $<

.PHONY: clean
clean:
	rm -rf bin/

DEFAULT := main
