CC = g++
LIBS = 
CPPFLAGS = 

ifeq ($(OPTIMIZE),no)
  CPPFLAGS += 
else
  CPPFLAGS += -O2
endif


lal: liquid_affinity.cpp
	$(CC) liquid_affinity.cpp $(LIBS) -o lal $(CPPFLAGS)
