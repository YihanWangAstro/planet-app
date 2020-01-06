all: explore mass ueq_long evn

VPATH =${MY_REPO_PATH}/SpaceX:${MY_REPO_PATH}/SpaceHub/src/
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall  
OBJS = dynamicSystem.o InitFileCreator.o \
	 chain.o chainSystem.o tools.o\
	 IO.o symplectic.o regularizator.o BSIterator.o\
	particleSystem.o integrator.o seqIterator.o Wrapper.o
	
INCL = dynamicSystem.h macros.h chain.h\
	regularizator.h particleSystem.h interaction.h chainSystem.h \
	integrator.h BSIterator.h errhand.h InitFileCreator.h\
	seqIterator.h odeIterator.h vector3d.h tools.h

explore : $(OBJS) explore_v_a.o
	$(CXX) -std=c++17 -O3 -o explore.ex $(OBJS) explore_v_a.o -pthread

mass : $(OBJS) explore_mass.o
	$(CXX) -std=c++17 -O3 -o explore_mass.ex $(OBJS) explore_mass.o -pthread

ueq_long: $(OBJS) ueq_longer.o
	$(CXX) -std=c++17 -O3 -o ueq_longer.ex $(OBJS) ueq_longer.o -pthread

evn: $(OBJS) evn-distro.o
	$(CXX) -std=c++17 -O3 -o evn-distro.ex $(OBJS) evn-distro.o -pthread

$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) *.o *.ex
