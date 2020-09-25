.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -Wall -D NDEBUG -m64
	
SOURCE		= main.cpp GetData.cpp Tools.cpp Prediction.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -lpthread -o GlycoPredictWebServer $(LIB)

%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $< 

clean:
			rm -f *.o *~
