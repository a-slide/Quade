############################# VARIABLES ################################

PYVERSION = 2.7

CY = cython
	# Cython compilater

CC = gcc
	# C Compiler

CFLAGS = -Os -I /usr/include/python$(PYVERSION)
	# Compilation options

LFLAGS = -lpython$(PYVERSION) -lpthread -lm -lutil -ldl
	# link edition options

CYFLAGS = --embed
	# cython compiler options

#################### INSTRUCTIONS DE COMPILATION #######################
# $@ =  Target # $^ = list of dependencies # $< First dependency #

all: Quade

Quade: Quade.o
	#Link editing
	$(CC) -o $@ $^ $(LFLAGS)
	
Quade.o: Quade.c
	#Compilation of source object
	$(CC) -c $^ -o $@ $(CFLAGS)

Quade.c: Quade.py
	#Cythonization
	$(CY) $^ -o $@ $(CYFLAGS)


##################### INSTRUCTIONS DE NETTOYAGE ########################

.PHONY: clean mrproper

clean:
	#Clean intermediate files
	rm -f Quade.c Quade.o

mrproper:
	#Clean everything
	rm -f Quade Quade.c Quade.o
