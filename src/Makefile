CC=gcc 
CFLAGS=-std=c99 -DUSE_DOUBLES -c

all: particle

particle: main.o electro_magnetic_field.o rk4.o simulation_parameters.o
	$(CC) main.o electro_magnetic_field.o rk4.o simulation_parameters.o -o particle.e -lm

main.o:
	$(CC) $(CFLAGS) main.c 

electro_magnetic_field.o:
	$(CC) $(CFLAGS) electro_magnetic_field.c

rk4.o:
	$(CC) $(CFLAGS) rk4.c

simulation_parameters.o:
	$(CC) $(CFLAGS) simulation_parameters.c

clean:
	rm -rf *.o particle.e
