#include "simulation_parameters.h"

struct simulation_parameters simulation;      // this will now be visible to all files that include "simulation_parameters.h"

void Register_Simulation_Parameters(struct simulation_parameters *simulation)
{
	// this function takes the pointer to the struct simulation and puts values in
	// when you call it you need to provide a pointer, not a variable itself, i.e.
	// Register_Simulation_Parameters(&simulation)

	simulation->dt=0.01;
	simulation->dt=simulation->dt*2*PI;
	simulation->NofTS=10000;
}