/* Described the structure that contains main simulation parameters 
Later in the code one needs to register (fill in the values) for 
variable "simulation" of type struct simulation_parameters.
To use this variable one needs to include "simulation_parameters.h" 
*/

#include "common.h"
#include <stddef.h>

#ifndef SIMULATION_PARAMETERS_H    // http://en.wikipedia.org/wiki/Include_guard 
#define SIMULATION_PARAMETERS_H

struct simulation_parameters
{
	my_float dt;			// time step
	size_t NofTS;		    // number of time steps	
};

extern struct simulation_parameters simulation;

void Register_Simulation_Parameters(struct simulation_parameters * simulation);

#endif   // SIMULATION_PARAMETERS_H
