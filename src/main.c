#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "simulation_parameters.h"
#include "electro_magnetic_field.h"
#include "rk4.h"

int main()
{
	printf("Simple one charged particle in EM-field RK4 code\n");

	Register_Simulation_Parameters(&simulation);
	printf("Registered simulation parameters:\n");
	Register_Electro_Magnetic_Fields(&em_field);
	printf("Registered electro-magnetic parameters:\n");
	
	printf("dt = %e\n", simulation.dt);
	printf("NofTS = %d\n", simulation.NofTS);

	// check which type we are using for floats: double of float?
	// at compilation time to use double precision use, for example, gcc -DUSE_DOUBLES *.c
	if(sizeof(my_float)==sizeof(float)) printf("Using float instead of double\n");
	else printf("Using double instead of float\n");
	
	my_float time=0.0;
	my_float Fields[6];
	
	my_float * x;
	my_float * y;
	my_float * z;
	my_float * px;
	my_float * py;
	my_float * pz;

	x = (my_float *)malloc(simulation.NofTS*sizeof(my_float));
	y = (my_float *)malloc(simulation.NofTS*sizeof(my_float));
	z = (my_float *)malloc(simulation.NofTS*sizeof(my_float));
	px = (my_float *)malloc(simulation.NofTS*sizeof(my_float));
	py = (my_float *)malloc(simulation.NofTS*sizeof(my_float));
	pz = (my_float *)malloc(simulation.NofTS*sizeof(my_float));

// initial conditions
	x[0]=0.0;
	y[0]=0.0;
	z[0]=0.0;
	px[0]=0.0;
	py[0]=0.0;
	pz[0]=0.0;

	my_float y_in[6];   /// needed for Runge-Kutta step
	my_float y_out[6];  /// needed for Runge-Kutta step

	for(int n=1; n<simulation.NofTS; n++)
	{
		time=n*simulation.dt;
		y_in[0]=x[n-1];
		y_in[1]=y[n-1];
		y_in[2]=z[n-1];
		y_in[3]=px[n-1];
		y_in[4]=py[n-1];
		y_in[5]=pz[n-1];

		rk4_step(y_in, 6, time, simulation.dt, y_out, derivatives_electron);

		x[n]=y_out[0];
		y[n]=y_out[1];
		z[n]=y_out[2];
		px[n]=y_out[3];
		py[n]=y_out[4];
		pz[n]=y_out[5];

	}
	

// Dump data to disk as binary
// Be careful to read correct amount of bytes depending on my_float
	FILE * SaveFile=fopen("px.dat", "wb");
	fwrite(px, simulation.NofTS, sizeof(my_float), SaveFile);
	fclose(SaveFile);

	SaveFile=fopen("py.dat", "wb");
	fwrite(py, simulation.NofTS, sizeof(my_float), SaveFile);
	fclose(SaveFile);

	SaveFile=fopen("z.dat", "wb");
	fwrite(z, simulation.NofTS, sizeof(my_float), SaveFile);
	fclose(SaveFile);
	
	SaveFile=fopen("pz.dat", "wb");
	fwrite(pz, simulation.NofTS, sizeof(my_float), SaveFile);
	fclose(SaveFile);
////////////

	free(x);
	free(y);
	free(z);
	free(px);
	free(py);
	free(pz);


	printf("End of simulation...\n");
	return 0;
}
