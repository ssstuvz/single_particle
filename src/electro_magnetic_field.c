#include "electro_magnetic_field.h"
#include <math.h>
#include <stdio.h>
#include "common.h"


struct electro_magnetic_field em_field;

//// Here one can describe different fields and use function pointer to return

void Gaussian_Plane_Wave(my_float t, my_float x, my_float y, my_float z, my_float * Fields)
{
	// Fields[0] - Ex
	// Fields[1] - Ey
	// Fields[2] - Ez
	// Fields[3] - Bx
	// Fields[4] - By
	// Fields[5] - Bz

	// gaussian e-m pulse parameters
	my_float a0=0.1;          // pulse amplitude
	my_float Duration=10;		// pulse Duration in cycles
	my_float Start=50;        // temporal offset (simulation starts at t=0)
	int Polarization=0;
	my_float phase=0;

	// convert to omega*t normalization
	Duration=Duration*2*PI;
	Start=Start*2*PI;


	switch(Polarization)
	{
		case 0: // ex, by polarization
			Fields[0]=a0*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase);
			Fields[1]=0.0;
			Fields[2]=0.0;
			Fields[3]=0.0;
			Fields[4]=Fields[0];
			Fields[5]=0.0;
			break;
		case 1: // ey, bx
			break;
		case 2:
			Fields[0]=a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*sin(t-z-Start+phase);
			Fields[1]=a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*sin(t-z-Start+phase+PI*0.5);
			Fields[2]=0.0;
			Fields[3]=-Fields[1];
			Fields[4]=Fields[0];
			Fields[5]=0.0;
			break;
	}
	//Fields[0]=1.0;
	//Fields[1]=3.0;

}

void Register_Electro_Magnetic_Fields(struct electro_magnetic_field * em_field)
{
	em_field->return_field=&Gaussian_Plane_Wave;   // function pointer to needed field
}
