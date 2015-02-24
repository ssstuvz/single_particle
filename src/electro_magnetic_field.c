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
	my_float a0=3.0;          // pulse amplitude
	my_float Duration=7*2.668223;		// pulse Duration in cycles
	my_float Start=40;        // temporal offset (simulation starts at t=0)
	int Polarization=0;
	my_float phase=1.57;

	// convert to omega*t normalization
	//Duration=Duration*2*PI;
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

void gaussianPulseLeft(my_float t, my_float x, my_float y, my_float z, my_float * Fields)
{
	my_float a0 = 3.0;
	my_float duration = 7*2.668223;
	my_float x0 = 100;
	my_float phase = 1.57;

    Fields[0]=a0*exp(-(t+z-x0)*(t+z-x0)/2.0/duration/duration)*sin(t+z+phase-x0); // Doesn't work as expected
	//Fields[0]=a0*exp(-(-t+z+x0)*(-t+z+x0)/2.0/duration/duration)*sin(-t+z+phase+x0); // Works
	Fields[1]=0.0;
	Fields[2]=0.0;
	Fields[3]=0.0;
    Fields[4]=-Fields[0];
    Fields[5]=0.0;
}

void Polarization_Gating(my_float t, my_float x, my_float y, my_float z, my_float * Fields)
{
	// Fields[0] - Ex
	// Fields[1] - Ey
	// Fields[2] - Ez
	// Fields[3] - Bx
	// Fields[4] - By
	// Fields[5] - Bz
	// Fields[5] - n_e

	// gaussian e-m pulse parameters
	my_float a0=5.0;          // pulse amplitude
	my_float Duration=7.0;		// pulse Duration in cycles
	my_float Start=-40.0;        // temporal offset (simulation starts at t=0)
	my_float phase=0;
	my_float Delay=2.0; // delay in periods

	// mirror parameters
	my_float Density = 200.0; // electron density
	my_float omega_p = sqrt(Density);
	my_float alpha = atan(omega_p);
	phase=alpha;

	// convert to omega*t normalization
	Duration=Duration*2.885;
	Start=Start*2*PI;
	Delay=Delay*2*PI;


/*	Fields[0]= 2*(PulseFunction(z, t, phase, Start, Duration, a0)+0*PulseFunction(z, t, phase, Start+Delay, Duration, a0))/sqrt(1+Density);
	Fields[1]= 2*(PulseFunction(z, t, phase, Start+PI*0.5, Duration, a0)+0*PulseFunction(z, t, phase, Start-PI*0.5+Delay, Duration, a0))/sqrt(1+Density);
	Fields[2]=Density*z;
	Fields[3]=-2*omega_p*(PulseFunction(z, t, phase+PI/2.0, Start, Duration, a0)+0*PulseFunction(z, t, phase+PI/2.0, Start+Delay, Duration, a0))/sqrt(1+Density);
	Fields[4]=2*omega_p*(PulseFunction(z, t, phase+PI/2.0, Start+PI*0.5, Duration, a0)+0*PulseFunction(z, t, phase+PI/2.0, Start-PI*0.5+Delay, Duration, a0))/sqrt(1+Density);
	Fields[5]=0.0;
	*/
	Fields[0]= 2*(PulseFunction(z, t, phase, Start, Duration, a0)+PulseFunction(z, t, phase, Start+Delay, Duration, a0))/sqrt(1+Density);
	Fields[1]= 2*(PulseFunction(z, t, phase+PI*0.5, Start, Duration, a0)+PulseFunction(z, t, phase-PI*0.5, Start+Delay, Duration, a0))/sqrt(1+Density);
	Fields[2]=Density*z;
	Fields[3]= -omega_p*2*(PulseFunction(z, t, phase, Start, Duration, a0)+PulseFunction(z, t, phase-PI, Start+Delay, Duration, a0))/sqrt(1+Density);;
	Fields[4]= omega_p*2*(PulseFunction(z, t, phase-PI*0.5, Start, Duration, a0)+PulseFunction(z, t, phase-PI*0.5, Start+Delay, Duration, a0))/sqrt(1+Density);;
	//Fields[3]=-2*omega_p*(PulseFunction(z, t, phase, Start, Duration, a0)+PulseFunction(z, t, phase, Start+Delay, Duration, a0))/sqrt(1+Density);
	//Fields[4]=2*omega_p*(PulseFunction(z, t, phase+PI*0.5, Start+PI*0.5, Duration, a0)+PulseFunction(z, t, phase+PI*0.5, Start-PI*0.5+Delay, Duration, a0))/sqrt(1+Density);
	Fields[5]=0.0;


/*	Fields[0]=2*a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*cos(t-z-Start+phase+alpha)/sqrt(1+Density);
				Fields[1]=2*a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*cos(t-z-Start+phase+alpha+PI*0.5)/sqrt(1+Density);
				Fields[2]=Density*z;
				Fields[3]=-2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha+PI*0.5)/sqrt(1+Density);
				Fields[4]=2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha)/sqrt(1+Density);
				Fields[5]=0.0;

*/
	//Fields[0]=1.0;
	//Fields[1]=3.0;

}

double PulseFunction(double x, double t, double phase, double start, double duration, double a0) {
	return a0*cos(t-x+phase+start)*exp(-(start-x+t)*(start-x+t)/(2*duration*duration));
}

void Oscillating_Mirror(my_float t, my_float x, my_float y, my_float z, my_float * Fields,  my_float gamma)
{
	// Fields[0] - Ex
	// Fields[1] - Ey
	// Fields[2] - Ez
	// Fields[3] - Bx
	// Fields[4] - By
	// Fields[5] - Bz
	// Fields[5] - n_e

	// gaussian e-m pulse parameters
	my_float a0=10.0;          // pulse amplitude
	my_float Duration=4;		// pulse Duration in cycles
	my_float Start=10;        // temporal offset (simulation starts at t=0)
	int Polarization=0;
	my_float phase=0;

	// mirror parameters
	my_float Density = 400.0; // electron density
	my_float omega_p = sqrt(Density)/sqrt(gamma);
	my_float alpha = atan(omega_p);

	// convert to omega*t normalization
	Duration=Duration*2.885;
	Start=Start*2*PI;


	switch(Polarization)
	{
		case 0: // ex, by polarization
			Fields[0]=2*a0*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*cos(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[1]=0.0;
			Fields[2]=Density*z;
			Fields[3]=0.0;
			Fields[4]=2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[5]=0.0;
		case 1: // ey, bx
			break;
		case 2:
			Fields[0]=2*a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*cos(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[1]=2*a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*cos(t-z-Start+phase+alpha+PI*0.5)/sqrt(1+Density);
			Fields[2]=Density*z;
			Fields[3]=-2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha+PI*0.5)/sqrt(1+Density);
			Fields[4]=2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[5]=0.0;
			break;
	}
	//Fields[0]=1.0;
	//Fields[1]=3.0;

}

void Intensity_Gating(my_float t, my_float x, my_float y, my_float z, my_float * Fields)
{
	// Fields[0] - Ex
	// Fields[1] - Ey
	// Fields[2] - Ez
	// Fields[3] - Bx
	// Fields[4] - By
	// Fields[5] - Bz
	// Fields[5] - n_e

	// gaussian e-m pulse parameters
	my_float a0=5.0;          // pulse amplitude
	my_float Duration=4;		// pulse Duration in cycles
	my_float Start=10;        // temporal offset (simulation starts at t=0)
	int Polarization=0;
	

	// mirror parameters
	my_float Density = 100.0; // electron density
	my_float omega_p = sqrt(Density);
	my_float alpha = atan(omega_p);
	
	my_float phase=-1.57;

	// convert to omega*t normalization
	Duration=Duration*2.885;
	Start=Start*2*PI;


	switch(Polarization)
	{
		case 0: // ex, by polarization
			Fields[0]=2*a0*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*cos(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[1]=0.0;
			Fields[2]=Density*z;
			Fields[3]=0.0;
			Fields[4]=2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[5]=0.0;
		case 1: // ey, bx
			break;
		case 2:
			Fields[0]=2*a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*cos(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[1]=2*a0*exp(-pow(t-z-Start,2)/2.0/Duration/Duration)*cos(t-z-Start+phase+alpha+PI*0.5)/sqrt(1+Density);
			Fields[2]=Density*z;
			Fields[3]=-2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha+PI*0.5)/sqrt(1+Density);
			Fields[4]=2*a0*omega_p*exp(-pow(t-z-Start,2)/2.0/pow(Duration,2))*sin(t-z-Start+phase+alpha)/sqrt(1+Density);
			Fields[5]=0.0;
			break;
	}
	//Fields[0]=1.0;
	//Fields[1]=3.0;

}

void Static_Magnetic_Field(my_float t, my_float x, my_float y, my_float z, my_float * Fields)
{
	my_float B0 = 0.5; // Field Amplitude
	my_float boundary_min_x=-100, boundary_min_y=-100;
	my_float boundary_max_x=100, boundary_max_y=100;

	if (x > boundary_min_x && x < boundary_max_x && y > boundary_min_y && y < boundary_max_y)
	{
		Fields[5] = B0;
	}
	else
	{
		Fields[5] = 0;
	}


}

void Register_Electro_Magnetic_Fields(struct electro_magnetic_field * em_field)
{
	//em_field->return_field=&Gaussian_Plane_Wave;   // function pointer to needed field
	em_field->return_field=&Oscillating_Mirror;
}
