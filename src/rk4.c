#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "rk4.h"
#include "electro_magnetic_field.h"

void derivatives_electron(my_float t, my_float y[], my_float dy_dt[])
{
	// y[0]=x		Fields[0]=Ex
	// y[1]=y 		Fields[1]=Ey
	// y[2]=z		Fields[2]=Ez
	// y[3]=px		Fields[3]=Bx
	// y[4]=py		Fields[4]=By
	// y[5]=pz		Fields[5]=Bz

	// dx/dt=px/gamma

	my_float Fields[6];

	my_float gamma=sqrt(1+pow(y[3],2) + pow(y[4],2) + pow(y[5], 2)); // relativistic gamma factor
	dy_dt[0]=y[3]/gamma;
	dy_dt[1]=y[4]/gamma;
	dy_dt[2]=y[5]/gamma;

	em_field.return_field(t, y[0], y[1], y[2], Fields);

	dy_dt[3]=-Fields[0]-(y[4]*Fields[5]-y[5]*Fields[4]);
	dy_dt[4]=-Fields[1]-(y[5]*Fields[3]-y[3]*Fields[5]);
	dy_dt[5]=-Fields[2]-(y[3]*Fields[4]-y[4]*Fields[3]);

}


void rk4_step(my_float y_in[], int number_of_eqns, my_float t,
 my_float dt, my_float y_out[], void (*derivatives)(my_float, my_float[], my_float[]))
 {
 	my_float t_temporary=t;
 	my_float y_temporary[6];
 	my_float dy_dt[6];
 	my_float OneSixth=1.0/6.0;  // saves some operations
 	my_float OneThird=1.0/3.0;
 	int n=0; // iterator for for loops

////////////

		derivatives(t, y_in, dy_dt); // gives us k1

		for(n=0; n<number_of_eqns; n++)
		{
			y_temporary[n]=y_in[n]+0.5*dy_dt[n]*dt;
			y_out[n]=y_in[n]+dt*dy_dt[n]*OneSixth;
		}
		
		derivatives(t+0.5*dt, y_temporary, dy_dt); // gives us k2

		for(n=0; n<number_of_eqns; n++)
		{
			y_temporary[n]=y_in[n]+0.5*dy_dt[n]*dt;
			y_out[n]+=dt*dy_dt[n]*OneThird;	
		}

		derivatives(t+0.5*dt, y_temporary, dy_dt);  // k3
		for(n=0; n<number_of_eqns; n++)
		{
			y_temporary[n]=y_in[n]+dy_dt[n]*dt;
			y_out[n]+=dt*dy_dt[n]*OneThird;	
		}

		derivatives(t+dt, y_temporary, dy_dt); // k4
		for(n=0; n<number_of_eqns; n++)
		{
			y_out[n]+=dt*dy_dt[n]*OneSixth;
		}
		
////////////
 }