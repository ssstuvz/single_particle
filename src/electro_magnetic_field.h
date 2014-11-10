#ifndef ELECTRO_MAGNETIC_FIELD_H
#define ELECTRO_MAGNETIC_FIELD_H
#include "common.h"

/* Define electromagnetic field parameters */

struct electro_magnetic_field
{
	my_float * Fields;
	void (*return_field)(my_float t, my_float x, my_float y, my_float z, my_float * Fields);      
	// This is a function pointer, we can put later different fields without 
    // changing the syntax http://www.cprogramming.com/tutorial/function-pointers.html
	// http://www.newty.de/fpt/fpt.html
};


extern struct electro_magnetic_field em_field;

void Register_Electro_Magnetic_Fields(struct electro_magnetic_field * em_field);
void Static_Magnetic_Field(my_float t, my_float x, my_float y, my_float z, my_float * Fields);
void Oscillating_Mirror(my_float t, my_float x, my_float y, my_float z, my_float * Fields);


#endif // ELECTRO_MAGNETIC_FIELD_H
