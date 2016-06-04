#include <stdio.h>
#include <stdlib.h>
#include "boundedness.h"
#include <time.h>

float stepSum(int numOfSteps, struct Step *steps, float t)
{
  float sum = 0;
  int i;
	for (i = 0; i < numOfSteps && steps[i].time < t; i++)
		sum += steps[i].length;
	return sum;
}

int main()
{
	const int numOfSteps = 10;       // Total number of steps used here.
	const float stepDuration = 0.5f; // The duration of a step in seconds.
	struct Step stepX[numOfSteps], stepY[numOfSteps];
	struct Step *steps[numOfDims] = {stepX, stepY};
	struct Vec pos;
	struct Robot robot;
  FILE *data;
  int i;
  float stepLenY, footPosY, footPosX;
  float t;
	float time_spent = 0;
	clock_t begin, end;

	// Height of center of mass over ground
	robot.zh = 0.26; // [m]
	
	// The gravity
	robot.g = 9.81; // [m/s^2]
	
	// Duration of one frame. You may also try here higher frequencies.
	robot.dt = 0.01; // [s]
	
	// Spring constant
	robot.k = 1000;
	
	// Damper constant
	robot.b = 200;
	
	// Weight of mass 1
	robot.M = 4.5; // [kg]
	 
	// Create file for data to be plotting using gnuplot
	data = fopen("data", "w"); // Desired y position of small mass
	fprintf(data, "Time c2dY pRefY\n");
	for (i = 0; i < numOfSteps; i++)
	{
		steps[X][i].length = 0.1;
		steps[X][i].time = stepDuration * (i + 1); // Note that the duration may
                                               // change from step to step.
    stepLenY = 0.1f;
    if (i == 0) stepLenY = 0.05;                                           
		steps[Y][i].length = (1 - 2 * (i%2)) * stepLenY;
		steps[Y][i].time = stepDuration * (i + 1);
	}
		
	// This is the main loop in on a walking robot. Depending on the
	// architecture you would not have a loop, just a function with the
	// inner of this loop that is executed every dt seconds.
	for (t = 0; t < numOfSteps * stepDuration; t += robot.dt)
	{	
		begin = clock();
		
		// The only call you need in fact.
		boundednessController(&pos, &robot, numOfSteps, steps, t);
		
		end = clock();
		time_spent += (float)(end - begin) / CLOCKS_PER_SEC;
			
		footPosY = stepSum(numOfSteps, steps[Y], t);
		footPosX = stepSum(numOfSteps, steps[X], t);
		fprintf(data, "%lf %lf %lf %lf %lf\n",
		        t, pos.y, footPosY, pos.x, footPosX);
	}
	printf("Time for calculation: %lf s\n", time_spent);
	
	fclose(data);
	
	printf("Starting gnuplot...\n");
	system("gnuplot plot");

	return 0;
}

