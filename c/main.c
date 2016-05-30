#include <stdio.h>
#include <stdlib.h>
#include "boundedness.h"

int main()
{
	const int numOfSteps = 10;       // Total number of steps used here.
	const float stepDuration = 0.5f; // The duration of a step in seconds.
	struct Step steps[numOfDims][numOfSteps];
	struct Vec pos;
	struct Robot robot;
	
	robot.zh = 0.26;
	robot.g = 9.81;
	robot.dt = 0.01;
	robot.k = 1000;
	robot.b = 200;
	robot.M = 4.5;
	 
	// Create file for data to be plotting using gnuplot
	FILE * c2dY = fopen("c2dY", "w"); // Desired y position of small mass
	
	for (int i = 0; i < numOfSteps; i++)
	{
		steps[X][i].length = 0.1;
		steps[X][i].time = stepDuration; // Note that the duration may
		                                 // change from step to step.
		steps[Y][i].length = (1 - 2 * ((i+1)%2)) * 0.1;
		steps[Y][i].time = stepDuration;
	}
	
	// This is the main loop in on a walking robot.
	for (float t = 0; t < numOfSteps * stepDuration; t += robot.dt)
	{
		// At first, on a physical robot you would need to make sure
		// that your control loop starts here with the right frequency,
		// e.g. that it starts here 10 ms after the last run if your
		// frequency is 100 Hz.
		
		boundednessController(&pos, &robot, numOfSteps, steps, t);
		fprintf(c2dY, "%lf %lf \n", t, pos.y); // Write data to plot
	}
	
	fclose(c2dY);
	
	printf("Starting gnuplot...\n");
	system("gnuplot plot");

	return 0;
}

