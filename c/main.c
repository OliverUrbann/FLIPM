#include <stdio.h>
#include <stdlib.h>
#include "boundedness.h"

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
	struct Step *steps[numOfDims];
   Step stepX[numOfSteps];
	struct Vec pos;
	struct Robot robot;
   FILE *data;
   int i;
   float stepLenY, footPosY;
   float t;

   steps[X] = stepX;
	
	robot.zh = 0.26;
	robot.g = 9.81;
	robot.dt = 0.01;
	robot.k = 1000;
	robot.b = 200;
	robot.M = 4.5;
	 
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
		
	// This is the main loop in on a walking robot.
	for (t = 0; t < numOfSteps * stepDuration; t += robot.dt)
	{
		// At first, on a physical robot you would need to make sure
		// that your control loop starts here with the right frequency,
		// e.g. that it starts here 10 ms after the last run if your
		// frequency is 100 Hz.
		
		boundednessController(&pos, &robot, numOfSteps, steps, t);
		
		footPosY = stepSum(numOfSteps, steps[Y], t);
		fprintf(data, "%lf %lf %lf\n", t, pos.y, footPosY); // Write data to plot
	}
	
	fclose(data);
	
	printf("Starting gnuplot...\n");
	system("gnuplot plot");

	return 0;
}

