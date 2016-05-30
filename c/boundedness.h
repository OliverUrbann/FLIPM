#pragma once
enum {X, Y, numOfDims};

struct Vec
{
	float x, y;
};

struct Robot
{
	float g, dt, zh, k, b, M;
};

struct Step
{
	float length, time;
};

void boundednessController(struct Vec *v, 
                           const struct Robot *r, 
                           int numOfSteps,
                           struct Step steps[numOfDims][numOfSteps],
                           float t);
													 
void boundednessCapturePoint(struct Vec *v, 
                             const struct Robot *r, 
                             int numOfSteps,
                             struct Step steps[numOfDims][numOfSteps],
                             float t);
