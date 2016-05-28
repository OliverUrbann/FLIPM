#pragma once
enum {X, Y};

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
                           const struct Step steps[2][numOfSteps],
                           float t);

