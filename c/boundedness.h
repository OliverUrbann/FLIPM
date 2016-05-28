#pragma once

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

void boundednessController(struct Vec *v, const struct Robot *r);