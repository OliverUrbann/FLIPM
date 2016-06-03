#include <math.h>
#include "boundedness.h"

int heaviside(float t)
{
  return t >= 0 ? 1 : -1;
}

float eta(const struct Step steps[], 
          int numOfSteps, 
          const struct Robot *r, 
          float t)
{
  float om = sqrt(r->g/r->zh);
  float om2 = om * om;
  float eta = 0;
  int i;
  float _t, _l;

  for (i = 0; i < numOfSteps; i++)
  {
    _t = steps[i].time;
    _l = steps[i].length;
    if (t < _t)
      eta += _l * om2 * 0.5 * r->M / (r->k + r->b*om) * (exp(om*(t-_t)) 
				     - exp(-(r->k/r->b*t+om*_t))); 
    else
      eta +=  _l * om2 * 0.5 * r->M / (r->k + r->b*om) * (exp(-r->k/r->b*(t-_t)) 
              - exp(-(r->k/r->b*t+om*_t)))
              - _l * om2 * 0.5 * r->M / (r->k - r->b*om) * (exp(-om*(t-_t))
              - exp(-r->k/r->b*(t-_t)));
  }
  return eta;
}

float dc1d(const struct Step steps[], 
           int numOfSteps, 
           const struct Robot *r,
           float t)
{
  float om = sqrt(r->g/r->zh);
  float dc1d = 0;
  int i;
  float _t, _l;
  for (i = 0; i < numOfSteps; i++)
  {
    _t = steps[i].time;
    _l = steps[i].length;
    if (t < _t)
      dc1d += _l * om * 0.5 * exp(om * (t-_t));
    else
      dc1d += _l * om * 0.5 * exp(-om * (t-_t));
  }
  return dc1d;
}

float c1d(const struct Step steps[], 
          int numOfSteps,
          const struct Robot *r,
          float t)
{
  float om = sqrt(r->g/r->zh);
  float c1d = 0;
  int i;
  float _t, _l;
  for (i = 0; i < numOfSteps; i++)
  {
    _t = steps[i].time;
    _l = steps[i].length;
    if (t < _t)
      c1d += _l * 0.5 * exp(om * (t-_t));
    else
      c1d += _l * 0.5 * (2 - exp(-om * (t-_t)));
  }
  return c1d;
}

float controller(const struct Robot *r, const struct Step steps[],
                 int numOfSteps,        float t)
{
  float _c1d = c1d(steps, numOfSteps, r, t);
  float _eta = eta(steps, numOfSteps, r, t);
  float p = _c1d + _eta;
  return p;
}

float capturePoint(const struct Robot *r, const struct Step steps[],
                   int numOfSteps,        float t)
{
  float curAbsPos = 0;
  int i;
  float om, _c1d, _dc1d, cpStepLen;
  for (i = 0; i < numOfSteps && steps[i].time < t; i++)
    curAbsPos += steps[i].length;	
  om = sqrt(r->g/r->zh);
  _c1d = c1d(steps, numOfSteps, r, t);
  _dc1d = c1d(steps, numOfSteps, r, t);
  cpStepLen = _c1d + 1/om * _dc1d - curAbsPos;
  return cpStepLen;									 
}

void boundednessController(struct Vec *v, 
                           const struct Robot *r, 
                           int numOfSteps,
                           struct Step *steps[2],
                           float t)
{
  v->x = controller(r, steps[X], numOfSteps, t);
  v->y = controller(r, steps[Y], numOfSteps, t);
}

void boundednessCapturePoint(struct Vec *v, 
                             const struct Robot *r, 
                             int numOfSteps,
                             struct Step *steps[2],
                             float t)										
{
  v->x = capturePoint(r, steps[X], numOfSteps, t);
  v->y = capturePoint(r, steps[Y], numOfSteps, t);
}
