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
  for (int i = 0; i < numOfSteps; i++)
  {
    float _t = steps[i].time;
    float _l = steps[i].length;
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
  for (int i = 0; i < numOfSteps; i++)
  {
    float _t = steps[i].time;
    float _l = steps[i].length;
    if (t < _t)
      dc1d += _l * om * 0.5 * exp(om * (t-_t));
    else
      dc1d += _l * om * 0.5 * exp(-om * (t-_t));
  }
  return dc1d;
}

float c1d(const std::vector<float> &T, const std::vector<float> &alpha, float t, float g, float z_h)
{
  float om = sqrt(g/z_h);
  float c1d = 0;
  for (int i = 0; i < T.size(); i++)
  {
    float curT = T[i];
    float curAlpha = alpha[i];
    if (t < curT)
      c1d += curAlpha * 0.5 * exp(om * (t-curT));
    else
      c1d += curAlpha * 0.5 * (2 - exp(-om * (t-curT)));
  }
  return c1d;
}

void boundednessController(struct Vec *v, const struct Robot *r)
{
	
}