# Author:       Oliver Urbann
# E-Mail:       oliver.urbann@tu-dortmund.de
# Requirements: Python 3.4+, scipy, numpy, matplotlib, tkinter
# Usage:        python model.py
#
# This script is the implementation of the FLIP model as proposed by
# Urbann et al. [1]. Details about the gains used in the GUI can be
# found in [1] and [2].
#
# Two functions can be executed:
#
# Menu Commands->Simulate for an example of the model. At 0.5s a
# force is applied such that both masses begin to oscillate. Nothing 
# else happens here.
#
# Menu Commands->Controller runs a walk with a preview controller.
# Gains for controller are calculate before, using the given values.
#
# Values like masses, gravity, height, etc. can be selected manually.
# For a quick simulation or walk default values can be loaded by
# clicking the corresponding menu item before "Simulate" or 
# "Controller" is clicked. 
#
# Additionally, the gains and the computed controller gains can be
# saved for an application of the controller implemented on a 
# physical robot.
#
# "python model.py -c" and "python model.py -r" creates figures
# shown in [1]. Before this can be done run 
# "mkdir paper && mkdir paper/build". 
#
# [1] Flexible Linear Inverted Pendulum Model for Cost-Effective 
#     Biped Robots
#     Oliver Urbann, Ingmar Schwarz, Matthias Hofmann
#     Humanoid Robots (Humanoids), 2015 15th IEEE-RAS International
#     Conference on, 2015, pp. 128–131 
# [2] Observer-based biped walking control, a sensor fusion approach
#     Oliver Urbann, Stefan Tasse
#     Autonomous Robots 35.1 (2013) pp. 37–49. Springer US, 2013

import numpy as np
import scipy as sp
import scipy.linalg
import matplotlib as mpl
import matplotlib.pyplot as plt
import tkinter
import tkinter.filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import sys
import time

__author__ = "Oliver Urbann"
__copyright__ = "Copyright 2015, Oliver Urbann"
__credits__ = ["Oliver Urbann"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Oliver Urbann"
__email__ = "oliver.urbann@tu-dortmund.de"
__status__ = "Production"

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# To set breakpoint: import pdb; pdb.set_trace()

# Calculate acc for mass 2
#       L
#     <====
#    O-----o
#   x1     x2


# x_{t+1} = A * x_t + B * u(t)

# Spring
# F = D * L, F = m * a
# D * L = m * a
# a = D * L / m
# a = D * (x1 - x2) / m
# a = D / m * x1 - D / m * x2

# Damping has the same derivation

simDefault = (0.1, 1, 9.81, 0.3, 0.01, 1, 0.1, 1, 1, 10**-10, 100, 5)

# Schwierig zu regeln für Preview, Qx 150, Qe 100
# Höhere Frequenz hilft hier, aber dauert länger, zb 0.0025 N=400
# b = 2 einfacher
lanariControllerDefault = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.001, 1000, 10**-10, 100, 3)
# tau Problem bei höheren Dämpfungen
# Qe bei Lanari von 0.0004 ist bei E=200 besser, höhere Frequenzen ebenfalls besser
controllerDefault = (0.1, 4.5, 9.81, 0.26, 0.01, 5000, 200, 1, 10, 10**-10, 100, 5)

# Für Tests mit CP
#controllerDefault = (0.1, 4.5, 9.81, 0.26, 0.005, 10000, 200, 10, 10, 10**-10, 200, 5)


def main():
  noGui = False
  for arg in sys.argv:
    if arg == '-c': # Plot figure as shown in [1]
      fig1 = plt.figure(figsize=(6, 2), dpi=150)
      ax1 = fig1.add_subplot(111)
      control(controllerDefault[4], int(sp.floor(controllerDefault[10])), controllerDefault[11], ax1, *flipm_gains(*controllerDefault[:11]))
      plt.savefig("paper/build/FLIPMex" + ".pdf", format='pdf', bbox_inches='tight')
      noGui = True
    if arg == '-r': # Plot another figure shown in [1]
      fig1 = plt.figure(figsize=(6, 2), dpi=150)
      ax1 = fig1.add_subplot(111)
      g = flipm_gains(*simDefault[:11])
      sim(simDefault[4],simDefault[11], g[0], g[1], g[2], ax1)
      plt.savefig("paper/build/FLIPMmodel" + ".pdf", format='pdf', bbox_inches='tight')
      noGui = True
    if arg == '-b':
      plotBounded()
      noGui = True
    if arg.startswith('-l='): # Print gain tables in latex
      m.dump(open(arg[3:] + "A.tex",'w'))
      noGui = True
  if noGui:
    exit()
  root = tkinter.Tk()
  app = FLIPMApp(root)
  app.mainloop()
 
def plotBounded():
  # Figure for working stopping step
  print("Plotting iterative.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 2, 0.001, 1000, 10**-10, 100, 3)
  A, b, c, *r = flipm_gains(*params[:11])
  bcontrol(A, b, c, params[11], params[0], params[1], params[2],
           params[3], params[4], params[5], params[6],
           (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
           (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
           params[7],
           ax1)
  plt.savefig("iterative" + ".pdf", format='pdf', bbox_inches='tight')
  
  print("Plotting b=20,tau=001.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 20, 0.001, 1000, 10**-10, 100, 2.5)
  A, b, c, *r = flipm_gains(*params[:11])
  bcontrol(A, b, c, params[11], params[0], params[1], params[2],
           params[3], params[4], params[5], params[6],
           (0.5, 1, 1.5, 20, 2.5, 3, 3.5, 4, 4.5),
           (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
           params[7],
           ax1)
  plt.savefig("b=20,tau=001" + ".pdf", format='pdf', bbox_inches='tight')
  
  print("Plotting b=200,tau=001.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.001, 1000, 10**-10, 100, 2.5)
  A, b, c, *r = flipm_gains(*params[:11])
  bcontrol(A, b, c, params[11], params[0], params[1], params[2],
           params[3], params[4], params[5], params[6],
           (0.5, 1, 1.5, 20, 2.5, 3, 3.5, 4, 4.5),
           (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
           params[7],
           ax1)
  plt.savefig("b=200,tau=001" + ".pdf", format='pdf', bbox_inches='tight')
  
  print("Plotting b=200,tau=0005.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.0005, 1000, 10**-10, 100, 2.5)
  A, b, c, *r = flipm_gains(*params[:11])
  bcontrol(A, b, c, params[11], params[0], params[1], params[2],
           params[3], params[4], params[5], params[6],
           (0.5, 1, 1.5, 20, 2.5, 3, 3.5, 4, 4.5),
           (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
           params[7],
           ax1)
  plt.savefig("b=200,tau=0005" + ".pdf", format='pdf', bbox_inches='tight')
  
  print("Plotting b=200,tau=00025.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.00025, 1000, 10**-10, 100, 2.5)
  A, b, c, *r = flipm_gains(*params[:11])
  bcontrol(A, b, c, params[11], params[0], params[1], params[2],
           params[3], params[4], params[5], params[6],
           (0.5, 1, 1.5, 20, 2.5, 3, 3.5, 4, 4.5),
           (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
           params[7],
           ax1)
  plt.savefig("b=200,tau=00025" + ".pdf", format='pdf', bbox_inches='tight')
  
  return
  
  print("preview.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.1, 0.001, 10**-10, 2000, 3)
  # Preview 2000 and dt = 0.0005
  A, b, c, *r = flipm_gains(*params[:11])
  control(0, pref_y, params[2], params[3], params[4], params[10],
          params[11], ax1, *flipm_gains(*params[:11]))
  plt.savefig("preview" + ".pdf", format='pdf', bbox_inches='tight')
  
  
  
def getValues(d):
  return (d["Small Mass"].get(),
          d["Large Mass"].get(),
          d["Gravity"].get(),
          d["Height"].get(),
          d["Frame Length"].get(),
          d["Spring Constant"].get(),
          d["Damper Constant"].get(),
          d["Qe"].get(),
          d["Qx"].get(),
          d["R"].get(),
          int(sp.floor(d["N"].get())),
          d["End"].get())
          
def lipm_gains(m, g, z_h, dt, Qe, Qx, R, N):
  gz = g/z_h*dt
  A = np.array(
      [[      1,     dt,       0],
       [     gz,      1,     -gz],
       [      0,      0,       1]])
  b = np.matrix(np.array([0, 0, dt])).transpose()
  c = np.matrix(np.array([0, 0, 1]))
  c_x1 = np.matrix(np.array([1, 0, 0]))
  c_x2 = np.matrix(np.array([0, 0, 0]))
  return A, b, c, c_x1, c_x2, gains(A, b, c, Qe, Qx, R, N)
  
def flipm_zmp(m, M, g, z_h, dt, D, E, Qe, Qx, R, N):
  
    DM = D/M
    EM = E/M
    Dm = D/m
    Em = E/m
    dt2 = dt**2/2
    gz = g/z_h*dt
    
    A = np.array(
  # From vector element        ------------->             to vector element
  #                                                                     ||
  #        x1       dx1        p       dp       x2       dx2   ddx2     \/
  [[        1,       dt,       0,       0,       0,        0,     0 ], # x1
   [       gz,        1,     -gz,       0,       0,        0,     0 ], # dx1
   [        0,        0,       1,      dt, dt**2/2,        0,     0 ], # p
   [        0,        0,  -DM*dt, 1-EM*dt,   DM*dt,    EM*dt,     0 ], # dp
   [        0,        0,       0,       0,       1,       dt,   dt2 ], # x2
   [        0,        0,   Dm*dt,   Em*dt,  -Dm*dt,  1-Em*dt,    dt ], # dx2
   [        0,        0,       0,       0,       0,        0,     1 ]])# ddx2
 

    b = np.matrix(np.array([0, 0, 0, 0, dt**3 / 6, dt**2 / 2, dt])).transpose()
    c = np.matrix(np.array([0, 0, 1, 0, 0, 0, 0]))
    c_x1 = np.matrix(np.array([1, 0, 0, 0, 0, 0, 0]))
    c_x2 = np.matrix(np.array([0, 0, 0, 0, 1, 0, 0]))
    return A, b, c, c_x1, c_x2, (*gains(A, b, c, Qe, Qx, R, N))

def flipm_gains4D(m, M, g, z_h, dt, D, E, Qe, Qx, R, N):
  # Calculation of controller flipm_gains as explained in [2].
  A = np.array(
      [[      1,      dt,        0,        0],
       [-D/M*dt,1-E/M*dt,   D/M*dt,   E/M*dt],
       [      0,       0,        1,       dt],
       [ D/m*dt,  E/m*dt,  -D/m*dt, 1-E/m*dt]])
  b = np.matrix(np.array([0, 0, dt**2 / 2, dt])).transpose()
  c = np.matrix(np.array([1+D/M*z_h/g,E/M*z_h/g, -D/M*z_h/g, -E/M*z_h/g]))
  c_x1 = np.matrix(np.array([1, 0, 0, 0]))
  c_x2 = np.matrix(np.array([0, 0, 1, 0]))
  return A, b, c, c_x1, c_x2, (*gains(A, b, c, Qe, Qx, R, N))

def flipm_gains(m, M, g, z_h, dt, D, E, Qe, Qx, R, N):
  # Calculation of controller flipm_gains as explained in [2].
  A = np.array(
      [[      1,     dt, dt**2/2,       0,        0,       0 ],
       [      0,      1,      dt,       0,        0,       0 ],
       [   -D/M,   -E/M,       0,     D/M,      E/M,       0 ],
       [      0,      0,       0,       1,       dt, dt**2/2 ],
       [ D/m*dt, E/m*dt,       0, -D/m*dt, 1-E/m*dt,      dt ],
       [      0,      0,       0,       0,        0,       1 ]])
  b = np.matrix(np.array([0, 0, 0, dt**3 / 6, dt**2 / 2, dt])).transpose()
  c = np.matrix(np.array([1, 0, -z_h/g, 0, 0, 0]))
  c_x1 = np.matrix(np.array([1, 0, 0, 0, 0, 0]))
  c_x2 = np.matrix(np.array([0, 0, 0, 1, 0, 0]))
  return A, b, c, c_x1, c_x2, (*gains(A, b, c, Qe, Qx, R, N))
  
def flipm_gains2(m, M, g, z_h, dt, D, E, Qe, Qx, R, N):
  # New try with different integration
  DM = D/M
  EM = E/M
  Dm = D/m
  Em = E/m
  dt2 = dt**2/2

  A = np.array(
  # From vector element      ------------->       to vector element
  #                                                              ||
  #       x1       dx1     ddx1       x2       dx2     ddx2      \/
  [[1-DM*dt2,dt-EM*dt2,     dt2,  DM*dt2,   EM*dt2,       0 ], # x1
   [  -DM*dt,  1-EM*dt,      dt,   DM*dt,    EM*dt,       0 ], # dx1
   [     -DM,      -EM,       1,      DM,       EM,       0 ], # ddx1
   [  Dm*dt2,   Em*dt2,       0,1-Dm*dt2,dt-Em*dt2,     dt2 ], # x2
   [   Dm*dt,    Em*dt,       0,  -Dm*dt,  1-Em*dt,      dt ], # dx2
   [      Dm,       Em,       0,     -Dm,      -Em,       0 ]])# ddx2

  b = np.matrix(np.array([0, 0, 0, dt**3 / 6, dt**2 / 2, dt])).transpose()
  c = np.matrix(np.array([1, 0, -z_h/g, 0, 0, 0]))
  c_x1 = np.matrix(np.array([1, 0, 0, 0, 0, 0]))
  c_x2 = np.matrix(np.array([0, 0, 0, 1, 0, 0]))
  return A, b, c, c_x1, c_x2, (*gains(A, b, c, Qe, Qx, R, N))
  
def FLIPM_CP(x, dx, g, z_h, t = 0):
  om = np.sqrt(g/z_h)
  return np.exp(om * t) * (x + 1/om * dx)
  
def gains(A, b, c, Qe, Qx, R, N):
  s = b.shape[0]
  Bt = np.matrix(np.zeros((s+1, 1)))
  Bt[0, 0] = np.dot(c, b)
  Bt[1:s+1, 0] = b
  It = np.zeros((s+1,1))
  It[0, 0] = 1
  #It = np.matrix(np.array([1, 0, 0, 0, 0, 0, 0])).transpose()
  Ft = np.matrix(np.zeros((s+1, s)))
  Ft[0, 0:s] = c * A
  Ft[1:s+1, 0:s] = A
  Qt = np.matrix(np.zeros((s+1, s+1)))
  Qt[0, 0] = Qe
  Qt[1:s+1, 1:s+1] = c.transpose() * Qx *c
  At = np.matrix(np.zeros((s+1, s+1)))
  At[0:s+1, 0] = It
  At[0:s+1, 1:s+1] = Ft
  Pt = np.matrix(sp.linalg.solve_discrete_are(At, Bt, Qt, np.ones((1,1)) * R))
  Gx = (R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * Pt * Ft
  Gi = (R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * Pt * It
  Ac = At - Bt * (R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * Pt * At
  M = -Ac.transpose() * Pt * It
  Gd = [-Gi]
  for i in range(1, N):
    Gd.append((R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * M)
    M = Ac.transpose() * M
  return Gi, Gx, Gd
  
def zsp(Gx, Gd, x, N):
  # 0 = -Gx * x - sum(Gd * z)
  # -Gx * x = z * sum(Gd)
  z = -Gx * x / sum(Gd)
  z = -Gx * x / sum(Gd)
  return z[0,0]
  
class BoundedAnalytical:  
  def c1d(self, T, α, t, g, z_h):
    om = np.sqrt(g/z_h)
    c1d = 0
    for (Ti, α1) in zip(T, α):
      c1d += α1 * 0.5 * (np.exp(om * (t - Ti)) * (1 - heaviside(t - Ti)) + \
             (2 - np.exp(-om * (t - Ti))) * heaviside(t - Ti))
    return c1d
  def dc1d(self, T, α, t, g, z_h):
    om = np.sqrt(g/z_h)
    dc1d = 0
    for (Ti, α1) in zip(T, α):
      dc1d += α1 * om * 0.5 * (np.exp(om * (t - Ti)) * (1 - heaviside(t - Ti)) + \
              np.exp(-om * (t - Ti)) * heaviside(t - Ti))
    return dc1d
  def ddc1d(self, T, α, t, g, z_h):
    om = np.sqrt(g/z_h)
    ddc1d = 0
    for (Ti, α1) in zip(T, α):
      ddc1d += α1 * om**2 * 0.5 * (np.exp(om * (t - Ti)) * (1 - heaviside(t - Ti)) - np.exp(-om * (t - Ti)) * heaviside(t - Ti))
    return ddc1d
  def dddc1d(self, T, α, t, g, z_h, dt):
    om = np.sqrt(g/z_h)
    dddc1d = 0
    for (Ti, α1) in zip(T, α):
      dddc1d += α1 * om**3 * 0.5 * (np.exp(om * (t - Ti)) * (1 - heaviside(t - Ti)) + \
                np.exp(-om * (t - Ti)) * heaviside(t - Ti)) - om**2 * impulse(Ti, t, dt)
    return dddc1d
  def η(self, T, α, t, g, z_h, k, b, M):
    om = np.sqrt(g/z_h)
    η = 0
    for (Ti, α1) in zip(T, α):
      η += α1 * om**2 * 0.5 * M / b / (om + k/b) * (np.exp(om * (t - Ti)) * (1 - heaviside(t - Ti)) + np.exp(-k/b * (t - Ti)) * heaviside(t - Ti) - np.exp(-k/b * t - om * Ti)) \
     - α1 * om**2 * 0.5 * M / b / (-om + k/b) * ((np.exp(-om * (t - Ti)) - np.exp(-k/b * (t - Ti))) * heaviside(t - Ti))  
    return η
  def ud(self, T, α, t, g, z_h, k, b, M, m, dt):
    η = self.η(T, α, t, g, z_h, k, b, M)
    ddc1d = self.ddc1d(T, α, t, g, z_h)
    dddc1d  = self.dddc1d(T, α, t, g, z_h, dt)  
    return m*k**2/b**2 * η - M*m/b * (k/b-b*(m+M)/M*m) * ddc1d + M*m/b * dddc1d
    
  def bacontrol(self, A, B, C, end, m, M, g, z_h, dt, D, E, T, α, τ, plotarea):
    prefout = list()
    c1d = list()
    c2d = list()
    dc1d = list()
    ud = list()
    x1out = list()
    x2out = list()
    zmpout = list()
    y = np.matrix(0)
    
    x   = np.matrix([[self.c1d(T, α, 0, g, z_h)], 
                     [self.dc1d(T, α, 0, g, z_h)],
                     [self.c1d(T, α, 0, g, z_h)],
                     [self.dc1d(T, α, 0, g, z_h)]])
    
    X = np.linspace(0, end - dt, end*(1/dt))
    ud.append(self.ud(T, α, 0, g, z_h, D, E, M, m, dt))
    for t in X:
      c1d.append(self.c1d(T, α, t, g, z_h))
      c2d.append(c1d[-1] + self.η(T, α, t, g, z_h, D, E, M))
      ud.append(self.ud(T, α, t, g, z_h, D, E, M, m, dt))
      prefout.append(pref_lanari(T, α, t))

      x   = A  * x   + B  * ud[-1] / m
      y   = C  * x
      
      x1out.append((x[0]).item())
      x2out.append((x[2]).item())
      zmpout.append(y.item())
      
      
    
    #plotarea.plot(X, c1d, label="$c_1^d$", linestyle="dashed")
    #plotarea.plot(X, c2d, label="$c_2^d$")
    plotarea.plot(X, x1out, label="$c_{1,y}$", linestyle="dashed")
    plotarea.plot(X, x2out, label="$c_{2,y}$")
    plotarea.plot(X, zmpout, linewidth=2, label="$p_y$")
    plotarea.plot(X, prefout, linewidth=2, label="$p^{ref}_y$", linestyle="dashed")
    plotarea.legend(prop={'size':11}, borderpad=0.1)
    plotarea.set_xlabel('Time [s]')
    plotarea.set_ylabel('Position (y) [m]')
  
# A controller based on Lanari's work, iterative version
def bcontrol(A, B, C, end, m, M, g, z_h, dt, D, E, T, α, τ, plotarea):
  prefout = list()
  x1out = list()
  x2out = list()
  zmpout = list()
  xcout = list()
  ycout = list()
  
  Mb = 1/M + 1/m
  om = np.sqrt(g/z_h)
  xu0 = 0
  for e in zip(T, α):
    xu0 += e[1] * np.exp(-om * e[0])
  xs0 = 0
  Xc0 = (xu0 + xs0) / 2
  Xcdot0 = (xu0 - xs0) * om / 2
  XcStar0 = np.matrix([[Xc0], [Xcdot0]])
  x0Flex = np.matrix([[0.], [-0.]])
  
  AFlex = np.matrix([[     0,         1],
                    [-D * Mb,  - E * Mb]]);
  BFlex = np.matrix([[0],
                    [1 / m]]);
  CFlex = np.matrix([D / M - Mb * E * τ / M, 
                      E / M + τ / M * (D - Mb * E**2)]);
  DFlex = τ * E / (M * m);
  
  Ai = AFlex - BFlex*CFlex/DFlex;
  Bi = BFlex/DFlex;
  Ci = -CFlex/DFlex;
  Di = 1/DFlex;
  
  Al = np.matrix([[0, 1], [om**2, 0]])
  Bl = np.matrix([[0], [-om**2]])
  Cl = np.matrix([om**2, 0])
  Dl = np.matrix([-om**2])
  
  Ac = np.matrix([[0, 1], [0, 0]])
  Bc = np.matrix([[0], [1]])
  Cc = np.matrix([[1, 0]])
  Dc = np.matrix([-1/om**2])
  
  x_l = XcStar0.copy()
  x_i = x0Flex.copy()
  x   = np.matrix([[XcStar0[0].item()], 
                   [XcStar0[1].item()],
                   [0],
                   [XcStar0[0].item()-x0Flex[0].item()],
                   [XcStar0[1].item()-x0Flex[1].item()],
                   [0]]) 
  x_f = x0Flex.copy()
  x_c = XcStar0.copy()
  
  X = np.linspace(0, end - dt, end*(1/dt))
  X2 = np.linspace(end, end + 5, 5*(1/dt))
  Xg = np.linspace(0, end + 5, (end+5)*(1/dt))
  u_c = 0
  y_i = 0
  y_l = 0
  y_f = 0
  y_c = 0
  start = time.time()
  for t in X:
    
    if t < 1.8:
      p   = pref_lanari(T, α, t)
    if t == 1.8:
      p = FLIPM_CP(x_l[0], x_l[1], g, z_h, dt)
    
    dx_l = Al * x_l + Bl * p
    y_l_t = Cl * x_l + Dl * p
    
    dx_i = Ai * x_i + Bi * y_l
    y_i_t = Ci * x_i + Di * y_l
  
    #dx_f = AFlex * x_f + BFlex * y_i
    #y_f_t = CFlex * x_f + DFlex * y_i
    
    #dx_c = Ac * x_c + Bc * y_f
    #y_c_t = Cc * x_c + Dc * y_f
    
    u_c = y_i_t  - y_i
    
    x_l += dt * dx_l
    y_l = y_l_t
  
    x_i += dt * dx_i
    y_i = y_i_t
    
    #x_f += dt * dx_f
    #y_f = y_f_t
    
    #x_c += dt * dx_c
    #y_c = y_c_t
    
    x   = A  * x   + B  * u_c / m / dt
    y   = C  * x

    x1out.append((x[0]).item())
    x2out.append((x[3]).item())
    zmpout.append(y.item())
    #xcout.append((x_c[0]).item())
    #ycout.append(y_c.item())
    prefout.append(p)
    
  end = time.time()
  print("Time taken for bounded controll:" + str(end - start))
  cp = FLIPM_CP(x[0], x[1], g, z_h, dt)
  #for t in X2:
    
    
  plotarea.plot(X, x1out, label="$c_{1,y}^d$", linestyle="dashed")
  plotarea.plot(X, x2out, label="$c_{2,y}^d$")
  #plotarea.plot(X, xcout, label="$c_{1,y}$", linestyle="dashed")
  #plotarea.plot(X, x2out, label="$c_{2,y}$")
  plotarea.plot(X, zmpout, linewidth=2, label="$x_{zmp}$")
  #plotarea.plot(X, ycout, linewidth=2, label="$p_y$")
  plotarea.plot(X, prefout, linewidth=2, label="$x_{zmp}^d$", linestyle="dashed")
  plotarea.legend(prop={'size':11}, borderpad=0.1)
  plotarea.set_xlabel('Time [s]')
  plotarea.set_ylabel('Position (y) [m]')
    
  return
            

def control(cp_stop, pref, g, z_h, dt, N, end, plotarea, A, b, c, c_x1, c_x2, Gi, Gx, Gd): # a walk with controller
  s = b.shape[0]          
  x1out = list()
  x2out = list()
  zmpout = list()
  prefout = list()
  zspout = list()
  x = np.matrix(np.zeros((s,1)))
  X = np.linspace(0, end - dt, end*(1/dt))
  Xg = X
  v = 0
  start = time.time()
  for t in X:
    s = 0
    for t2 in range(0, N):
      s += Gd[t2] * pref(t+t2*dt)
    u = -Gi * v - Gx * x - s
    x = A * x + b * u
    v = v + c * x - pref(t)
    x1out.append((c_x1 * x).item())
    x2out.append((c_x2 * x).item())
    zmpout.append((c * x).item())
    prefout.append(pref(t))
  end = time.time()
  
  print("Time taken for preview controll:" + str(end - start))

  if cp_stop:
    X2 = np.linspace(end, end + 5, 5*(1/dt))
    Xg = np.linspace(0, end + 5, (end+5)*(1/dt))
    
    # Another loop for stopping with CP
    v = 0
    lp = prefout[-1]
    init_x = x[0] #- lp
    print("Initial position:" + str(init_x) + ", speed:" + str(x[1]))
    cp = FLIPM_CP(init_x, x[1], g, z_h, dt)#+ lp
    print("CP set to:" + str(cp))
    x[3] = x[0]
    x[4] = x[1]
    x[5] = 0  
    x[2] = 0
    for t in X2:
      s = 0
      for t2 in range(0, N):
        s += Gd[t2] * cp
      u = -Gi * v - Gx * x - s
      x = A * x + b * u
      v = v + c * x - cp
      x1out.append((c_x1 * x).item())
      x2out.append((c_x2 * x).item())
      zmpout.append((c * x).item())
      prefout.append(cp)    
    
  plotarea.plot(Xg, x1out, label="$c_{1,y}$", linestyle="dashed")
  plotarea.plot(Xg, x2out, label="$c_{2,y}$")
  plotarea.plot(Xg, zmpout, linewidth=2, label="$x_{zmp}$")
  plotarea.plot(Xg, prefout, linewidth=2, label="$x_{zmp}^d$", linestyle="dashed")
  plotarea.legend(prop={'size':11}, borderpad=0.1)
  plotarea.set_xlabel('Time [s]')
  plotarea.set_ylabel('Position (y) [m]')
  
def sim(dt, end, A, b, c, plotarea): # Just a simple demo
  x = np.matrix(np.zeros((6,1)))
  X = np.linspace(0, end - dt, end*(1/dt))
  u = np.zeros(end*(1/dt))
  u[0.5/dt] = 100 / dt
  u[0.5/dt+1] = -100 / dt
  x1 = list() # positions of mass 1
  x2 = list() # positions of mass 2

  for t in X:
    x = A * x + b * u[t*(1/dt)]
    x1.append(x[0].item())
    x2.append(x[3].item())

  plotarea.plot(X, x1, linewidth=2, label="$c_{1}$")
  plotarea.plot(X, x2, linewidth=2, label="$c_{2}$")
  plotarea.legend(prop={'size':11}, borderpad=0.1, loc="lower right")
  plotarea.set_xlabel('Time [s]')
  plotarea.set_ylabel('Position (y) [m]')
  

def pref_y(t):
  if t < 1:
    return 0
  if t - sp.floor(t) < 0.5:
    return -0.05
  else:
    return 0.05

def pref_x(t):
  if t < 1:
    return 0
  return (t // 0.5) / 10
  
def pref_lanari(T, α, t):
  s = 0
  for e in zip(T, α):
    s += e[1] * heaviside(t - e[0])
  return s
  
def heaviside(t):
  return 0.5 * (np.sign(t) + 1)




def impulse(t1, t2, dt):
  if abs(t1-t2) < dt/4:
    return 100
  else:
    return 0
  
class FLIPMApp(tkinter.Frame):
  def __init__(self, master=None):
    tkinter.Frame.__init__(self, master)
    self.createWidgets()
    self.menuBar = tkinter.Menu(master)
    master.config(menu = self.menuBar)
    self.fillMenuBar()
    self.pack()
    self.onLanariControllerDef()
    self.onAnalyticalBoundednessController()
  def addValue(self, text, min, max, inc):
    frame = tkinter.Frame(self.controlframe)
    frame.pack(fill = "x")
    txt = tkinter.Label(frame)
    txt["text"] = text
    txt.pack(side = "left", fill = "x")
    v = tkinter.DoubleVar()
    s = tkinter.Spinbox(frame, textvariable = v)
    s["to"] = max
    s["from"] = min
    s["increment"] = inc
    s.pack(side = "right", fill = "x")
    self.values[text] = v

  def fillMenuBar(self):
    # File menu
    self.menuFile = tkinter.Menu(self.menuBar)
    self.menuFile.add_command(label = "Save Gains...", command = self.onSave)
    
    # Commands menu
    self.menuCommand = tkinter.Menu(self.menuBar)
    self.menuCommand.add_command(label = "Simulate", command = self.onSim)
    self.menuCommand.add_command(label = "Load Simulation Default", command = self.onSimDef)
    self.menuCommand.add_command(label = "Controller", command = self.onController)
    self.menuCommand.add_command(label = "Load Controller Default", command = self.onControllerDef)
    self.menuCommand.add_command(label = "Load Lanari Controller Default", command = self.onLanariControllerDef)
    self.menuCommand.add_command(label = "Lanari Controller", command = self.onLanariController)
    self.menuCommand.add_command(label = "Boundedness Controller (analytical)", command = self.onAnalyticalBoundednessController)
    
    # Plotting menu
    self.menuPlot = tkinter.Menu(self.menuBar)
    self.menuPlot.add_command(label = "Plot Bounded", command = self.onPlotBounded)
    # Add everything
    self.menuBar.add_cascade(label = "File", menu = self.menuFile)
    self.menuBar.add_cascade(label = "Commands", menu = self.menuCommand)
    self.menuBar.add_cascade(label = "Plot", menu = self.menuPlot)
  def createWidgets(self):
    self.values = {}
    self.controlframe = tkinter.Frame(self)
    self.controlframe.pack(side = "left", fill = "both")
    self.addValue("Small Mass", 0.000000001, 100, 0.1)
    self.addValue("Large Mass", 0.00001, 1000000000, 0.1)
    self.addValue("Gravity", 0.1, 100, 0.1)
    self.addValue("Height", 0.01, 100, 0.1)
    self.addValue("Frame Length", 0.00001, 1, 0.001)
    self.addValue("Spring Constant", 0.001, 1000000, 0.1)
    self.addValue("Damper Constant", 0.001, 1000000, 0.1)
    self.addValue("Qx", 10**-10, 10**10, 1)
    self.addValue("Qe", 10**-10, 10**10, 1)
    self.addValue("R", 10**-10, 10**10, 1)
    self.addValue("N", 1, 1000, 1)
    self.addValue("End", 0, 100, 0.5)    
    f = Figure(figsize=(10,8), dpi=50, facecolor='white')
    self.a = f.add_subplot(111)
    self.canvas = FigureCanvasTkAgg(f, master=self)
    self.canvas.show()
    self.canvas.get_tk_widget().pack(fill=tkinter.BOTH, expand=1)
    toolbar = NavigationToolbar2TkAgg( self.canvas, self )
    toolbar.update()
    self.canvas._tkcanvas.pack(side=tkinter.TOP,  expand=1)
  def onSave(self): # m, M, g, z_h, dt, D, E, Qe, Qx, R, N
    f = tkinter.filedialog.asksaveasfile()
    v = getValues(self.values)
    g = flipm_gains(*v[:11]) # A, b, c, Gi, Gx, Gd
    s = \
    """
    A = {{ content = [{}]; }};
    b = {{ content = [{}]; }};
    c = {{ content = [{}]; }};
    Gi = {};
    Gx = {{ content = [{}]; }};
    Gd = {{ content = [{}]; }};
    m = {};
    M = {};
    g = {};
    z_h = {};
    dt = {};
    D = {};
    E = {};
    Qe = {};
    Qx = {};
    R = {};
    N = {};
    """.format(",".join(map(str, map(float, g[0].flatten()))),
               ",".join(map(str, map(float, g[1]))),
               ",".join(map(str, map(float, g[2].tolist()[0]))),
               float(g[3]),
               ",".join(map(str, map(float, g[4].flatten().tolist()[0]))),
               ",".join(map(str, map(float, g[5]))),
               *(map(str, v)))
    f.write(s)
    f.close()
  def onSim(self):
    self.a.cla()
    A, b, c, *r = flipm_gains(*getValues(self.values)[:11])
    sim(self.values["Frame Length"].get(), self.values["End"].get(), A, b, c, self.a)
    self.canvas.show()
  def onControl3D(self):
    self.startgl()
    
  def onController(self, export = False):
    self.a.cla();
    control(0, pref_y, self.values["Gravity"].get(),
            self.values["Height"].get(), 
            self.values["Frame Length"].get(),
            int(sp.floor(self.values["N"].get())),
            self.values["End"].get(), self.a,
            *flipm_gains(*getValues(self.values)[:11]))
    self.canvas.show()
  def onLanariController(self, export = False):
    self.a.cla();
    A, b, c, *r = flipm_gains(*getValues(self.values)[:11])
    bcontrol(A, b, c, 
             self.values["End"].get(),
             self.values["Small Mass"].get(),
             self.values["Large Mass"].get(),
             self.values["Gravity"].get(),
             self.values["Height"].get(),
             self.values["Frame Length"].get(),
             self.values["Spring Constant"].get(),
             self.values["Damper Constant"].get(),
             (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
             (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
             self.values["Qe"].get(),
             self.a)
    self.canvas.show()
    
  def onAnalyticalBoundednessController(self, export = False):
    self.a.cla();
    ba = BoundedAnalytical()
    A, b, c, *r = flipm_gains4D(*getValues(self.values)[:11])
    ba.bacontrol(A, b, c, 
             self.values["End"].get(),
             self.values["Small Mass"].get(),
             self.values["Large Mass"].get(),
             self.values["Gravity"].get(),
             self.values["Height"].get(),
             self.values["Frame Length"].get(),
             self.values["Spring Constant"].get(),
             self.values["Damper Constant"].get(),
             (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
             (-0.05, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1),
             #list([0.5]),
             #list([-0.05]),
             self.values["Qe"].get(),
             self.a)
    self.canvas.show()
  def onControllerDef(self):
    self.setValues(*controllerDefault)
  def onLanariControllerDef(self):
    self.setValues(*lanariControllerDefault)
    
  def onPlotBounded(self):
      plotBounded()
  def onSimDef(self):
    self.setValues(*simDefault)
  def setValues(self, m, M, g, z_h, dt, D, E, Qe, Qx, R, N, end):
    self.values["Small Mass"].set(m)
    self.values["Large Mass"].set(M)
    self.values["Gravity"].set(g)
    self.values["Height"].set(z_h)
    self.values["Frame Length"].set(dt)
    self.values["Spring Constant"].set(D)
    self.values["Damper Constant"].set(E)
    self.values["Qe"].set(Qe)
    self.values["Qx"].set(Qx)
    self.values["R"].set(R)
    self.values["N"].set(N)
    self.values["End"].set(end)
    
main()