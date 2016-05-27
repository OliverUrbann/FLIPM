# Author:       Oliver Urbann
# E-Mail:       oliver.urbann@tu-dortmund.de
# Requirements: Python 3.4+, scipy, numpy, matplotlib, tkinter
# Usage:        python model.py
#
# Note:         It can be quite tricky to get this script fully 
#               running as most distributions use buggy packages
#               and versions. In case of crashes use the recommended
#               distro version.
#               Recommended for OS X: Anaconda 3.18.8 x64
#               Recommended for Windows: Anaconda 
#
# This script is the implementation of the FLIP model utilized by a 
# preview controller as proposed by Urbann et al. [1]. Details about 
# the gains used in the GUI can be found in [1] and [2].
#
# Additionally, it contains also the Boundedness controller as 
# proposed by Lanari et al. [3] that utilizes the same FLIPM. This
# controller can be executed as an iterative version including a
# capture step and as an analytical version.
#
# With this tool, you can run a simulation of the flexible system
# (Commands->Simulate), run a preview controller 
# (Commands->Preview Controller) and run the Boundedness Controller
# in the iterative version (Commands->Boundedness Controller (iterative))
# and in the analytical closed form 
# (Commands->Boundedness Controller (analytical)). 
# You can chose your own gains on the left side but I would recommend
# to load some corresponding default values using Commands->Load...
#
# All figures except experiments in simulation or on physical robot
# in [1] and [3] where done using this python tool. You can plot
# them on your own with Plot->Plot for paper [1] and
# Plot->Plot for paper [3]. This way it should be easier for you to
# see how the controller and model was used to get these figures.
# The result will be pdf documents in the working folder.
# 
# Additionally, the gains and the computed controller gains can be
# saved for an application of the controller implemented on a 
# physical robot.
#
#
# [1] Flexible Linear Inverted Pendulum Model for Cost-Effective 
#     Biped Robots
#     Oliver Urbann, Ingmar Schwarz, Matthias Hofmann
#     Humanoid Robots (Humanoids), 2015 15th IEEE-RAS International
#     Conference on, 2015, pp. 128–131 
# [2] Observer-based biped walking control, a sensor fusion approach
#     Oliver Urbann, Stefan Tasse
#     Autonomous Robots 35.1 (2013) pp. 37–49. Springer US, 2013
# [3] Boundedness Approach to Gait Planning for the Flexible Linear 
#     Inverted Pendulum Model
#     Leonardo Lanari, Oliver Urbann, Seth Hutchinson, Ingmar Schwarz
#     RoboCup 2016: Robot World Cup XX, 2017, to appear

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
__copyright__ = "Copyright 2016, Oliver Urbann"
__credits__ = ["Oliver Urbann"]
__license__ = "GPL"
__version__ = "2.0.0"
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

# Some default values. Format is:
# 
simDefault = (0.1, 1, 9.81, 0.3, 0.01, 1, 0.1, 1, 1, 10**-10, 100, 5)
boundednessControllerDefault = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.001, 1000, 10**-10, 100, 3)
previewControllerDefault = (0.1, 4.5, 9.81, 0.26, 0.01, 5000, 200, 1, 10, 10**-10, 100, 5)

# For tests with CP
# previewControllerDefault = (0.1, 4.5, 9.81, 0.26, 0.005, 10000, 200, 10, 10, 10**-10, 200, 5)


def main():
  noGui = False
  for arg in sys.argv:
    if arg.startswith('-l='): # Print gain tables in latex
      m.dump(open(arg[3:] + "A.tex",'w'))
      noGui = True
  if noGui:
    exit()
  root = tkinter.Tk()
  app = FLIPMApp(root)
  app.mainloop()
 

# This plots the figures of [3] that are obtained using this python code. 
# The data of the simulations using ODE can be found in the sub directory
# BoundednessController.
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
  
  # This will take a lot of time. High frequency of 2 kHz ensures low tracking
  # error but preview contains a lot of elements then.
  
  print("preview.pdf")
  fig1 = plt.figure(figsize=(5, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  params = (1, 5, 9.81, 0.5, 0.0005, 1000, 200, 0.1, 0.001, 10**-10, 2000, 3)
  # Preview 2000 and dt = 0.0005
  A, b, c, *r = flipm_gains(*params[:11])
  control(0, pref_y, params[2], params[3], params[4], params[10],
          params[11], ax1, *flipm_gains(*params[:11]))
  plt.savefig("preview" + ".pdf", format='pdf', bbox_inches='tight')
  
  
  

# Plot the figures of [1] without the figures of the experiments with
# physical robot.
def plotPreview():
  fig1 = plt.figure(figsize=(6, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  control(0, pref_y, previewControllerDefault[2], 
                     previewControllerDefault[3], 
                     previewControllerDefault[4],                 
                     int(sp.floor(previewControllerDefault[10])),
                     previewControllerDefault[11], ax1, 
                     *flipm_gains(*previewControllerDefault[:11]))
  plt.savefig("FLIPMex" + ".pdf", format='pdf', bbox_inches='tight')
  
  fig1 = plt.figure(figsize=(6, 2), dpi=150)
  ax1 = fig1.add_subplot(111)
  g = flipm_gains(*simDefault[:11])
  sim(simDefault[4],simDefault[11], g[0], g[1], g[2], ax1)
  plt.savefig("FLIPMmodel" + ".pdf", format='pdf', bbox_inches='tight')
  
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
          

# State space matrices for a usual LIPM
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
  
# State space matrices for FLIPM based on the LIP. In contrast to flipm_gains()
# here the ZMP is part of the state and actively controlled (experimental).
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



# Reduced version of the FLIPM with just 4 dimensions. Used for the analytical
# version of the Boundedness Controller
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


# State space matrices for FLIPM as proposed in [1]. Used for iterative version
# of Boundedness Controller and preview controller.
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
  
# Like flipm_gains() but with different integration method (experimental).
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
  

# Determine the position for a capture step based on FLIPM and the
# Boundedness Controller [3].
def FLIPM_CP(x, dx, g, z_h, t = 0):
  om = np.sqrt(g/z_h)
  return np.exp(om * t) * (x + 1/om * dx)
  

# Calculates the gain matrices for the preview controller and observer.
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
  
  """ 
    This function shows the analytical version of the Boundedness Controller. The output c2d of it is used for the simulation experiments in [3]. Here additionally ud is calculated which is usually used for state space representations x=Ax+Bu. As can be seen, using ud to control the system is not appropriate as the output y reveals the typical form of integrational errors. Thus the approximate system is required to avoid the impulse in dddc1d. Note that for controlling the robot in simulation ud is not required. As it is a position controlled robot, only c2d is required. In contrast to preview control, it is therefore not required to do the integration step x=Ax+Bu using Boundedness Control.
  """
  def bacontrol(self, A, B, C, end, m, M, g, z_h, dt, D, E, T, α, τ, plotarea):
    # Some empty lists, later filled for plotting
    prefout = list() # The used ZMP reference, here not used, just for plotting
    c1d = list() # Position of the large mass, we expect that this will be
                 # the position of the CoM of the robot if we use c2d as reference
    c2d = list() # Position of the small mass that is actively controller and
                 # applied to the robot
    dc1d = list() # Speed of c1d
    ud = list() # The output of the controller. On our physical robot we do
                # not use this value, we directly apply c2d 
    x1out = list() # Position of large CoM  of a simulation based on FLIPM
    x2out = list() # Simulated position of small CoM
    zmpout = list() # The resulting ZMP
    y = np.matrix(0)
    
    # Set the initial conditions
    x   = np.matrix([[self.c1d(T, α, 0, g, z_h)], 
                     [self.dc1d(T, α, 0, g, z_h)],
                     [self.c1d(T, α, 0, g, z_h)],
                     [self.dc1d(T, α, 0, g, z_h)]])
    
    # Begin a loop over the time to be simulated
    X = np.linspace(0, end - dt, end*(1/dt))
    ud.append(self.ud(T, α, 0, g, z_h, D, E, M, m, dt))
    for t in X:
      # Add data for later plotting
      c1d.append(self.c1d(T, α, t, g, z_h))
      c2d.append(c1d[-1] + self.η(T, α, t, g, z_h, D, E, M))
      
      # This appends the output of the actual tun of the controller
      ud.append(self.ud(T, α, t, g, z_h, D, E, M, m, dt))
      
      # Just for plotting
      prefout.append(pref_boundedness(T, α, t))

      # ud[-1] is the controller value of the current time frame. Next two lines
      # simulate a walk of the robot that is modelled here using FLIPM. It is
      # not comparable to a dynamics or rigid body simulation like ODE but a
      # good way to test the system. 
      x   = A  * x   + B  * ud[-1] / m
      y   = C  * x
      
      x1out.append((x[0]).item())
      x2out.append((x[2]).item())
      zmpout.append(y.item())
      
    plotarea.plot(X, c1d, label="$c_1^d$", linestyle="dashed")
    plotarea.plot(X, c2d, label="$c_2^d$")
    plotarea.plot(X, x1out, label="$c_{1,y}$", linestyle="dashed")
    plotarea.plot(X, x2out, label="$c_{2,y}$")
    plotarea.plot(X, zmpout, linewidth=2, label="$p_y$")
    plotarea.plot(X, prefout, linewidth=2, label="$p^{ref}_y$", linestyle="dashed")
    plotarea.legend(prop={'size':11}, borderpad=0.1)
    plotarea.set_xlabel('Time [s]')
    plotarea.set_ylabel('Position (y) [m]')
  
# A controller based on [3], iterative version
def bcontrol(A, B, C, end, m, M, g, z_h, dt, D, E, T, α, τ, plotarea):
  prefout = list() # The used ZMP reference, here not used, just for plotting
  x1out = list() # Position of large CoM  of a simulation based on FLIPM
  x2out = list() # Simulated position of small CoM
  zmpout = list() # The resulting ZMP

  # Some calculations to set up the controller
  Mb = 1/M + 1/m
  om = np.sqrt(g/z_h)
  xu0 = 0
  for e in zip(T, α):
    xu0 += e[1] * np.exp(-om * e[0])
  xs0 = 0
  
  # Set initial conditions for controller (state etc.)
  Xc0 = (xu0 + xs0) / 2
  Xcdot0 = (xu0 - xs0) * om / 2
  XcStar0 = np.matrix([[Xc0], [Xcdot0]])
  x0Flex = np.matrix([[0.], [-0.]])
  
  # The flexible part of the flexible-cart-series.
  # We start with the description of the flexible system.
  AFlex = np.matrix([[     0,         1],
                    [-D * Mb,  - E * Mb]]);
  BFlex = np.matrix([[0],
                    [1 / m]]);
  CFlex = np.matrix([D / M - Mb * E * τ / M, 
                      E / M + τ / M * (D - Mb * E**2)]);
  DFlex = τ * E / (M * m);
  
  # Now the inversion to retrieve to controller for the flexible part
  Ai = AFlex - BFlex*CFlex/DFlex;
  Bi = BFlex/DFlex;
  Ci = -CFlex/DFlex;
  Di = 1/DFlex;
  
  # Here the LIP system is set up. It is the first part of the controller.
  Al = np.matrix([[0, 1], [om**2, 0]])
  Bl = np.matrix([[0], [-om**2]])
  Cl = np.matrix([om**2, 0])
  Dl = np.matrix([-om**2])

  x_l = XcStar0.copy()
  x_i = x0Flex.copy()
  
  # Initialization of the state for the simulation. It must be initialized
  # with the same values used for the initial state of the controller.
  x   = np.matrix([[XcStar0[0].item()], 
                   [XcStar0[1].item()],
                   [0],
                   [XcStar0[0].item()-x0Flex[0].item()],
                   [XcStar0[1].item()-x0Flex[1].item()],
                   [0]]) 
  x_f = x0Flex.copy()
  x_c = XcStar0.copy()
  
  X = np.linspace(0, end - dt, end*(1/dt))
  u_c = 0
  y_i = 0
  y_l = 0
  y_f = 0
  start = time.time()
  
  # Now the control/simulation loop
  for t in X:

    # This simulates an unforeseen and instant stop at 1.8s.
    if t < 1.8:
      p   = pref_boundedness(T, α, t)
    if t == 1.8:
      p = FLIPM_CP(x_l[0], x_l[1], g, z_h, dt)
    
    # First: the LIP part of the controller
    dx_l = Al * x_l + Bl * p
    y_l_t = Cl * x_l + Dl * p
    
    # Second: the flexible part of the controller
    dx_i = Ai * x_i + Bi * y_l
    y_i_t = Ci * x_i + Di * y_l
   
    # Determine the control output. It is the derivative of y_i
    u_c = y_i_t  - y_i
    
    # Control steps are done, now integrate for next loop.
    x_l += dt * dx_l
    y_l = y_l_t
  
    x_i += dt * dx_i
    y_i = y_i_t
      
    # Next two lines
    # simulate a walk of the robot that is modelled here using FLIPM. It is
    # not comparable to a dynamics or rigid body simulation like ODE but a
    # good way to test the system. 
    x   = A  * x   + B  * u_c / m / dt
    y   = C  * x

    x1out.append((x[0]).item())
    x2out.append((x[3]).item())
    zmpout.append(y.item())
    prefout.append(p)
    
  end = time.time()
  print("Time taken for bounded controll:" + str(end - start))
  cp = FLIPM_CP(x[0], x[1], g, z_h, dt)

  plotarea.plot(X, x1out, label="$c_{1,y}^d$", linestyle="dashed")
  plotarea.plot(X, x2out, label="$c_{2,y}^d$")
  plotarea.plot(X, zmpout, linewidth=2, label="$x_{zmp}$")
  plotarea.plot(X, prefout, linewidth=2, label="$x_{zmp}^d$", linestyle="dashed")
  plotarea.legend(prop={'size':11}, borderpad=0.1)
  plotarea.set_xlabel('Time [s]')
  plotarea.set_ylabel('Position (y) [m]')
    
  return
            

# A FLIPM walk with a preview controller
def control(cp_stop, pref, g, z_h, dt, N, end, plotarea, A, b, c, c_x1, c_x2, Gi, Gx, Gd): 
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
  


# Returns the referene ZMP at time t for y axis.
def pref_y(t):
  if t < 1:
    return 0
  if t - sp.floor(t) < 0.5:
    return -0.05
  else:
    return 0.05

# Returns the referene ZMP at time t for x axis.
def pref_x(t):
  if t < 1:
    return 0
  return (t // 0.5) / 10
  

# Converts the given steps to a ZMP trajectory. Better for plotting.
def pref_boundedness(T, α, t):
  s = 0
  for e in zip(T, α):
    s += e[1] * heaviside(t - e[0])
  return s
  
def heaviside(t):
  return 0.5 * (np.sign(t) + 1)





# Note that this is not a valid implementation of an impulse and as the
# example of the analytical boundedness controller shows, it is also not
# a useful approximation.
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
    
    self.menuCommand.add_command(label = "Load Simulation Default", command = self.onSimDef)
    self.menuCommand.add_command(label = "Simulate", command = self.onSim)
    self.menuCommand.add_separator()
    
    self.menuCommand.add_command(label = "Load Preview Controller Default", command = self.onControllerDef)
    self.menuCommand.add_command(label = "Preview Controller", command = self.onController)
    self.menuCommand.add_separator()
    
    self.menuCommand.add_command(label = "Load Boundedness Controller Default", command = self.onBoundednessControllerDef)
    self.menuCommand.add_command(label = "Boundedness Controller (iterative)", command = self.onBoundednessController)
    self.menuCommand.add_command(label = "Boundedness Controller (analytical)", command = self.onAnalyticalBoundednessController)
    
    # Plotting menu
    self.menuPlot = tkinter.Menu(self.menuBar)
    self.menuPlot.add_command(label = "Plot for paper [3] (Boundedness Control)", command = self.onPlotBounded)
    self.menuPlot.add_command(label = "Plot for paper [1] (FLIPM with Preview Control)", command = self.onPlotPreview)
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
       
  def onController(self, export = False):
    self.a.cla();
    control(0, pref_y, self.values["Gravity"].get(),
            self.values["Height"].get(), 
            self.values["Frame Length"].get(),
            int(sp.floor(self.values["N"].get())),
            self.values["End"].get(), self.a,
            *flipm_gains(*getValues(self.values)[:11]))
    self.canvas.show()
  def onBoundednessController(self, export = False):
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
             self.values["Qe"].get(),
             self.a)
    self.canvas.show()
  def onControllerDef(self):
    self.setValues(*previewControllerDefault)
  def onBoundednessControllerDef(self):
    self.setValues(*boundednessControllerDefault)
    
  def onPlotBounded(self):
      plotBounded()
  def onPlotPreview(self):
      plotPreview()
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