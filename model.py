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
from tkinter import ttk
import tkinter.filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import sys

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

simDefault = (0.1, 1, 9.81, 0.3, 0.01, 1, 0.1, 1, 1, 10**-10,
              np.matrix(np.array([[1, 0, 0, 0 ,0 ,0 ],
                                  [0, 1, 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0],
                                  [0, 0, 0, 1, 0, 0],
                                  [0, 0, 0, 0, 1, 0],
                                  [0, 0, 0, 0 ,0 ,1]])),
              np.matrix(np.array([[1,  0, 0],
                                  [0,  1, 0],
                                  [0,  0, 1]])),              
              100, 5)
controllerDefault = (0.1, 4.5, 9.81, 0.26, 0.01, 5000, 200, 1, 10, 10**-10,
                     np.matrix(np.array([[10**-0, 0, 0, 0 ,0 ,0 ],
                                         [0, 1, 0, 0, 0, 0],
                                         [0, 0, 10**2, 0, 0, 0],
                                         [0, 0, 0, 10**2, 0, 0],
                                         [0, 0, 0, 0, 1, 0],
                                         [0, 0, 0, 0 ,0 ,1]])),
                     np.matrix(np.array([[10**17,  0, 0],
                                         [0,  10**4, 0],
                                         [0,  0, 10**4]])),
                     100, 5)
error = np.matrix([0.0, 0.0, 0.0]).transpose()

def main():
  noGui = False
  for arg in sys.argv:
    if arg == '-c': # Plot figure as shown in [1]
      fig1 = plt.figure(figsize=(6, 2), dpi=150)
      ax1 = fig1.add_subplot(111)
      control(controllerDefault[4], int(sp.floor(controllerDefault[10])), controllerDefault[11], ax1, *gains(*controllerDefault[:11]))
      plt.savefig("paper/build/FLIPMex" + ".pdf", format='pdf', bbox_inches='tight')
      noGui = True
    if arg == '-r': # Plot another figure shown in [1]
      fig1 = plt.figure(figsize=(6, 2), dpi=150)
      ax1 = fig1.add_subplot(111)
      g = gains(*simDefault[:11])
      sim(simDefault[4],simDefault[11], g[0], g[1], g[2], ax1)
      plt.savefig("paper/build/FLIPMmodel" + ".pdf", format='pdf', bbox_inches='tight')
      noGui = True
    if arg.startswith('-l='): # Print gain tables in latex
      m.dump(open(arg[3:] + "A.tex",'w'))
      noGui = True
  if noGui:
    exit()
  root = tkinter.Tk()
  app = FLIPMApp(root)
  app.mainloop()
 
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
          d["Ql"].getMatrix(),
          d["RO"].getMatrix(),
          int(sp.floor(d["N"].get())),
          d["End"].get())
def gains(m, M, g, z_h, dt, D, E, Qe, Qx, R, Ql, RO, N):
  # Calculation of controller gains as explained in [2].
  A = np.array(
      [[      1,     dt, dt**2/2,       0,        0,       0 ],
       [      0,      1,      dt,       0,        0,       0 ],
       [   -D/M,   -E/M,       0,     D/M,      E/M,       0 ],
       [      0,      0,       0,       1,       dt, dt**2/2 ],
       [ D/m*dt, E/m*dt,       0, -D/m*dt, 1-E/m*dt,      dt ],
       [      0,      0,       0,       0,        0,       1 ]])
  b = np.matrix(np.array([0, 0, 0, dt**3 / 6, dt**2 / 2, dt])).transpose()
  c = np.matrix(np.array([1, 0, -z_h/g, 0, 0, 0]))
  Bt = np.matrix(np.zeros((7, 1)))
  Bt[0, 0] = np.dot(c, b)
  Bt[1:7, 0] = b
  It = np.matrix(np.array([1, 0, 0, 0, 0, 0, 0])).transpose()
  Ft = np.matrix(np.zeros((7, 6)))
  Ft[0, 0:6] = c * A
  Ft[1:7, 0:6] = A
  Qt = np.matrix(np.zeros((7, 7)))
  Qt[0, 0] = Qe
  Qt[1:7, 1:7] = c.transpose() * Qx *c
  At = np.matrix(np.zeros((7, 7)))
  At[0:7, 0] = It
  At[0:7, 1:7] = Ft
  Pt = np.matrix(sp.linalg.solve_discrete_are(At, Bt, Qt, np.ones((1,1)) * R))
  Gx = (R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * Pt * Ft
  Gi = (R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * Pt * It
  Ac = At - Bt * (R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * Pt * At
  X = -Ac.transpose() * Pt * It
  Gd = [-Gi]
  for i in range(1, N):
    Gd.append((R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * X)
    X = Ac.transpose() * X

  ##########################
  # Observer implementation
  ##########################
  Cm = np.matrix(np.zeros((3, 6)))  # Was ich von dem Zustand messen kann
  Cm[0,0] = 1
  Cm[1,2] = 1    
  Cm[2,3] = 1  
  
  
  L,S,e = dlrq(A.transpose(),Cm.transpose(), Ql, RO)
  L = L.transpose()
  return A, b, c, Gi, Gx, Gd, L

def dlrq(A,B,Q,R) :                                         # http://de.mathworks.com/help/control/ref/dlqr.html
  S = np.matrix(sp.linalg.solve_discrete_are(A, B, Q, R))   # solve the discrete-time Riccati equation
  K = (B.transpose()*S*B + R) ** -1 * (B.transpose()*S*A)   # calculate gain matrix k
  e = np.matrix(sp.linalg.eigvals(A - B*K))                 # calculate eigenvalues
  return K,S,e

def control(dt, N, end, error, plotarea, A, b, c, Gi, Gx, Gd, L): # a walk with controller
  x1out = list()
  x2out = list()
  zmpout = list()
  prefout = list()
  x = np.matrix(np.zeros((6,1)))
  X = np.linspace(0, end - dt, end*(1/dt))
  v = 0
  for t in X:
    s = 0
    for t2 in range(0, N):
      s += Gd[t2] * pref(t+t2*dt)
    u = -Gi * v - Gx * x - s
    x = A * x + b * u + L*e(error,t)
    v = v + c * x - pref(t)
    x1out.append(x[0].item())
    x2out.append(x[3].item())
    zmpout.append((c * x).item())
    prefout.append(pref(t))
  plotarea.plot(X, x1out, label="$c_{1,y}$", linestyle="dashed")
  plotarea.plot(X, x2out, label="$c_{2,y}$")
  plotarea.plot(X, zmpout, linewidth=2, label="$p_y$")
  plotarea.plot(X, prefout, linewidth=2, label="$p^{ref}_y$", linestyle="dashed")
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

def e(error, t): #Error
  if t >=2 and t<2.1:
    return error
  else:
    return np.matrix([0.0, 0.0, 0.0]).transpose()

def pref(t):
  if t < 1:
    return 0
  if t - sp.floor(t) < 0.5:
    return -0.05
  else:
    return 0.05
    
class StepControl:
  def __init__(self, N, A, b, c, Gi, Gx, Gd):
    self.N = N
    self.A = A
    self.b = b
    self.c = c
    self.Gi = Gi
    self.Gx = Gx
    self.Gd = Gd
    x = np.matrix(np.zeros((6,1)))
  def step(pref):
    for i in range(0, N):
      s += Gd[i] * pref[i]
    u = -Gi * v - Gx * x - s
    x = A * x + b * u
    v = v + c * x - pref(0)

class MatrixInput(tkinter.Frame):
    def __init__(self, observerframe, text, m, n, min, max, inc):
        frame = tkinter.Frame(observerframe)
        self.matrix = {}
        self.m = m
        self.n = n
        
        for m in range(self.m+1):
            if m==0:
              tkinter.Label(frame, text=text, width= 8,  font = ('bold') ).grid(row = m,  column=0)
            else:
              tkinter.Label(frame, text=m, width= 8).grid(row = m, column=0)
        for n in range(self.n+1):
            if n!=0:
              tkinter.Label(frame, text=n, width= 8).grid(row = 0, column=n)      
              
        for m in range(self.m):
            for n in range(self.n):
                v = tkinter.DoubleVar()
                s = tkinter.Spinbox(frame, textvariable = v, width = 10)
                s["to"] = max
                s["from"] = min
                s["increment"] = inc
                s.grid(row=m+1, column=n+1, stick="nsew")
                self.matrix[(m,n)] = v
        frame.pack(fill = "both")

    def getMatrix(self):
        result = []
        for m in range(self.m):
            temp_row = []
            for n in range(self.n):
                temp_row.append(self.matrix[(m,n)].get())
            result.append(temp_row)
        return np.matrix(result)
    
    def setMatrix(self, setMatrix):
        for m in range(self.m):
            for n in range(self.n):
                self.matrix[(m,n)].set(setMatrix[m,n])
        
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
    self.menuFile = tkinter.Menu(self.menuBar)
    self.menuFile.add_command(label = "Save Gains...", command = self.onSave)
    self.menuCommand = tkinter.Menu(self.menuBar)
    self.menuCommand.add_command(label = "Simulate", command = self.onSim)
    self.menuCommand.add_command(label = "Load Simulation Default", command = self.onSimDef)
    self.menuCommand.add_command(label = "Controller", command = self.onController)
    self.menuCommand.add_command(label = "Load Controller Default", command = self.onControllerDef)
    self.menuBar.add_cascade(label = "File", menu = self.menuFile)
    self.menuBar.add_cascade(label = "Commands", menu = self.menuCommand)
  def createWidgets(self):
    notebook = ttk.Notebook(self)
    notebook.pack()
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
    self.observerframe = tkinter.Frame(self)
    self.observerframe.pack(side = "left", fill = "both")
    self.values["Ql"] = MatrixInput(self.observerframe,"Ql", 6, 6, 10**10, 10**10, 100)
    self.values["RO"] = MatrixInput(self.observerframe,"RO", 3, 3, 10**-10, 10**20, 100)
    self.addValue("N", 1, 1000, 1)
    self.addValue("End", 0, 100, 0.5)
    self.values["e"] = MatrixInput(self.observerframe,"Error", 3, 1, -10**3, 10**3, 0.001)
    ###########
    # Buttons #
    ###########
    frame = tkinter.Frame(self.observerframe)
    frame.pack(fill = "x")
    b = tkinter.Button(frame, text="Load Controller Default", command=self.onControllerDef)
    b.pack(fill = "x")
    b = tkinter.Button(frame, text="Simulate", command=self.onController)
    b.pack(fill = "x")
    frame = tkinter.Frame(self.controlframe)
    frame.pack(fill = "x")
    b = tkinter.Button(frame, text="Load Controller Default", command=self.onControllerDef)
    b.pack(fill = "x")
    b = tkinter.Button(frame, text="Simulate", command=self.onController)
    b.pack(fill = "x")
    ###########
    # Buttons #
    ###########       
    f = Figure(figsize=(16,12), dpi=50, facecolor='white')
    self.a = f.add_subplot(111)
    self.canvas = FigureCanvasTkAgg(f, master=self)
    self.canvas.show()
    self.canvas.get_tk_widget().pack(fill=tkinter.BOTH, expand=1)
    toolbar = NavigationToolbar2TkAgg( self.canvas, self )
    toolbar.update()
    self.canvas._tkcanvas.pack(side=tkinter.TOP,  expand=1)
    notebook.add(self.controlframe, text="Controller", state="normal")
    notebook.add(self.observerframe, text="Observer", state="normal")
    notebook.pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=1)
  def onSave(self): # m, M, g, z_h, dt, D, E, Qe, Qx, R, Ql, R0, N
    f = tkinter.filedialog.asksaveasfile()
    v = getValues(self.values)
    g = gains(*v[:13]) # A, b, c, Gi, Gx, Gd, L
    
    s = \
"""
A = {{ content = [{}]; }};
b = {{ content = [{}]; }};
c = {{ content = [{}]; }};
Gi = {};
Gx = {{ content = [{}]; }};
Gd = {{ content = [{}]; }};
L = {{ content = [{}]; }};
error = {{ content = [{}]; }};
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
Ql = {{ content = [{}]; }};
RO = {{ content = [{}]; }};
N = {};
""".format(",".join(map(str, map(float, g[0].flatten()))),
           ",".join(map(str, map(float, g[1]))),
           ",".join(map(str, map(float, g[2].tolist()[0]))),
           float(g[3]),
           ",".join(map(str, map(float, g[4].flatten().tolist()[0]))),
           ",".join(map(str, map(float, g[5]))),
           ",".join(map(str, map(float, np.array(g[6]).flatten()))),
           ",".join(map(str, map(float, self.values["e"].getMatrix()))),
           *(map(str, v)))
    f.write(s)
    f.close()
  def onSim(self):
    self.a.cla()
    A, b, c, *r = gains(*getValues(self.values)[:13])
    sim(self.values["Frame Length"].get(), self.values["End"].get(), A, b, c, self.a)
    self.canvas.show()
  def onControl3D(self):
    self.startgl()
    
  def onController(self, export = False):
    self.a.cla();
    control(self.values["Frame Length"].get(), int(sp.floor(self.values["N"].get())),  self.values["End"].get(), self.values["e"].getMatrix(), self.a, *gains(*getValues(self.values)[:13]))
    self.canvas.show()
  def onControllerDef(self):
    self.setValues(*controllerDefault)
  def onSimDef(self):
    self.setValues(*simDefault)
  def setValues(self, m, M, g, z_h, dt, D, E, Qe, Qx, R, Ql, RO, N, end):
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
    self.values["Ql"].setMatrix(Ql)
    self.values["RO"].setMatrix(RO)
    self.values["N"].set(N)
    self.values["End"].set(end)
    
main()