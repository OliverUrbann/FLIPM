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
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *

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
controllerDefault = (0.1, 4.5, 9.81, 0.26, 0.01, 5000, 200, 1, 10, 10**-10, 100, 5)

def main():
  noGui = False
  for arg in sys.argv:
    if arg == '-c':
      fig1 = plt.figure(figsize=(6, 2), dpi=150)
      ax1 = fig1.add_subplot(111)
      control(controllerDefault[4], int(sp.floor(controllerDefault[10])), controllerDefault[11], ax1, *gains(*controllerDefault[:11]))
      plt.savefig("paper/build/FLIPMex" + ".pdf", format='pdf', bbox_inches='tight')
      noGui = True
    if arg == '-r':
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
          int(sp.floor(d["N"].get())),
          d["End"].get())
def gains(m, M, g, z_h, dt, D, E, Qe, Qx, R, N):
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
  M = -Ac.transpose() * Pt * It
  Gd = [-Gi]
  for i in range(1, N):
    Gd.append((R + Bt.transpose() * Pt * Bt) ** -1 * Bt.transpose() * M)
    M = Ac.transpose() * M
  return A, b, c, Gi, Gx, Gd

def control(dt, N, end, plotarea, A, b, c, Gi, Gx, Gd): # a walk with controller
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
    x = A * x + b * u
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
        
class FLIPMApp(tkinter.Frame):
  def timer(self, v):
    global rotx,roty
    #roty += 1
    glutPostRedisplay()
    glutTimerFunc(10,self.timer,10)
    return
    
  def display(self):
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    gluLookAt(1,1,1,0,0,0,0,0,1)
    #glRotatef(roty,0,1,0)
    #glRotatef(rotx,1,0,0)
    glCallList(1)
    glutSwapBuffers()
    return
    
  def startgl(self):
    name = "FLIPM"
    height = 400
    width = 400

    #glutInit(name)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(height,width)
    glutCreateWindow(name)
    glClearColor(0.,0.,0.,1.)
    
    # setup display list
    glNewList(1,GL_COMPILE)
    glPushMatrix()
    glScalef(2.,1.,1.)
    glutSolidCube(0.05)
    glPopMatrix()
    glTranslatef(0.,0,0.27) #move to where we want to put object
    glutSolidSphere(0.05,10,10) # make radius 1 sphere of res 10x10
    glEndList()
    
    #setup lighting
    glEnable(GL_CULL_FACE)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    lightZeroPosition = [10.,4.,10.,1.]
    lightZeroColor = [0.8,1.0,0.8,1.0] # greenish
    glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
    glEnable(GL_LIGHT0)
    
    #setup cameras
    glMatrixMode(GL_PROJECTION)
    gluPerspective(40.,1.,1.,40.)
    glMatrixMode(GL_MODELVIEW)
    gluLookAt(0,0,10,0,0,0,0,1,0)
    glPushMatrix()
    
    #setup callbacks
    glutDisplayFunc(self.display)
    glutTimerFunc(10,self.timer,1)
    self.controller = [StepControl(int(sp.floor(self.values["N"].get())), *gains(*getValues(self.values)[:11]))]
    
    glutMainLoop()
    
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
    #self.menuCommand.add_command(label = "Controller 3D", command = self.onControl3D)
    self.menuCommand.add_command(label = "Load Controller Default", command = self.onControllerDef)
    self.menuBar.add_cascade(label = "File", menu = self.menuFile)
    self.menuBar.add_cascade(label = "Commands", menu = self.menuCommand)
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
    g = gains(*v[:11]) # A, b, c, Gi, Gx, Gd
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
    A, b, c, *r = gains(*getValues(self.values)[:11])
    sim(self.values["Frame Length"].get(), self.values["End"].get(), A, b, c, self.a)
    self.canvas.show()
  def onControl3D(self):
    self.startgl()
    
  def onController(self, export = False):
    self.a.cla();
    control(self.values["Frame Length"].get(), int(sp.floor(self.values["N"].get())),  self.values["End"].get(), self.a, *gains(*getValues(self.values)[:11]))
    self.canvas.show()
  def onControllerDef(self):
    self.setValues(*controllerDefault)
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