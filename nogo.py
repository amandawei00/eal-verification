import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def solp(b):
   return -b*b/(-b*b+2*b-1)

print(solp(-1/3))
nx = np.linspace(-1.,1.,200)
psol = solp(nx)

fig, ax = plt.subplots()
ax.plot(nx, psol, linewidth=2.0)
ax.set(xlim=(-1., 1.))
ax.set(ylim=(0,1))

plt.show()
'''
n_min = 0
n_max = 1
n_init = 0.5

x = np.linspace(0, 1,500)
fig = plt.figure(figsize=(8,3))

x_ax = plt.axes([0., np.pi/2, np.pi, 3*np.pi/2])
slider_ax = plt.axes([0.1, 0.3, 0.5, 0.9])

plt.axes(x_ax)
c_plot, = plt.plot(x, cfunc(n_init, x), 'b')
plt.xlim(0, 1.0)
plt.ylim(0., 0.5)

n_slider=Slider(slider_ax,'a',n_min,n_max,valinit=n_init)

def update(a):
    c_plot.set_ydata(cfunc(x, a))
    fig.canvas.draw_idle()

n_slider.on_changed(update)

plt.show()
'''