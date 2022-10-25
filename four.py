from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


xmin = 0
xmax = 1
 
def func(x, a0, a1, a2, a3, a4):
    f = a0 + \
        a1*np.sin( 1.*np.pi*(x+1.)/2. ) + \
        a2*np.sin( 2.*np.pi*(x+1.)/2. ) + \
        a3*np.sin( 3.*np.pi*(x+1.)/2. ) + \
        a4*np.sin( 4.*np.pi*(x+1.)/2. )
    return f  

def cdf(x, a0, a1, a2, a3, a4):

    f = (x+1)/2 + \
        a1*np.sin( 1.*np.pi*(x+1.)/2. ) + \
        a2*np.sin( 2.*np.pi*(x+1.)/2. ) + \
        a3*np.sin( 3.*np.pi*(x+1.)/2. ) + \
        a4*np.sin( 4.*np.pi*(x+1.)/2. )
    return f  
 
def dfdx(x, a0, a1, a2, a3, a4):
    f = a0
    for a, alpha in zip([a1, a2, a3, a4], [1, 2, 3, 4]):
        s = (a*alpha*np.pi/2) * np.cos(alpha*np.pi*(x+1)/2)
        f = f + s
    return f  


t = np.linspace(xmin, xmax, 100)
theta = np.arccos(t)*180./np.pi

a0 = 0.5
a1 = 0
a2 = 0
a3 = 0
a4 = 0

axis_color = 'seashell'
bar_color = 'powderblue'
sty = {'fontsize': 10, 'fontweight': 'bold', 'color': 'darkslategrey'}
textsets = {'rotation':'vertical', 'fontstyle':'italic',
            'fontweight':800, 'fontfamily':'serif',
            'backgroundcolor':(1,1,1,0.5),
            'fontsize':8,
            'color':'rebeccapurple',
            'horizontalalignment':'center'}

fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(left=0.09, bottom=0.2, right=0.98, top=0.95)
 
# PANEL - - - -
ax = fig.add_subplot(221)
 
[line01] = ax.plot(theta, dfdx(t, a0, a1, a2, a3, a4), linewidth=3,
                  color='salmon')
ax.plot([0, 180],[a0, a0], linestyle='--', linewidth=1, color='slategrey')

ax.set_xlim([0, 180])
ax.set_ylim([-0.2, 1.2])
ax.get_xaxis().get_major_formatter().labelOnlyBase = False
ax.grid(color='whitesmoke', linestyle='-', linewidth=0.5)
ax.set_title('Angle distribution function  ', y=0.9, loc='right',
        fontdict=sty)
ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'f($\theta$)')
ax.text(170, 0., 'aligned', **textsets)
ax.text(10, 0., 'antialigned', **textsets)
ax.text(90, 0., 'perpedicular', **textsets)
                                                    
 
# PANEL - - - -
ax = fig.add_subplot(222)
 
[line02] = ax.plot(t, dfdx(t, a0, a1, a2, a3, a4), linewidth=3,
                  color='salmon')
ax.plot([-1, 1],[a0, a0], linestyle='--', linewidth=1, color='slategrey')

ax.set_xlim([xmin, xmax])
ax.set_ylim([-0.2, 1.2])
ax.get_xaxis().get_major_formatter().labelOnlyBase = False
ax.grid(color='whitesmoke', linestyle='-', linewidth=0.5)
ax.set_title(r'cos($\theta$) distribution function  ', y=0.9, loc='right',
        fontdict=sty)
ax.set_xlabel(r'cos($\theta$)')
ax.set_ylabel(r'f(cos($\theta$))')
ax.text(0.75, 0., 'aligned', **textsets)
ax.text(-0.75, 0., 'antialigned', **textsets)
ax.text(0., 0., 'perpedicular', **textsets)
 

# PANEL - - - -
ax = fig.add_subplot(223)
 
[line03] = ax.plot(t, cdf(t, a0, a1, a2, a3, a4), linewidth=3,
                  color='salmon')
ax.plot([-1, 1],[0, 1], linestyle='--', linewidth=1, color='slategrey')

ax.set_xlim([xmin, xmax])
ax.set_ylim([-0.2, 1.2])
ax.get_xaxis().get_major_formatter().labelOnlyBase = False
ax.grid(color='whitesmoke', linestyle='-', linewidth=0.5)
ax.set_title('Cumulative distribution function  ', y=0.9, loc='right',
        fontdict=sty)
ax.set_xlabel(r'cos($\theta$)')
ax.set_ylabel(r'F(cos($\theta$))')

 

# PANEL - - - -
ax = fig.add_subplot(224)

[line04] = ax.plot(t, func(t, a0, a1, a2, a3, a4), linewidth=3,
                 color='salmon')
ax.plot([-1, 1],[a0, a0], linestyle='--', linewidth=1, color='slategrey')

ax.set_xlim([xmin, xmax])
ax.set_ylim([0.3, 0.8])
ax.get_xaxis().get_major_formatter().labelOnlyBase = False
ax.grid(color='whitesmoke', linestyle='-', linewidth=0.5)
ax.set_title('residues w.r.t. uniform dist. CDF  ', y=0.9, loc='right',
        fontdict=sty)
ax.set_xlabel(r'cos($\theta$)')
ax.set_ylabel(r'f(cos($\theta$)) - (cos($\theta$)+1)/2')


# Add two sliders for tweaking the parameters

# Define an axes area and draw sliders:
a0_slider_ax  = fig.add_axes([0.1, 0.12, 0.65, 0.01], 
                clip_on=False,
                facecolor=axis_color)
a0_slider = Slider(a0_slider_ax, 'a0', 0., 1,
                   edgecolor='white', facecolor=bar_color, valinit=a0)

a1_slider_ax = fig.add_axes([0.1, 0.10, 0.65, 0.01], facecolor=axis_color)
a1_slider = Slider(a1_slider_ax, 'a1', -0.3, 0.3, 
                   edgecolor='white', facecolor=bar_color, valinit=a1)

a2_slider_ax = fig.add_axes([0.1, 0.08, 0.65, 0.01], facecolor=axis_color)
a2_slider = Slider(a2_slider_ax, 'a2', -0.15, 0.15,
                   edgecolor='white', facecolor=bar_color, valinit=a2)

a3_slider_ax = fig.add_axes([0.1, 0.06, 0.65, 0.01], facecolor=axis_color)
a3_slider = Slider(a3_slider_ax, 'a3', -0.1, 0.1,
                   edgecolor='white', facecolor=bar_color, valinit=a3)

a4_slider_ax = fig.add_axes([0.1, 0.04, 0.65, 0.01], facecolor=axis_color)
a4_slider = Slider(a4_slider_ax, 'a4', -0.05, 0.05,
                   edgecolor='white', facecolor=bar_color, valinit=a4)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line01.set_ydata(dfdx(t, a0_slider.val, a1_slider.val,  
                        a2_slider.val,  a3_slider.val,  a4_slider.val))
    line02.set_ydata(dfdx(t, a0_slider.val, a1_slider.val,  
                        a2_slider.val,  a3_slider.val,  a4_slider.val))
    line03.set_ydata(cdf(t, a0_slider.val, a1_slider.val,  
                        a2_slider.val,  a3_slider.val,  a4_slider.val))
    line04.set_ydata(func(t, a0_slider.val, a1_slider.val,  
                        a2_slider.val,  a3_slider.val,  a4_slider.val))
    fig.canvas.draw_idle()



a0_slider.on_changed(sliders_on_changed)
a1_slider.on_changed(sliders_on_changed)
a2_slider.on_changed(sliders_on_changed)
a3_slider.on_changed(sliders_on_changed)
a4_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
#reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button_ax = fig.add_axes([0.85, 0.08, 0.08, 0.025])
reset_button = Button(reset_button_ax, 'reset', color=axis_color, 
        #hovercolor='0.975'
        hovercolor='gold'
        )
def reset_button_on_clicked(mouse_event):
    a0_slider.reset()
    a1_slider.reset()
    a2_slider.reset()
    a3_slider.reset()
    a4_slider.reset()
    rn = 39

reset_button.on_clicked(reset_button_on_clicked)

#plt.tight_layout()
plt.show() 
