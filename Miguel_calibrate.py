import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

# Initial parameters
G12 = 4400.
tau_0 = 23.
p = -0.6
g_f = 0.6

p_0 = 60.
mu = 0.4

# Precompression
s22 = -0.


def solve_curve(g_f, G12, s22, p, p_0, mu, tau_0):
    g_0 = tau_0 / G12
    g12 = np.linspace(0., g_f, 2000)
    s12 = []
    dmg = []
    for g_ in g12:
        if g_ < g_0:
            s12.append(g_ * G12)
            dmg.append(0)
        else:
            # Compute damage
            d = (g_ ** p - g_0 ** p) / (g_f ** p - g_0 ** p)
            # Friction
            s12.append((1. - d) * G12 * g_ - d * np.sign(g_) * mu * (s22 - p_0))
            dmg.append(d)

    return g12, s12, dmg


fig, axes = plt.subplots(ncols=2, figsize=(8, 6), sharex=True)
ax, ax2 = axes
plt.subplots_adjust(left=0.10, bottom=0.35, top=0.95)
g12, s12, dmg = solve_curve(g_f, G12, s22, p, p_0, mu, tau_0)
# strain12 = np.tan(g12)
l_, = ax.plot(g12, s12, 'b-', lw=2)
ax.plot(g12, s12, '-', color='grey', lw=0.5, label='T700/E445 (Costa, 2021)')
ax.set_xlabel(r'$\gamma \, (\%)$')
ax.set_ylabel(r'$\tau \, (MPa)$')
ax.set_xlim(0, 0.30)
ax.set_ylim(0, 140)

# Damage axis
# ax2 = ax.twinx()
ax2.set_ylim(0, 1)
l_2, = ax2.plot(g12, dmg, 'b--', lw=2)
# plt.sca(ax)

axcolor = 'lightgoldenrodyellow'
ax_S22 = plt.axes([0.35, 0.02, 0.55, 0.03], facecolor='white')
ax_p0 = plt.axes([0.35, 0.07, 0.55, 0.03], facecolor=axcolor)
ax_p = plt.axes([0.35, 0.12, 0.55, 0.03], facecolor=axcolor)
ax_tau0 = plt.axes([0.35, 0.17, 0.55, 0.03], facecolor=axcolor)
ax_gf = plt.axes([0.35, 0.22, 0.55, 0.03], facecolor=axcolor)

s_S22 = Slider(ax_S22, r'$\sigma_{22} \, (MPa)$', -100., 100.0, valinit=s22, valstep=1.0)
s_p0 = Slider(ax_p0, r'$p_0 \, (MPa)$', 0.0, 200.0, valinit=p_0, valstep=1.0)
s_p = Slider(ax_p, r'$p$', -1.0, 1.0, valinit=p)
s_tau0 = Slider(ax_tau0, r'$\tau_0 \, (MPa)$', 1, 150, valinit=tau_0, valstep=1.0)
s_gf = Slider(ax_gf, r'$\gamma_f$', 0, 2.0, valinit=g_f, valstep=0.05)


def update(val):
    # print(val)  # Value of the slider changed
    gf_i = s_gf.val
    tau0_i = s_tau0.val
    p_i = s_p.val
    p0_i = s_p0.val
    s22_i = s_S22.val

    mu_i = float(r_mu.value_selected)

    g12, s12, dmg = solve_curve(gf_i, G12, s22_i, p_i, p0_i, mu_i, tau0_i)
    strain12 = np.tan(g12)
    l_.set_xdata(g12)
    l_.set_ydata(s12)
    l_2.set_xdata(g12)
    l_2.set_ydata(dmg)


s_S22.on_changed(update)
s_p0.on_changed(update)
s_p.on_changed(update)
s_tau0.on_changed(update)
s_gf.on_changed(update)

ax_mu = plt.axes([0.10, 0.02, 0.12, 0.23], facecolor=axcolor)
r_mu = RadioButtons(ax_mu, ('0.6', '0.4', '0.2', '0.1', '0'), active=1, activecolor='grey')

colors = {'0.6': 'red', '0.4': 'black', '0.2': 'blue', '0.1': 'yellow', '0': 'green'}


def colorfunc(label):
    l_.set_color(colors[label])
    l_2.set_color(colors[label])
    # r_mu.activecolor = colors[label]
    update(0)
    fig.canvas.draw_idle()


r_mu.on_clicked(colorfunc)

colorfunc(r_mu.value_selected)

# Ramberg-Osgood
alpha = 2.86e-11
eta = 6.49
G = 4830.
tau_ = np.linspace(0, 130.)
gamma_ = (tau_ + alpha * tau_ ** eta) / G
ax.plot(gamma_, tau_, '-', color='cyan', label='AS4/8552 (Bergan, 2018)')

leg = ax.legend()
#leg._draggable(True)

plt.show()
