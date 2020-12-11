# -*- coding: utf-8 -*-
"""
Script to visualise the lamination parameters
of a population of stacking sequences in 4 subplots
    - (lampam_1, lampam_2)
    - (lampam_1, lampam_3)
    - (lampam_9, lampam_10)
    - (lampam_9, lampam_11)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Set of constraints:
Const = 'C0'
#Const = 'C1'

# Number of plies
n_plies = 40
#n_plies = 80
#n_plies = 200

data_filename = 'pop_sym_' + Const + '_' + str(n_plies) + 'plies'

data = pd.read_excel(data_filename + '.xlsx', 'stacks')

print(len(data))

## Displaying the results on graphs

lampam_1 = list(data['lampam[1]'])
lampam_2 = list(data['lampam[2]'])
lampam_3 = list(data['lampam[3]'])
lampam_9 = list(data['lampam[9]'])
lampam_10 = list(data['lampam[10]'])
lampam_11 = list(data['lampam[11]'])

fig1, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
axes = axes.flatten()
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
ax4 = axes[3]

#my_labelsize = 15
#my_titlesize = 20
#my_axissize = 18

my_labelsize = 30
my_titlesize = 30
my_axissize = 40
my_font = 'Times new Roman'

for ax in axes:
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontweight('normal') for label in labels]
    [label.set_fontname(my_font) for label in labels]


# FIRST SUBPLOT

# boundaries
x = np.linspace(-1, 1, 1000)
y1 = np.ones(len(x))
y2 = 2*x**2-1
x3 = np.linspace(-0.6, 0.6, 1000)
y3 = 0.6*np.ones(len(x3))
x4 = np.linspace(-0.6, 0, 1000)
x5 = np.linspace(0, 0.6, 1000)
y4 = -2*x4 -0.6
y5 = 2*x5 -0.6

# plot
ax1.plot(lampam_1, lampam_2,'bx',linewidth = 2, markersize=8)
ax1.plot(x,y1,'k-',linewidth = 4)
ax1.plot(x,y2,'k-', linewidth = 2)

if Const == 'C1':
    ax1.plot(x3, y3,'k-', linewidth = 2)
    ax1.plot(x4, y4,'k-', linewidth = 2)
    ax1.plot(x5, y5,'k-', linewidth = 2)

# formatting
#ax1.set_title('(a)',size = my_titlesize, fontname = my_font, fontweight="normal") # Title # fontweight="bold"
#ax1.set_xlabel(r'$\mathcal{\xi}_1$', fontsize = my_axissize, fontweight="normal") # X label
#ax1.set_ylabel(r'$\mathcal{\xi}\mathcal{2}$', fontweight = "normal", fontsize = my_axissize, fontname=my_font) # Y label
ax1.set_xlim([-1, 1])
ax1.set_ylim([-1, 1])
ax1.grid(True)
#ax1.xaxis.set_ticks(np.linspace(-1, 1, 5, endpoint=True))
#ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#ax1.yaxis.set_ticks(np.linspace(-1, 1, 5, endpoint=True))
#ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#ax1.tick_params(direction='out', length = 6, width = 2, colors = 'black', labelsize = my_labelsize)
ax1.xaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax1.yaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax1.tick_params(direction='out', length = 6, width = 1, colors = 'black', labelsize = my_labelsize)
#ax1.legend(handletextpad=0,loc='upper right') # Plot legend

# SECOND SUBPLOT

# boundaries
th = np.linspace(0, 2*np.pi, 10_000)
xunit = np.cos(th)
yunit = np.sin(th)

# plot
ax2.plot(lampam_1, lampam_3,'bx',linewidth = 1, markersize=8)
ax2.plot(xunit, yunit,'k-', linewidth = 2)

# formatting
#ax2.set_title('(b)', size= my_titlesize, fontname=my_font) # Title
#ax2.set_ylabel(r'$\mathcal{\xi}_\mathcal{3}$', fontweight="normal", fontsize = my_axissize, fontname=my_font) # Y label
#ax2.set_xlabel(r'$\mathcal{\xi}_\mathcal{1}$', fontsize =my_axissize , fontname=my_font) # X label
ax2.set_xlim([-1, 1])
ax2.set_ylim([-1, 1])
ax2.grid(True)
ax2.xaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax2.yaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax2.tick_params(direction='out', length = 6, width = 1, colors = 'black', labelsize = my_labelsize)
#ax2.legend(handletextpad=0,loc='upper right') # Plot legend

# THIRD SUBPLOT

# plot
ax3.plot(lampam_9, lampam_10,'bx',linewidth = 1, markersize=8)
ax3.plot(x,y1,'k-',linewidth = 4)
ax3.plot(x,y2,'k-', linewidth = 2)

# formatting
#ax3.set_title('(c)',size= my_titlesize, fontname=my_font) # Title # fontweight="bold"
#ax3.set_ylabel(r'$\mathcal{\xi}_\mathcal{10}$', fontweight="normal", fontsize = my_axissize, fontname=my_font) # Y label
#ax3.set_xlabel(r'$\mathcal{\xi}_\mathcal{9}$', fontsize =my_axissize , fontname=my_font) # X label
ax3.set_xlim([-1, 1])
ax3.set_ylim([-1, 1])
ax3.grid(True)
ax3.xaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax3.yaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax3.tick_params(direction='out', length = 6, width = 1, colors = 'black', labelsize = my_labelsize)
#ax3.legend(handletextpad=0,loc='upper right') # Plot legend

# FOURTH SUBPLOT

# plot
ax4.plot(lampam_9, lampam_11,'bx',linewidth = 1, markersize=8)
ax4.plot(xunit, yunit,'k-', linewidth = 2)

# formatting
#ax4.set_title('(d)', size= my_titlesize, fontname=my_font) # Title
#ax4.set_ylabel(r'$\mathcal{\xi}_\mathcal{11}$', fontweight="normal", fontsize = my_axissize, fontname=my_font) # Y label
#ax4.set_xlabel(r'$\mathcal{\xi}_\mathcal{9}$', fontsize =my_axissize , fontname=my_font) # X label
ax4.set_xlim([-1, 1])
ax4.set_ylim([-1, 1])
ax4.grid(True)
ax4.xaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax4.yaxis.set_ticks(np.linspace(-1, 1, 3, endpoint=True))
ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
ax4.tick_params(direction='out', length = 6, width = 1, colors = 'black', labelsize = my_labelsize)
#ax4.legend(handletextpad=0,loc='upper right') # Plot legend

plt.tight_layout(pad=2, w_pad=10, h_pad=10)
plt.savefig(data_filename + '.svg')

