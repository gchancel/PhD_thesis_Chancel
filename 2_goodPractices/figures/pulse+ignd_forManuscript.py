# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 20:00:21 2023.

@author: geoff2
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter

plt.style.use('ggplot')

form_volt = EngFormatter(unit='V')
form_amp = EngFormatter(unit='A')
form_sec = EngFormatter(unit='s')

basePath = "C:/Users/geoff2/Seafile/Thèse/PUBLI FDTC 2023"
pathBad = f"{basePath}/badGND/"
# pathGood = f"{basePath}/goodGND_ONLY/"
pathGoodImp = f"{basePath}/goodGND+50Ohms/"

badGND = {}
# goodGND = {}
impGND = {}

badGND['iTime'] = np.load(pathBad + "iGNDTime.npy", allow_pickle = True)
badGND['vTime'] = np.load(pathBad + "vPulseTime.npy", allow_pickle = True)
badGND['iGND'] = np.load(pathBad + "iGNDData.npy", allow_pickle = True) / 5.0
badGND['vPulse'] = np.load(pathBad + "vPulseData.npy", allow_pickle = True)

# goodGND['iTime'] = np.load(pathGood + "iGNDTime.npy", allow_pickle = True)
# goodGND['vTime'] = np.load(pathGood + "vPulseTime.npy", allow_pickle = True)
# goodGND['iGND'] = np.load(pathGood + "iGNDData.npy", allow_pickle = True) / 5.0
# goodGND['vPulse'] = np.load(pathGood + "vPulseData.npy", allow_pickle = True)

impGND['iTime'] = np.load(pathGoodImp + "iGNDTime.npy", allow_pickle = True)
impGND['vTime'] = np.load(pathGoodImp + "vPulseTime.npy", allow_pickle = True)

impGND['iGND'] = np.load(pathGoodImp + "iGNDData.npy", allow_pickle = True) / 5.0
impGND['vPulse'] = np.load(pathGoodImp + "vPulseData.npy", allow_pickle = True)

###############################################################################################
# # Reducing time axis for more space
# # beginImpGndiTime = np.where(np.isclose(impGND['iTime'], 30e-9, atol = 100e-12))[0][0]
# beginImpGndvTime = np.where(np.isclose(impGND['vTime'], 1e-9, atol = 100e-12))[0][0]
# impGND['iTime'] = impGND['iTime'][beginImpGndvTime:]
# # impGND['iTime'] = impGND['iTime'] - impGND['iTime'].min()

# impGND['vTime'] = impGND['vTime'][beginImpGndvTime:]
# # impGND['vTime'] = impGND['vTime'] - impGND['vTime'].min()

# impGND['iGND'] = impGND['iGND'][beginImpGndvTime:]
# impGND['vPulse'] = impGND['vPulse'][beginImpGndvTime:]
###############################################################################################

idealPulse = np.zeros(impGND['vTime'].shape)
borneGauche = np.where(impGND['vTime'] > 14e-9)[0]
borneDroite = np.where(impGND['vTime'] < 36e-9)[0]
# bornesAll = np.where(borneGauche == borneDroite)[0]
bornesAll = np.intersect1d(borneGauche, borneDroite)
idealPulse[bornesAll] = -140

beginBad = np.where(np.isclose(badGND['iTime'], impGND['iTime'][0], atol = 100e-12))[0][0]
endBad = np.where(np.isclose(badGND['iTime'], impGND['iTime'][-1], atol = 100e-12))[0][-1]

badGND['iTime'] = badGND['iTime'][beginBad:endBad]
badGND['vTime'] = badGND['vTime'][beginBad:endBad]

badGND['iTime'][:] -= badGND['iTime'].min()
badGND['vTime'][:] -= badGND['vTime'].min()

impGND['iTime'][:] -= impGND['iTime'].min()
impGND['vTime'][:] -= impGND['vTime'].min()

rmsBad = np.sqrt(np.mean(np.square(badGND['iGND'][beginBad:endBad])))
# rmsGood = np.sqrt(np.mean(np.square(goodGND['iGND'])))
rmsImp = np.sqrt(np.mean(np.square(impGND['iGND'])))

# Original scale
# scale = 1.0
# lw = 2
# labSize = 10
# fSize = 12 * 0.9
# ySize = 13 * scale
# titleSize = 14 * scale
# tPad = 15 * scale
# tLoc0 = 'center'
# tLoc1 = 'center'
# tLoc2 = 'center'
# yLabelPad = 15 * scale
# tBox = dict(facecolor = 'yellow', pad = 4 * scale, alpha = 0.2)
# sBox = dict(facecolor = 'blue', pad = 4 * scale, alpha = 0.2)
# xBox = dict(facecolor = 'green', pad = 4 * scale, alpha = 0.2)
# supSize = 20
# yLabelRot = 90
# textposX = 0.99
# textposY = 0.01
# offset = 0.18 / 2

# New scale
scale = 1.0
lw = 2
labSize = 10
fSize = 12 * 1.05
ySize = 13 * scale
titleSize = 14 * scale
tPad = 15 * scale
tLoc0 = 'center'
tLoc1 = 'center'
tLoc2 = 'center'
yLabelPad = 15 * scale
tBox = dict(facecolor = 'yellow', pad = 4 * scale, alpha = 0.2)
sBox = dict(facecolor = 'blue', pad = 4 * scale, alpha = 0.2)
xBox = dict(facecolor = 'green', pad = 4 * scale, alpha = 0.2)
supSize = 20
yLabelRot = 90
textposX = 0.99
textposY = 0.01
offset = 0.25 / 2

# subAdj = {'top': 0.915,
#           'bottom': 0.031,
#           'left': 0.155,
#           'right': 0.863,
#           'hspace': 0.745,
#           'wspace': 0.200}

# fig, ax = plt.subplots(4, figsize = (5.0 * scale, 11 * 0.67 * scale), sharex = True)
fig, ax = plt.subplots(2, 2, figsize = (11.5 * scale, 5 * scale), sharex = True)
# fig.suptitle("BBI practices comparison", size = supSize * scale, bbox = sBox, y = 0.99)

ax[0, 1].plot(badGND['iTime'], badGND['iGND'][beginBad:endBad],
              lw = lw * scale)
ax[0, 1].yaxis.set_major_formatter(form_amp)
ax[0, 1].xaxis.set_major_formatter(form_sec)
ax[0, 1].text(textposX,
              textposY,
              r'RMS $I_{GND}$ = ' + f'{rmsBad * 1000:.2f} mA',
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[0, 1].transAxes,
              fontsize = fSize * scale)
ax[0, 1].set_ylabel('$I_{GND}(t)$', fontsize = ySize, labelpad = yLabelPad, bbox = tBox,
                    rotation = yLabelRot)
# ax[0, 1].set_title("State of the art grounding",
#                    fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc0)
ax[0, 1].set_title("Default BBI platform (IC current)",
                   fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc0)

ax[0, 0].plot(badGND['vTime'], badGND['vPulse'][beginBad:endBad],
              lw = lw * scale)
ax[0, 0].yaxis.set_major_formatter(form_volt)
ax[0, 0].xaxis.set_major_formatter(form_sec)
# ax[0, 0].set_title("State of the art grounding",
#                    fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc0)
ax[0, 0].set_title("Default BBI platform (Generator voltage)",
                   fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc0)
ax[0, 0].text(textposX,
              textposY,
              f"Max. abs. voltage = {np.amin(badGND['vPulse'][beginBad:endBad]):.1f} V",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[0, 0].transAxes,
              fontsize = fSize * scale)
ax[0, 0].text(textposX,
              textposY + offset,
              "Set-point voltage = -140 V",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[0, 0].transAxes,
              fontsize = fSize * scale)
ax[0, 0].text(textposX,
              textposY + offset * 2.0,
              r"Rise time $\approx$ 65 ns",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[0, 0].transAxes,
              fontsize = fSize * scale)
ax[0, 0].text(textposX,
              textposY + offset * 3.0,
              r"Fall time $\approx$ 16 ns",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[0, 0].transAxes,
              fontsize = fSize * scale)
ax[0, 0].text(textposX,
              textposY + offset * 4.0,
              r"Pulse width $\approx$ 75 ns",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[0, 0].transAxes,
              fontsize = fSize * scale)
ax[0, 0].set_ylabel('$V_{PULSE}(t)$', fontsize = ySize, labelpad = yLabelPad, bbox = tBox,
                    rotation = yLabelRot)
# ax[0].get_shared_x_axes().join(ax[0], ax[1])
# ax[1].get_shared_x_axes().join(ax[1], ax[0])

# ax[0].vlines(x = 12e-9,
#              ymin = badGND['vPulse'][beginBad:endBad].min() * 1.0,
#              ymax = badGND['vPulse'][beginBad:endBad].max() * 1.0,
#              color = 'red',
#              lw = 1)

# ax[3].plot(goodGND['iTime'], goodGND['iGND'],
#            lw = lw * scale)
# ax[3].yaxis.set_major_formatter(form_amp)
# ax[3].xaxis.set_major_formatter(form_sec)
# ax[3].text(textposX,
#            textposY,
#            r'RMS $I_{GND}$ = ' + f'{round(rmsGood * 1000, 2)} mA',
#            verticalalignment = 'bottom', horizontalalignment = 'right',
#            transform = ax[3].transAxes,
#            fontsize = fSize * scale)
# ax[3].set_ylabel('$I_{GND}(t)$', fontsize = ySize, labelpad = yLabelPad, bbox = tBox,
#                  rotation = yLabelRot)
# ax[3].vlines(x = 20e-9,
#              ymin = goodGND['iGND'].min() * 1.0,
#              ymax = goodGND['iGND'].max() * 1.0,
#              color = 'red',
#              lw = 1)

# ax[2].plot(goodGND['vTime'], goodGND['vPulse'],
#            lw = lw * scale)
# ax[2].yaxis.set_major_formatter(form_volt)
# ax[2].xaxis.set_major_formatter(form_sec)
# ax[2].set_title("Scenario 2: Ground  bypass",
#                 fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc1)
# ax[2].text(textposX,
#            textposY,
#            f"Max. abs. voltage = {np.amin(goodGND['vPulse']):.1f} V",
#            verticalalignment = 'bottom', horizontalalignment = 'right',
#            transform = ax[2].transAxes,
#            fontsize = fSize * scale)
# ax[2].text(textposX,
#            textposY + offset,
#            "Set-point voltage = -140 V",
#            verticalalignment = 'bottom', horizontalalignment = 'right',
#            transform = ax[2].transAxes,
#            fontsize = fSize * scale)
# ax[2].text(textposX,
#            textposY + offset * 2.0,
#            r"Rise time $\approx$ 65 ns",
#            verticalalignment = 'bottom', horizontalalignment = 'right',
#            transform = ax[2].transAxes,
#            fontsize = fSize * scale)
# ax[2].text(textposX,
#            textposY + offset * 3.0,
#            r"Fall time $\approx$ 16 ns",
#            verticalalignment = 'bottom', horizontalalignment = 'right',
#            transform = ax[2].transAxes,
#            fontsize = fSize * scale)
# ax[2].set_ylabel('$V_{PULSE}(t)$', fontsize = ySize, labelpad = yLabelPad, bbox = tBox,
#                  rotation = yLabelRot)
# ax[2].vlines(x = 12e-9,
#              ymin = goodGND['vPulse'].min() * 1.0,
#              ymax = goodGND['vPulse'].max() * 1.0,
#              color = 'red',
#              lw = 1)

ax[1, 1].plot(impGND['iTime'], impGND['iGND'],
              lw = lw * scale)
ax[1, 1].yaxis.set_major_formatter(form_amp)
ax[1, 1].xaxis.set_major_formatter(form_sec)
ax[1, 1].text(textposX,
              textposY,
              r'RMS $I_{GND}$ = ' + f'{round(rmsImp * 1000, 2)} mA',
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[1, 1].transAxes,
              fontsize = fSize * scale)
ax[1, 1].set_ylabel('$I_{GND}(t)$', fontsize = ySize, labelpad = yLabelPad, bbox = tBox,
                    rotation = yLabelRot)
# ax[1, 1].set_title("Ground bypass and impedance matching",
#                    fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc2)
ax[1, 1].set_title("Enhanced BBI platform (IC current)",
                   fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc2)

ax[1, 0].plot(impGND['vTime'], impGND['vPulse'],
              lw = lw * scale)
ax[1, 0].yaxis.set_major_formatter(form_volt)
ax[1, 0].xaxis.set_major_formatter(form_sec)
# ax[1, 0].set_title("Ground bypass and impedance matching",
#                    fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc2)
ax[1, 0].set_title("Enhanced BBI platform (Generator voltage)",
                   fontsize = titleSize, bbox = tBox, pad = tPad, loc = tLoc2)
ax[1, 0].text(textposX,
              textposY,
              f"Max. abs. voltage = {np.amin(impGND['vPulse']):.1f} V",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[1, 0].transAxes,
              fontsize = fSize * scale)
ax[1, 0].text(textposX,
              textposY + offset,
              "Set-point voltage = -140 V",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[1, 0].transAxes,
              fontsize = fSize * scale)
ax[1, 0].text(textposX,
              textposY + offset * 2.0,
              r"Rise and fall time $\approx$ 16 ns",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[1, 0].transAxes,
              fontsize = fSize * scale)
ax[1, 0].text(textposX,
              textposY + offset * 3.0,
              r"Pulse width $\approx$ 20 ns",
              verticalalignment = 'bottom', horizontalalignment = 'right',
              transform = ax[1, 0].transAxes,
              fontsize = fSize * scale)
ax[1, 0].set_ylabel('$V_{PULSE}(t)$', fontsize = ySize, labelpad = yLabelPad, bbox = tBox,
                    rotation = yLabelRot)
# ax[2].vlines(x = 12e-9,
#              ymin = impGND['vPulse'].min() * 1.0,
#              ymax = impGND['vPulse'].max() * 1.0,
#              color = 'red',
#              lw = 1)

counter = 0
for axf in ax.ravel():
    axf.set_xlabel("Time", fontsize = ySize, bbox = xBox)
    axf.tick_params(axis = 'both', which = 'both', reset = True,
                    labelsize = labSize * scale)
    axf.yaxis.set_label_position("right")
    axf.text(0.05,
             0.05,
             f"S{counter//2 + 1}{'P' if counter % 2 == 0 else 'G'}",
             verticalalignment = 'bottom', horizontalalignment = 'left',
             transform = axf.transAxes,
             fontsize = fSize * 1.5 * scale,
             color = 'red')
    counter += 1

fig.tight_layout()
fig.show()

try:
    fig.subplots_adjust(top = subAdj['top'],
                        bottom = subAdj['bottom'],
                        left = subAdj['left'],
                        right = subAdj['right'],
                        hspace = subAdj['hspace'],
                        wspace = subAdj['wspace'])
except NameError:
    pass

saveDir = "C:/Users/geoff2/Dropbox/Applications/Overleaf/\
PhD_thesis_Chancel/2_goodPractices/figures/"

# plt.savefig(f"{saveDir}realPulsesComparisons.svg",
#             transparent = True)

plt.savefig(f"{saveDir}realPulsesComparisons.pdf",
            transparent = False)

fig.canvas.manager.window.raise_()
