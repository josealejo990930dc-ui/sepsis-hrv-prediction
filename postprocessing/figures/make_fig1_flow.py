"""Fig 1 (flow of patients and windows) as an image, for PLOS ONE (figures must be separate TIFF files).

Same numbers as figures/fig1_flow.tex (package v1.5 sections 5.4 and 9.7). Arithmetic is asserted.
Usage: python make_fig1_flow.py <figures dir>
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = sys.argv[1]
import os
os.makedirs(OUT, exist_ok=True)
assert 6393 - 5532 == 861 and 861 - 51 == 810 and 810 - (136 + 129 + 2) == 543 and 543 - 22 == 521
assert 3100 - 2555 == 545 and 521 - 487 == 34

INK, INK2, FILL = '#1f1f1e', '#5c5b55', '#f1f0ec'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8})
fig, ax = plt.subplots(figsize=(7.2, 8.4))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

def box(x, y, w, h, text, bold_last=True, fill='white', lw=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=0.8',
                                fc=fill, ec=INK, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', color=INK, linespacing=1.35)

def arrow(x0, y0, x1, y1):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9, mutation_scale=9))

MX, MW = 3, 44          # main column
EX, EW = 53, 45         # exclusion column
H = 6.6
ys = [95, 84, 73, 62, 51, 40]
main = ['First ICU stays meeting Sepsis-3 criteria\nMIMIC-III, MetaVision 2008–2012\n6393 patients',
        'Linked to an ECG waveform record\n(start 0–48 h from anchor; shift <12 h)\n861 patients',
        '≥3 h of signal after hour 0\n810 patients',
        'HRV extracted (lead II, 125 Hz)\n543 patients; 18 054 valid HRV hours',
        '≥1 candidate prediction window\n521 patients; 3100 windows',
        'Analytic cohort: windows with ≥6 valid HRV hours\n487 patients; 2555 windows']
excl = ['No usable waveform record: 5532',
        '<3 h of signal: 51',
        '<12 h of signal after hour 0: 136\nNo valid HRV hour: 129\nTechnical failure: 2 (no lead II; server error)',
        'Valid hours only before hour 11 (no window): 22',
        'Windows with <6 valid HRV hours: 545 (34 patients)\nSensitivity: ≥4 h 2734 windows/501 patients;\n≥8 h 2367/465']
for i, (y, t) in enumerate(zip(ys, main)):
    box(MX, y - H / 2, MW, H, t, lw=1.6 if i == 5 else 1.0)
    if i < 5:
        arrow(MX + MW / 2, y - H / 2 - 0.6, MX + MW / 2, ys[i + 1] + H / 2 + 0.6)
        ym = (y + ys[i + 1]) / 2
        eh = 7.2 if i in (2, 4) else 4.2
        box(EX, ym - eh / 2, EW, eh, excl[i], fill=FILL)
        ax.plot([MX + MW / 2, EX - 0.6], [ym, ym], color=INK, lw=0.9)
        ax.annotate('', xy=(EX - 0.5, ym), xytext=(EX - 2.5, ym),
                    arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9, mutation_scale=9))

# discharge exclusion, observable-horizon windows, and outcome-specific evaluability
assert 2555 - 188 == 2367 and 2367 - 686 == 1681 and 2367 - 2050 == 317
Y7 = 24
arrow(MX + MW / 2, ys[5] - H / 2 - 0.6, MX + MW / 2, Y7 + H / 2 + 0.6)
ym = (ys[5] + Y7) / 2
box(EX, ym - 2.1, EW, 4.2, 'Discharged alive during the 6-h horizon:\n188 windows (164 patients)', fill=FILL)
ax.plot([MX + MW / 2, EX - 0.6], [ym, ym], color=INK, lw=0.9)
ax.annotate('', xy=(EX - 0.5, ym), xytext=(EX - 2.5, ym),
            arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9, mutation_scale=9))
box(MX, Y7 - H / 2, MW, H, 'Windows with an observable horizon\n468 patients; 2367 windows\n'
    '(9 ended in death: classified as negative,\ndeterioration)')
# split to the two outcome boxes
yb = Y7 - H / 2 - 0.6
ax.plot([MX + MW / 2, MX + MW / 2], [yb, yb - 2.2], color=INK, lw=0.9)
ax.plot([MX + 22, EX + 22.5], [yb - 2.2, yb - 2.2], color=INK, lw=0.9)
arrow(MX + 22, yb - 2.2, MX + 22, 11.6)
arrow(EX + 22.5, yb - 2.2, EX + 22.5, 11.6)
box(MX, 1.0, 44, 10.0, 'SOFA improvement and SOFA deterioration\n1681 windows; 405 patients\n'
    'Improved 23.6%; deterioration 94 (5.6%)\nNot evaluable: 686 windows\n(paired SOFA components missing)')
box(EX, 1.0, 45, 10.0, 'Haemodynamic stability or recovery\n317 windows; 111 patients\nStable or recovered 45.4%\n'
    'Not evaluable: 2050 windows\n(MAP or lactate missing in a window)')
fig.savefig(f'{OUT}/fig1_flow.png', dpi=300, bbox_inches='tight', pad_inches=0.03)
print('written', f'{OUT}/fig1_flow.png')
