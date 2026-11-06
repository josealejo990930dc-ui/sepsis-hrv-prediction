"""Convert main-text figures to PLOS ONE TIFF (RGB, LZW, 300 dpi, no alpha) and check size limits."""
import sys
from PIL import Image
src, dst = sys.argv[1], sys.argv[2]
for png, tif in [('fig1_flow.png', 'Fig1.tif'), ('fig2_forest.png', 'Fig2.tif')]:
    im = Image.open(f'{src}/{png}').convert('RGB')
    w, h = im.size
    ok = 789 <= w <= 2250 and h <= 2625
    im.save(f'{dst}/{tif}', compression='tiff_lzw', dpi=(300, 300))
    print(tif, w, 'x', h, 'px', 'OK' if ok else 'OUT OF PLOS LIMITS')
