import matplotlib.colors as mc
import numpy as np
# lighter blue -> light amethyst -> lighter red
cmap=mc.LinearSegmentedColormap.from_list('bluered',
      [(0.00,'#5A86D6'),(0.50,'#9E6CC0'),(1.00,'#E86A6A')])
order=['cornerT','innerT','quad','pent','hex']
samples=np.linspace(0.05,0.95,5)
COL={k:mc.to_hex(cmap(s)) for k,s in zip(order,samples)}
print("COL =", COL)
