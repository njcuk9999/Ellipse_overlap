# Ellipse_overlap
Using a primary ellipse, find all other ellipses that overlap with primary ellipse

### function ```ellipse_overlap(center1, center2, radius1, radius2, angle1=0.0, angle2=0.0):```

Calculates whether secondary ellipse overlaps with a primary ellipse

* __:param center1:__ numpy array of floats, center x and y of primary ellipse

* __:param center2:__ numpy array of floats, center x and y of secondary ellipse

* __:param radius1:__ numpy array of floats, radius in x and y of primary ellipse

* __:param radius2:__ numpy array of floats, radius in x and y of secondary 
                ellipse
                
* __:param angle1:__ float, angle of rotation, for now has to be zero

* __:param angle2:__ float, angle of rotation, for now has to be zero

* __:return cond:__ boolean, True if there is overlap between primary ellipse and secondary ellipse


### Example

```python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Ellipse_overlap import OverlappingEllipses as ellipse_overlap_multi

# define a set of ellipses randomly
center_1, radius_1 = (2.0, 5.0), (2.0, 3.0)
num = 200

centers = np.random.rand(num, 2) * 10
radii = list(zip(0.5 * np.random.rand(num), 0.5 * np.random.rand(num)))

mask, distances = ellipse_overlap_multi(center_1, centers, radius_1, radii)

# plot to test
plt.close()
fig, frame = plt.subplots(ncols=1, nrows=1)

for e in range(num):
    ellipse = patches.Ellipse(xy=centers[e], width=2 * radii[e][0],
                              height=2 * radii[e][1], angle=0.0)
    frame.add_artist(ellipse)
    ellipse.set_clip_box(frame.bbox)
    if mask[e]:
        ellipse.set_edgecolor('r')
        ellipse.set_facecolor('none')
    else:
        ellipse.set_edgecolor('k')
        ellipse.set_facecolor('none')
    pass

# plt.scatter(centers[mask][:,0], centers[mask][:,1], color='r')
# plt.scatter(centers[~mask][:,0], centers[~mask][:,1], color='k')

ellipse1 = patches.Ellipse(xy=center_1, width=2 * radius_1[0],
                           height=2 * radius_1[1], angle=0.0)
frame.add_artist(ellipse1)
ellipse1.set_clip_box(frame.bbox)
ellipse1.set_edgecolor('b')
ellipse1.set_facecolor('none')

frame.set_title('Num ellipses in red = '
                '{0}'.format(np.sum(np.array(distances) <= 1.0)))

frame.set(xlim=(0, 10), ylim=(0, 10))

plt.show()
plt.close()

```

