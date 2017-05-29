#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/05/17 at 3:02 PM

@author: neil

Program description here

Version 0.0.0
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# =============================================================================
# Define functions
# =============================================================================
def bisect(f, t_0, t_1, err=0.0001, max_iter=100):
    """
    Bisect fitting function
    
    from: https://stackoverflow.com/a/24710526/7858439
    
    :param f: 
    :param t_0: 
    :param t_1: 
    :param err: 
    :param max_iter: 
    :return: 
    """
    iteration = 0
    ft_0 = f(t_0)
    ft_1 = f(t_1)
    try:
        assert ft_0 * ft_1 <= 0.0
    except AssertionError:
        return np.nan
    while True:
        t = 0.5 * (t_0 + t_1)
        ft = f(t)
        if iteration >= max_iter or ft < err:
            return t
        if ft * ft_0 <= 0.0:
            t_1 = t
            ft_1 = ft
        else:
            t_0 = t
            ft_0 = ft
        iteration += 1


class Ellipse(object):
    """
    
    Ellipse class
    
    from:     from: https://stackoverflow.com/a/24710526/7858439
    
    """

    def __init__(self, center, radius, angle=0.0):
        """
        
        :param center: list/numpy array of floats, shape(2), x and y position
                       of ellipse
                       
        :param radius: list/numpy array of floats, shape(2), x and y radius
                       of ellipse
                       
        :param angle: float, angle with respect to the y axis (only works for
                      0.0
        
        Ellipse equation:                                                             
        (x-center_x)^2/radius_x^2 + (y-center_y)^2/radius_y^2 = 1                     
        x = center_x + radius_x * cos(t)                                              
        y = center_y + radius_y * sin(t)  
        
        """
        assert len(center) == 2
        assert len(radius) == 2
        self.center = np.array(center)
        self.radius = np.array(radius)
        self.angle = angle
        # No rotation included yet TODO: Need rotation
        if angle != 0.0:
            raise NotImplementedError("Ellipses must be aligned with axis")

    def distance_from_origin(self):
        """                                                                           
        Ellipse equation:                                                             
        (x-center_x)^2/radius_x^2 + (y-center_y)^2/radius_y^2 = 1                     
        x = center_x + radius_x * cos(t)                                              
        y = center_y + radius_y * sin(t)                                              
        """
        cc = self.center
        rr = self.radius
        # rotate ellipse of -angle to become axis aligned
        c, s = np.cos(self.angle), np.sin(self.angle)
        cc = (c * cc[0] + s * cc[1], -s * cc[0] + c * cc[1])
        f = lambda t: (rr[1] * (cc[1] + rr[1] * np.sin(t)) * np.cos(t) -
                       rr[0] * (cc[0] + rr[0] * np.cos(t)) * np.sin(t))
        if cc[0] > 0.0:
            if cc[1] > 0.0:
                t_0, t_1 = -np.pi, -np.pi / 2
            else:
                t_0, t_1 = np.pi / 2, np.pi
        else:
            if cc[1] > 0.0:
                t_0, t_1 = -np.pi / 2, 0
            else:
                t_0, t_1 = 0, np.pi / 2
        t = bisect(f, t_0, t_1)
        if not np.isfinite(t):
            return np.nan
        else:
            x = cc[0] + rr[0] * np.cos(t)
            y = cc[1] + rr[1] * np.sin(t)
            return np.sqrt(x ** 2 + y ** 2)


def distance_from_primary(center1, center2, radius1, radius2):
    # get centers and radii
    old_center1, old_center2 = center1, center2
    old_radius1, old_radius2 = radius1, radius2
    # choose ellipse 1 as the primary ellipse
    # ellipse 1 --> ellipse 3, ellipse 2 --> ellipse 4

    # 1. move centers so ellipse 1 is at origin (0, 0)
    new_center1 = old_center1 - old_center1
    new_center2 = old_center2 - old_center1

    # 2. rescale ellipse 1 so it is a circle (of radius 1)
    new_radius1 = old_radius1 / old_radius1
    new_radius2 = old_radius2 / old_radius1
    # need to rescale distance between objects too
    new_center2 = new_center2 / old_radius1

    # 3. TODO: Rotate ellipse 2 so is it inline with the origin

    # 4. create new ellipses for e2
    e2 = Ellipse(new_center2, new_radius2)
    # 5. return distances (a distance of 1 is classed as inside)
    return np.array(e2.distance_from_origin())


def ellipse_overlap(center1, center2, radius1, radius2, angle1=0.0, angle2=0.0):
    """
    Calculates whether secondary ellipse overlaps with a primary ellipse

    :param center1: numpy array of floats, center x and y of primary ellipse

    :param center2: numpy array of floats, center x and y of secondary ellipse

    :param radius1: numpy array of floats, radius in x and y of primary ellipse

    :param radius2: numpy array of floats, radius in x and y of secondary 
                    ellipse
                    
    :param angle1: float, angle of rotation, for now has to be zero
    
    :param angle2: float, angle of rotation, for now has to be zero

    :return cond: boolean, True if there is overlap between primary ellipse
                  and secondary ellipse
    """

    # make sure centers and radii are numpy arrays
    center1, center2 = np.array(center1), np.array(center2)
    radius1, radius2 = np.array(radius1), np.array(radius2)
    assert len(center1) == 2
    assert len(center2) == 2
    assert len(radius1) == 2
    assert len(radius2) == 2

    # No rotation included yet TODO: Need rotation
    if angle1 != 0.0 or angle2 != 0.0:
        raise NotImplementedError("Ellipses must be aligned with axis")

    # find distance from primary
    distance2 = distance_from_primary(center1, center2, radius1, radius2)
    distance1 = distance_from_primary(center2, center1, radius2, radius1)

    # a distance of less than or equal to 1 is classed as inside
    # thus overlap is equivalent to distance 1 or distance 2 being less than
    # 1. We need both distances as primary could be completely contained within
    # secondary
    overlap = (distance2 <= 1.0) | (distance1 <= 1.0)
    # overlap
    return overlap


def ellipse_overlap_multi(center1, centers, radius1, radii):
    """
    Calculates which ellipses overlap with a primary ellipse
    
    :param center1: numpy array of floats, center x and y of primary ellipse
    
    :param centers: numpy array of floats, shape(Nx2), array of centers of
                    ellipses where each row has a center x and y position 
                    for each ellipse (relative to primary ellipse)
                    
    :param radius1: numpy array of floats, radius in x and y of primary ellipse
    
    :param radii: numpy array of floats, shape(Nx2), array of radii of
                  ellipses where each row has a radius in x and y 
                  for each ellipse (relative to primary ellipse)
                  
    :return mask: numpy array of boolean shape(N), True if there is overlap 
                  between primary ellipse and ellipse from centers/radii row
    """
    assert len(center1) == 2
    assert len(radius1) == 2

    mask = np.zeros(len(centers), dtype=bool)
    for row in range(len(centers)):
        center2, radius2 = centers[row], radii[row]
        mask[row] = ellipse_overlap(center1, center2, radius1, radius2)

    return mask


def test_plot(old_center1, old_center2, old_radius1, old_radius2,
              new_center1, new_center2, new_radius1, new_radius2):
    """
    Plot to test transformations
    :param old_center1: 
    :param old_center2: 
    :param old_radius1: 
    :param old_radius2: 
    :param new_center1: 
    :param new_center2: 
    :param new_radius1: 
    :param new_radius2: 
    :return: 
    """

    cs, rs = [old_center1, old_center2], [old_radius1, old_radius2]
    cs1, rs1 = [new_center1, new_center2], [new_radius1, new_radius2]
    mask = [True, False]
    plt.close()
    fig, frame = plt.subplots(ncols=2, nrows=1)
    for e in range(len(cs)):
        ellipse = patches.Ellipse(xy=cs[e], width=2 * rs[e][0],
                                  height=2 * rs[e][1], angle=0.0)
        frame[0].add_artist(ellipse)
        ellipse.set_clip_box(frame[0].bbox)
        if mask[e]:
            ellipse.set_edgecolor('r')
            ellipse.set_facecolor('none')
        else:
            ellipse.set_edgecolor('k')
            ellipse.set_facecolor('none')

    for e in range(len(cs1)):
        ellipse = patches.Ellipse(xy=cs1[e], width=2 * rs1[e][0],
                                  height=2 * rs1[e][1], angle=0.0)
        frame[1].add_artist(ellipse)
        ellipse.set_clip_box(frame[1].bbox)
        if mask[e]:
            ellipse.set_edgecolor('r')
            ellipse.set_facecolor('none')
        else:
            ellipse.set_edgecolor('k')
            ellipse.set_facecolor('none')

    title1 = ('centers 1=({0:.3f}, {1:.3f}) 2=({2:.3f}, {3:.3f})'.format(*(
        list(old_center1) + list(old_center2))))
    title1 += ("\n ")
    title1 += ('radii 1=({0:.3f}, {1:.3f}) 2=({2:.3f}, {3:.3f})'.format(*(
        list(old_radius1) + list(old_radius2))))
    title1 += ("\n ")

    title2 = ('centers 1=({0:.3f}, {1:.3f}) 2=({2:.3f}, {3:.3f})'.format(*(
        list(new_center1) + list(new_center2))))
    title2 += ("\n ")
    title2 += ('radii 1=({0:.3f}, {1:.3f}) 2=({2:.3f}, {3:.3f})'.format(*(
        list(new_radius1) + list(new_radius2))))
    title2 += ("\n ")

    frame[0].set(xlim=(0, 10), ylim=(0, 10), title=title1)
    frame[1].set(xlim=(0, 10), ylim=(0, 10), title=title2)

    plt.show()
    plt.close()


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    num = 200
    # define a set of ellipses randomlw
    center_1, radius_1 = (2.0, 5.0), (2.0, 3.0)

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

# =============================================================================
# End of code
# =============================================================================
