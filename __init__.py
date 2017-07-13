#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/05/17 at 4:43 PM

@author: neil

Program description here

Version 0.0.0
"""

from . import ellipse_intersection

__author__ = "Neil Cook"
__email__ = 'neil.james.cook@gmail.com'
__version__ = '0.1'
__all__ = ['Add_buttons', 'Rectangle_Selector']

# =============================================================================
#  Functions
# =============================================================================
OverlappingEllipse = ellipse_intersection.ellipse_overlap
OverlappingEllipses = ellipse_intersection.ellipse_overlap_multi

