'''
This script contains some functions to add colorbars for various parameters
'''
import matplotlib.pyplot as plt
def addtempcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for CAPE

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0
    bottom = axes_bbox.y0 - 0.015
    width = 0.38
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('2m Temperature (F)', size=8)  # MODIFY THIS for other fields!! MODIFY THIS for other fields!!

def addcloudcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.39
    bottom = axes_bbox.y0 - 0.015
    width = 0.38
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Total Cloud Cover (%)', size=8)  # MODIFY THIS for other fields!!
