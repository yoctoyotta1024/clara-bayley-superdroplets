### Generic plotting functions usually for subplot axes
import numpy as np
import matplotlib.pyplot as plt


def save_figure(savedir, savename):

    plt.savefig(savedir+savename, 
            dpi=400, bbox_inches="tight", 
            facecolor='w', format="png")

    

def linear_twinax(ax, lnr, eps):
    ''' linear x axis for lognormal 
    radius distribution plots'''
     
    axb = ax.twiny()
   
    axb.plot(np.e**lnr*1e6, eps, alpha=0)
    axb.set_xscale('log')
    axb.set_xlabel('radius, r /\u03BCm)')

    axb.xaxis.tick_bottom() 
    axb.xaxis.set_label_position('bottom') 

    # xlims = [np.exp(l)*1e6 for l in ax.get_xlim()]
    # axb.set_xlim(xlims)

    return axb



def axplt(ax, x, y, 
    xlab=None, ylab=None, 
    lab="_nolegend_", c=0, l='-'):
    
    if type(c)==type(0):
        c= 'C'+str(c)
    
    ax.plot(x,y, label=lab, color=c, linestyle=l)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)


