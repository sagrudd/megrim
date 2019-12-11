#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:36:50 2019

on a plane easyjet TXL-LGW

# this module requires access to font-awesome
# this should be included in the same directory and in this case has been
# sourced from https://github.com/FortAwesome/Font-Awesome

# cheatsheet can be found at
# https://fontawesome.com/icons/?s=solid

@author: srudd
"""

import fontawesome as fa
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fp = FontProperties(fname=r"Font Awesome 5 Free-Solid-900.otf") 

class InfographicNode:
        
    def __init__(self, legend=None, value=None, graphic=None):
        """
        method initialising the infographic_node

        Parameters
        ----------
        legend : TYPE, optional
            the infographic legend. The default is None.
        value : TYPE, optional
            DESCRIPTION. The default is None.
        graphic : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        self.set_legend(legend)
        self.set_value(value)
        self.set_graphic(graphic)
        
    def set_legend(self, legend):
        """
        Set the legend used in infographic display

        Parameters
        ----------
        legend : string
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.legend = legend
        
    def set_value(self, value):
        self.value = value
        
    def set_graphic(self, graphic):
        self.graphic = graphic
        
    def get_legend(self):
        return self.legend
    
    def get_value(self):
        return self.value
    
    def get_graphic(self):
        return self.graphic
    
    def __main__(self):
        return self.toString()
    
    def to_string(self):
        return None
    






class InfographicPlot:
    
    def __init__(self, plot_content, rows=1, columns=3):
        self.plot_content = plot_content
        self.rows = rows
        self.columns = columns
        
    def plot_infographic(self):
        plt.figure(figsize=(self.rows*5, self.columns*5))
        f, axarr = plt.subplots(self.rows, self.columns, sharex=True, sharey=True)
        x = 0
        y = 0
        fontA = 48 / (self.columns*1.15)
        fontB = 32 / (self.columns*1.15)
        ImageC = 96 / (self.columns*1.15)
        for infographic_node in self.plot_content: 
            rect = patches.Rectangle((0,0),1,1,linewidth=1,edgecolor='#2171b5',facecolor='#2171b5')
            print("x=%s, y=%s" % (x, y))
            if self.rows > 1:
                axarr[y, x].axis('off')
                axarr[y, x].axis([0,1,0,1])
                axarr[y, x].add_patch(rect)
                axarr[y, x].text(0.075, 0.45, infographic_node.get_value(), fontsize=fontA, color='#DEEBF7')
                axarr[y, x].text(0.075, 0.075, infographic_node.get_legend(), fontsize=fontB, color='#C6DBEF')
                axarr[y, x].text(.850, .875, infographic_node.get_graphic(), fontproperties=fp, size=ImageC, 
                             color="#6BAED6", ha="center", va="center")
                x += 1
                if x>=self.columns:
                    x = 0
                    y += 1
            else:
                axarr[x].axis('off')
                axarr[x].axis([0,1,0,1])
                axarr[x].add_patch(rect)
                axarr[x].text(0.075, 0.45, infographic_node.get_value(), fontsize=fontA, color='#DEEBF7')
                axarr[x].text(0.075, 0.075, infographic_node.get_legend(), fontsize=fontB, color='#C6DBEF')
                axarr[x].text(.850, .875, infographic_node.get_graphic(), fontproperties=fp, size=ImageC, 
                             color="#6BAED6", ha="center", va="center")
                x += 1
        plt.axis('off')
        #f.subplots_adjust(hspace=0.1)
        #plt.show()
        plt.savefig(fname="x.png", dpi=300)
        
    def __str__(self):
        return("InfographicPlot[%s,%s]" % (self.columns, self.rows))
        
        

ip = InfographicPlot(infographic_data)  
ip.plot_infographic()   
        



flowcell_node = InfographicNode(legend="flowcell", 
                                value="FAK85195", 
                                graphic=fa.icons['fingerprint'])

readcount_node = InfographicNode(legend="Reads produced", 
                                value="1,234,567",
                                graphic=fa.icons['filter'])

gb_seq_val = InfographicNode(legend="Gigabases called", 
                                value="13.1",
                                graphic=fa.icons['flag-checkered'])


infographic_data = [flowcell_node, readcount_node, gb_seq_val]






infographic_node = readcount_node

fig, ax = plt.subplots(facecolor='#2171b5')
ax.axis([0,1,0,1])

ax.text(0.075, 0.45, infographic_node.get_value(), fontsize=48, color='#DEEBF7')
ax.text(0.075, 0.075, infographic_node.get_legend(), fontsize=32, color='#C6DBEF')

ax.text(.850, .875, infographic_node.get_graphic(), fontproperties=fp, size=96, 
         color="#6BAED6", ha="center", va="center")

plt.axis('off')
plt.show()

