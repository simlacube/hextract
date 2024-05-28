#!/usr/bin/env python
# coding: utf-8

# In[5]:


import hexGridShapelyDec
import matplotlib.pyplot as plt


# In[10]:


def saveRegionFile(worldoverlapRA,worldoverlapDec,regionFilename, rad):
    '''
    worldoverlapRA: first output from hexGridShapelyDec.multiFileApCut()
    worldoverlapDec: second output from hexGridShapelyDec.multiFileApCut()
    regionFilename: the name of the region file you want to create
    rad: radius in arcsec
    '''
    ### create region file
    with open(regionFilename, 'a') as f:
        top = '# Region file format: DS9 version 4.1'
        middle = 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
        bottom = 'fk5'
        f.write(top)
        f.write('\n')
        f.write(middle)
        f.write('\n')
        f.write(bottom)

        if isinstance(worldoverlapRA,float):
            f.write('circle(nan.nan,nan)')
        else:
            for i in range(len(worldoverlapRA)):
                coord = 'circle('+str(worldoverlapRA[i])+','+str(worldoverlapDec[i])+','+str(rad)+'")'
                f.write('\n')
                f.write(coord) 


# In[7]:


### here i am setting path+file name equal to a variable just to keep things a bit neater
ngc6822_LL1 = '/home/lhands/Gold/fitsfiles/159/NGC6822/cube_gold_NGC6822_LL1.fits'
ngc6822_LL2 = '/home/lhands/Gold/fitsfiles/159/NGC6822/cube_gold_NGC6822_LL2.fits'
ngc6822_LL1_unc = '/home/lhands/Gold/fitsfiles/159/NGC6822/cube_gold_NGC6822_LL1_unc.fits'
ngc6822_LL2_unc = '/home/lhands/Gold/fitsfiles/159/NGC6822/cube_gold_NGC6822_LL2_unc.fits'


### make a list of the file names that you want to find overlap between and cast the next of apertures on  
### it could contain any number of fits files
filelist = [ngc6822_LL1,ngc6822_LL2,ngc6822_LL1_unc,ngc6822_LL2_unc]

### here you can set the radius of the circular apertures in arcseconds
radius = 25 #arcsec 


### here you call the function that returns the hexagonal grid 
### in a list of RA coordinates and a list of Dec coordinates 
worldoverlapRA,worldoverlapDec = hexGridShapelyDec.multiFileApCut(filelist,radius)


### here you can save the hexagonal grid as a ds9 region file
regionFilename = 'hexgrid.reg'
saveRegionFile(worldoverlapRA,worldoverlapDec,regionFilename, radius)


# In[ ]:





# In[ ]:




