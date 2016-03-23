"""
PyQt seam cells analysis GUI

NB: the python package tifffile (Gohlke) needs to be installed.

author: Nicola Gritti
last edited: June 2015
"""

import sys
# from tifffile import *
from timelapseFun import *
import pickle
import os
from PyQt4 import QtGui, QtCore
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends import qt_compat
import glob
import pandas as pd
use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE

class GUI(QtGui.QWidget):
    
    #-----------------------------------------------------------------------------------------------
    # INITIALIZATION OF THE WINDOW - DEFINE AND PLACE ALL THE WIDGETS
    #-----------------------------------------------------------------------------------------------

    def __init__(self):

        super(GUI, self).__init__()
        
        self.setWindowTitle('Seam Cells Analysis')
        self.scaleFactor = 4
        self.side = 'L'
        self.lbltxt = '"wheel" press: change side, currently %s\n"i" or "u" press: change cell sides'
        self.initUI()
        
    def initUI(self):
        
        # SET THE GEOMETRY
        
        mainWindow = QtGui.QVBoxLayout()
        mainWindow.setSpacing(15)
        
        fileBox = QtGui.QHBoxLayout()
        spaceBox1 = QtGui.QHBoxLayout()
        rawDataBox = QtGui.QHBoxLayout()
        spaceBox2 = QtGui.QHBoxLayout()
        straightBox = QtGui.QVBoxLayout()
        
        mainWindow.addLayout(fileBox)
        mainWindow.addLayout(spaceBox1)
        mainWindow.addLayout(rawDataBox)
        mainWindow.addLayout(spaceBox2)
        mainWindow.addLayout(straightBox)
        
        Col1 = QtGui.QGridLayout()
        Col2 = QtGui.QHBoxLayout()
        Col3 = QtGui.QVBoxLayout()
        Col4 = QtGui.QVBoxLayout()
        
        rawDataBox.addLayout(Col1)
        rawDataBox.addLayout(Col2)
        rawDataBox.addLayout(Col3)
        rawDataBox.addLayout(Col4)
        
        Raw1 = QtGui.QHBoxLayout()
        Raw2 = QtGui.QHBoxLayout()
        Raw3 = QtGui.QHBoxLayout()

        straightBox.addLayout(Raw1)
        straightBox.addLayout(Raw2)
        straightBox.addLayout(Raw3)

        self.setLayout(mainWindow)

        # DEFINE ALL WIDGETS AND BUTTONS
        
        loadBtn = QtGui.QPushButton('Load DataSet')
        saveBtn = QtGui.QPushButton('Save data (F12)')
        
        tpLbl = QtGui.QLabel('Timepoint:')
        slLbl = QtGui.QLabel('Slice:')
        hatchLbl = QtGui.QLabel('Hatching Time:')
        
        self.tp = QtGui.QSpinBox(self)
        self.tp.setValue(0)
        self.tp.setMaximum(100000)
        
        self.sl = QtGui.QSpinBox(self)
        self.sl.setValue(0)
        self.sl.setMaximum(100000)
        
        self.hatch = QtGui.QSpinBox(self)
        self.hatch.setValue(0)
        self.hatch.setMaximum(100000)
        
        self._488nmBtn = QtGui.QRadioButton('488nm')
        self._561nmBtn = QtGui.QRadioButton('561nm')
        self.CoolLEDBtn = QtGui.QRadioButton('CoolLED')
        
        automaticOutlineBtn = QtGui.QPushButton('Automatic Outline')
        straightenBtn = QtGui.QPushButton('Straighten images')
        
        self.sld1 = QtGui.QSlider(QtCore.Qt.Vertical, self)
        self.sld1.setMaximum(2**16-1)
        self.sld1.setValue(0)
        self.sld2 = QtGui.QSlider(QtCore.Qt.Vertical, self)
        self.sld2.setMaximum(2**16)
        self.sld2.setValue(2**16-1)

        self.fig1 = Figure((8.0, 8.0), dpi=100)
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        self.ax1 = self.fig1.add_subplot(111)
        self.canvas1 = FigureCanvas(self.fig1)
        self.canvas1.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas1.setFocus()
        self.canvas1.setFixedSize(QtCore.QSize(500,500)) ############################# change this value to set the figure size
        self.canvas1.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )

        self.cellTbl = QtGui.QTableWidget()

        self.sideLbl = QtGui.QLabel(self.lbltxt % self.side)
        autonameBtn = QtGui.QPushButton('Restore Automatic Cell Names')
        wrt2Btn = QtGui.QPushButton('Compute wrt-2 Expression')
        wrt2FullBtn = QtGui.QPushButton('Compute wrt-2 Expression on the raw data (SLOW!)')        
        
        self.fig2 = Figure((4.0, 2.0))        
        self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        self.ax2 = self.fig2.add_subplot(111)
        self.canvas2 = FigureCanvas(self.fig2)
        self.canvas2.setFocusPolicy( QtCore.Qt.StrongFocus )
        self.canvas2.setFocus()
        self.canvas2.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )

        self.fig3 = Figure((4.0, 2.0), dpi=100)
        self.fig3.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        self.ax3 = self.fig3.add_subplot(111)
        self.canvas3 = FigureCanvas(self.fig3)
        self.canvas3.setFocusPolicy( QtCore.Qt.ClickFocus )

        self.canvas3.setFocus()
        self.canvas3.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )
        
        # PLACE ALL THE WIDGET ACCORDING TO THE GRIDS

        fileBox.addWidget(loadBtn)
        fileBox.addWidget(saveBtn)

        spaceBox1.addWidget(self.HLine())

        Col1.addWidget(tpLbl, 0, 0)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self.tp, 0, 1)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(slLbl, 1, 0)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self.sl, 1, 1)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(hatchLbl, 2, 0)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self.hatch, 2, 1)#, 1, 1, Qt.AlignTop)
        Col1.addWidget(self._488nmBtn, 3, 0 )
        Col1.addWidget(self._561nmBtn, 4, 0 )
        Col1.addWidget(self.CoolLEDBtn, 5, 0 )
        Col1.addWidget(automaticOutlineBtn, 6, 0 )        
        Col1.addWidget(straightenBtn, 7, 0 )
        
        Col2.addWidget(self.sld1)
        Col2.addWidget(self.sld2)
        Col2.addWidget(self.canvas1)
        
        Col3.addWidget(self.VLine())
        
        Col4.addWidget(self.cellTbl)

        spaceBox2.addWidget(self.HLine())

        Raw1.addWidget(self.sideLbl)
        Raw1.addWidget(autonameBtn)
        Raw1.addWidget(wrt2Btn)
        Raw1.addWidget(wrt2FullBtn)

        Raw2.addWidget(self.canvas2)
        
        Raw3.addWidget(self.canvas3)
        # Raw3.addWidget(QtGui.QDockWidget())

        self.setFocus()
        self.show()
        
        # BIND BUTTONS TO FUNCTIONS
        
        loadBtn.clicked.connect(self.selectWorm)
        saveBtn.clicked.connect(self.saveData)

        self.tp.valueChanged.connect(self.loadNewStack)
        self.sl.valueChanged.connect(self.updateAllCanvas)
        self.sld1.valueChanged.connect(self.updateAllCanvas)
        self.sld2.valueChanged.connect(self.updateAllCanvas)
        self.hatch.valueChanged.connect(self.updateTidxDataFrame)

        self._488nmBtn.toggled.connect(self.radioClicked)
        self._561nmBtn.toggled.connect(self.radioClicked)
        self.CoolLEDBtn.toggled.connect(self.radioClicked)

        automaticOutlineBtn.clicked.connect(self.automaticOutline)
        straightenBtn.clicked.connect(self.straightenWorm)
        autonameBtn.clicked.connect(self.automaticCellNames)

        self.fig1.canvas.mpl_connect('button_press_event',self.onMouseClickOnCanvas1)        
        self.fig2.canvas.mpl_connect('button_press_event',self.onMouseClickOnCanvas2)        
        self.fig1.canvas.mpl_connect('scroll_event',self.wheelEvent)        
        self.fig2.canvas.mpl_connect('scroll_event',self.wheelEvent)        
        self.fig3.canvas.mpl_connect('scroll_event',self.wheelEvent)        
        
    #-----------------------------------------------------------------------------------------------
    # FORMATTING THE WINDOW
    #-----------------------------------------------------------------------------------------------

    def center(self):
        
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    def HLine(self):
        
        toto = QtGui.QFrame()
        toto.setFrameShape(QtGui.QFrame.HLine)
        toto.setFrameShadow(QtGui.QFrame.Sunken)
        return toto

    def VLine(self):
        
        toto = QtGui.QFrame()
        toto.setFrameShape(QtGui.QFrame.VLine)
        toto.setFrameShadow(QtGui.QFrame.Sunken)
        return toto

    def heightForWidth(self, width):
        
        return width
    
    #-----------------------------------------------------------------------------------------------
    # BUTTON FUNCTIONS
    #-----------------------------------------------------------------------------------------------

    def selectWorm(self):

        self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'Y:\\Images')
        self.worm = self.pathDial.split("\\")[-1]
        self.path = self.pathDial[:-len(self.worm)]
        
        self.setWindowTitle('Seam Cells Analysis - ' + self.pathDial)
        
        if 'downsized' not in self.worm:
            QtGui.QMessageBox.about(self,'Warning!','These are not downsized (512x512) images! Convert the images first!')
            return

        # print(path + worm + '\\z*.txt')
        self.fList = {}
        if os.path.isfile(self.path + self.worm + '\\z001_488nm.tif'):
            self.fList['488nm'] = glob.glob(self.path + self.worm + '\\z*488nm.tif')
            self.fList['488nm'].sort()
        if os.path.isfile(self.path + self.worm + '\\z001_561nm.tif'):
            self.fList['561nm'] = glob.glob(self.path + self.worm + '\\z*561nm.tif')
            self.fList['561nm'].sort()
        if os.path.isfile(self.path + self.worm + '\\z001_CoolLED.tif'):
            self.fList['CoolLED'] = glob.glob(self.path + self.worm + '\\z*CoolLED.tif')
            self.fList['CoolLED'].sort()

        self.fMetaList = glob.glob(self.path + self.worm + '\\z*.txt')
        
        self.channel = list(self.fList.keys())[0]
        self.tp.setMaximum(len(self.fList[self.channel])-1)
        
        if not os.path.isfile( self.path + '\\worm' + self.worm.split('_')[0] + '.pickle' ):

            create_worm( self.path, self.worm.split('_')[0], 40, len(self.fList[self.channel]) - self.hatch.value(), self.hatch.value()+1 )

        self.df = pickle.load( open(self.path + '\\worm' + self.worm.split('_')[0] + '.pickle','rb') )

        # if there are more timepoints, add the body rows to the dataframe
        if np.max(self.df.tidx)-np.min(self.df.tidx) < len(self.fMetaList):

            times = extract_times( self.path+'\\'+self.worm, -np.min(self.df.tidx) )[int(np.max(self.df.tidx)-np.min(self.df.tidx)+1):]
            df1 = pd.DataFrame( { 'rowtype': 'body',
                          'tidx': np.arange(np.max(self.df.tidx)+1,len(self.fMetaList)+np.min(self.df.tidx)),
                          'times': times,
                          'outline': np.nan,
                          'spline': np.nan} )

            df1.outline = df1.outline.astype(object)
            df1.spline = df1.spline.astype(object)
            for idx, row in df1.iterrows():    
                df1.outline.values[idx] = np.array([[np.nan,np.nan,np.nan]])
                df1.spline.values[idx] = np.array([[np.nan,np.nan,np.nan]])

            self.df = pd.concat([self.df,df1]).sort(['tidx','rowtype','cside','cXpos']).reset_index(drop=True)

        self.hatch.setValue(-np.min(self.df.tidx))
        self.tp.setValue( self.hatch.value() )

        self.loadNewStack()

        # self.pathDial.show()
        self.setFocus()

    def saveData(self):
        
        pickle.dump( self.df, open(self.path+'\\worm'+self.worm.split('_')[0]+'.pickle','wb'), protocol=2 )        
        
    def loadNewStack(self):
        
        # print(self.fList['gfp'][self.tp.value()])
        
        self.stacks = {}
        for key in self.fList.keys():
            self.stacks[key] = loadstack(self.fList[key][self.tp.value()])
        self.scaleFactor = 2048 / self.stacks[self.channel][0].shape[0]

        self.stacksStraight = {}
        self.straightFile = self.path+self.worm.split('_')[0]+'_straighten\\straight%.3d_%s.tif'%(self.tp.value()+1,self.channel)

        if os.path.isfile(self.straightFile):

            for key in self.fList.keys():
                self.stacksStraight[key] = loadstack(self.path+self.worm.split('_')[0]+'_straighten\\straight%.3d_%s.tif'%(self.tp.value()+1,key))
        # print(self.stacks.keys(), self.stacksStraight)
        self.sl.setMaximum(self.stacks[self.channel].shape[0]-1)

        self.setBCslidersMinMax()
        
        self.updateTable()
        self.updateAllCanvas()

    def updateAllCanvas(self):
        self.updateRadioBtn()
        self.updateCanvas1()
        self.updateCanvas2()
        self.updateCanvas3()
        
    def updateTidxDataFrame(self):
        
        times = extract_times(self.path+self.worm, self.hatch.value()+1)
        self.df.ix[self.df.rowtype=='body','tidx'] = np.arange(len(times))-self.hatch.value()
        self.df.ix[self.df.rowtype=='body','times'] = times

    def radioClicked(self):
        if self._488nmBtn.isChecked():
            if '488nm' in self.fList.keys():
                self.channel = '488nm'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No 488nm channel!')
        elif self._561nmBtn.isChecked():
            if '561nm' in self.fList.keys():
                self.channel = '561nm'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No 561nm channel!')
        elif self.CoolLEDBtn.isChecked():
            if 'CoolLED' in self.fList.keys():
                self.channel = 'CoolLED'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No CoolLED channel!')
        self.setBCslidersMinMax()
        self.resetBC()
        self.setFocus()
        self.updateAllCanvas()

    def automaticOutline(self):
        
        tp = self.tp.value()
        sl = self.sl.value()
        
        rowmask = self.df['rowtype']=='body'
        outline = self.df.ix[ rowmask, 'outline' ].values
        spline = self.df.ix[ rowmask, 'spline' ].values
        
        for idx, o in enumerate( outline ):

            if len(o) == 2:
                
                stack = loadstack( self.fList['488nm'][idx] )
                
                # automatically create the outline
                outline[idx] = automatic_outline_michalis( stack, o, sl, idx )

                # create the spline interpolation
                spline[idx] = interp_spline( outline[idx] )

        # update the dataframe with the new outline and spline
        self.df.ix[ rowmask, 'outline' ] = outline
        self.df.ix[ rowmask, 'spline' ] = spline
        self.updateCanvas1()
    
    def straightenWorm(self):

        raw_worm = self.worm.split('_')[0]
        path = self.path+raw_worm
        if not os.path.exists(path):
            QtGui.QMessageBox.about(self, 'Warning! No raw data found in the standard directory. No straightening will be performed', 'Y:\\Images')
            return
        
        spline = []
        for idx in self.df.ix[ pd.notnull(self.df.spline)].index:
            spline.append( np.array( [ self.df.spline.values[idx][:,0]*self.scaleFactor,
                                    self.df.spline.values[idx][:,1]*self.scaleFactor,
                                    self.df.spline.values[idx][:,2] ] ).T )

        fList = {}
        if os.path.isfile(path + '\\z001_488nm.tif'):
            fList['488nm'] = glob.glob(path + '\\z*488nm.tif')
            fList['488nm'].sort()
        if os.path.isfile(path+'\\z001_561nm.tif'):
            fList['561nm'] = glob.glob(path + '\\z*561nm.tif')
            fList['561nm'].sort()
        if os.path.isfile(path + '\\z001_CoolLED.tif'):
            fList['CoolLED'] = glob.glob(path + '\\z*CoolLED.tif')
            fList['CoolLED'].sort()

#        print(25*self.scaleFactor)
        for key in list( fList.keys() ):
            straighten( path, fList[str(key)], spline )

        self.loadNewStack()
        self.updateCanvas2()
        self.canvas1.setFocus()
         
    def automaticCellNames(self):

        tidxNow = self.tp.value() - self.hatch.value()

        tidxNowMask = self.df.tidx == tidxNow
        cellMask = self.df.rowtype == 'cell'
        sideMask = self.df.cside == self.side

        # find the previous timepoint with labeled cells
        if np.sum( sideMask & (self.df.tidx < tidxNow) ) == 0:
            QtGui.QMessageBox.about(self, 'Warning', 'No cells are labeled in the %s side yet!' % self.side)
            return   
        else:
            tidxPrev = np.max( self.df.ix[ sideMask & ( self.df.tidx<tidxNow ), 'tidx' ] )
            # print(self.side, tidxPrev)
            tidxPrevMask = self.df['tidx'] == tidxPrev

        # filter the cells from the dataframe
        cellsNow = self.df[cellMask&tidxNowMask&sideMask]
        cellsPrev = self.df[cellMask&tidxPrevMask&sideMask]

        # find the relative positions and sort the cells
        wNowLen = np.max(cellsNow.cXpos)-np.min(cellsNow.cXpos)
        wPrevLen = np.max(cellsPrev.cXpos)-np.min(cellsPrev.cXpos)

        for idx in cellsNow.index:
            cellsNow.ix[idx,'relPos'] = ( cellsNow.ix[idx,'cXpos'] - np.min(cellsNow.cXpos) )/wNowLen
        for idx in cellsPrev.index:
            cellsPrev.ix[idx,'relPos'] = ( cellsPrev.ix[idx,'cXpos'] - np.min(cellsPrev.cXpos) )/wPrevLen

        cellsNow = cellsNow.sort(['relPos'])
        cellsPrev = cellsPrev.sort(['relPos'])

        # assign all the names according to closest cell in previous timepoint
        # if there are more cells now, two cells will have the same name
        for idx, cell in cellsNow.iterrows():
            closestCell = self.closestCell( cell, cellsPrev )
            cellsNow.ix[idx,'cname'] = closestCell.cname.values[0]

        # if more cells, a division happened.
        if len(cellsNow) > len(cellsPrev):
            # print('I am fucked...')

            # for each cell (but first the most left one), find the first one to the left
            for idx, cell in cellsNow.drop( np.min(cellsNow.index.values) ).iterrows():
                cellToTheLeft = cellsNow.ix[ idx - 1 ]

                # if it has the same name, a division happened
                if cellsNow.ix[idx,'cname'] == cellToTheLeft.cname:
                    cellsNow.ix[idx-1,'cname'] += 'a'
                    cellsNow.ix[idx,'cname'] += 'p'

        # update dataframe
        for idx, cell in cellsNow.iterrows():
            self.df.ix[idx,'cname'] = cell.cname

        # update canvas
        self.updateTable()
        self.updateCanvas2()

    #-----------------------------------------------------------------------------------------------
    # DEFAULT FUNCTION FOR KEY AND MOUSE PRESS ON WINDOW
    #-----------------------------------------------------------------------------------------------

    def keyPressEvent(self, event):
        
        # print(event.key())

        # change timepoint
        if event.key() == QtCore.Qt.Key_Right:
            self.changeSpaceTime( self.tp, +1 )

        if event.key() == QtCore.Qt.Key_Left:
            self.changeSpaceTime( self.tp, -1 )

        # change slice
        if event.key() == QtCore.Qt.Key_Up:
            self.changeSpaceTime( self.sl, +1 )
            
        if event.key() == QtCore.Qt.Key_Down:
            self.changeSpaceTime( self.sl, -1 )

        # change channel
        if event.key() == QtCore.Qt.Key_Space:
            currentidx = list(self.fList.keys()).index(self.channel) 
            nextidx = np.mod(currentidx+1,len(self.fList.keys()))
            self.channel = list(self.fList.keys())[nextidx]
            self.setBCslidersMinMax()
            self.resetBC()
            self.updateAllCanvas()
            
        # key press on straighten worm
        if self.canvas2.underMouse():
            self.onKeyPressOnCanvas2(event)
            
        self.setFocus()
        
    def wheelEvent(self,event):
        if any([ self.canvas1.underMouse(), self.canvas2.underMouse(), self.canvas3.underMouse() ]):
            step = event.step
        else:          
            step = event.delta()/abs(event.delta())
        self.sl.setValue( self.sl.value() + step) 

    #-----------------------------------------------------------------------------------------------
    # ADDITIONAL FUNCTIONS FOR KEY AND MOUSE PRESS ON CANVASES
    #-----------------------------------------------------------------------------------------------

    def onKeyPressOnCanvas2(self, event):

        # print(event.key())
        
        cellsname = [ QtCore.Qt.Key_A, QtCore.Qt.Key_B, QtCore.Qt.Key_C, QtCore.Qt.Key_1,
                      QtCore.Qt.Key_2, QtCore.Qt.Key_3, QtCore.Qt.Key_4, QtCore.Qt.Key_5,
                      QtCore.Qt.Key_6, QtCore.Qt.Key_T ]
        cellspos = [ QtCore.Qt.Key_Q, QtCore.Qt.Key_W ]

        # find the position of the cursor relative to the image in pixel
        imgshape = self.stacksStraight[self.channel][self.sl.value()].shape
        arimg = imgshape[1]/imgshape[0]
        canshape = self.canvas2.size()
        arcan = canshape.width()/canshape.height()
        cf = 1
        if arimg>arcan:
            cf = imgshape[1]/canshape.width()
            origin = ( canshape.height()*cf - imgshape[0] ) / 2.
            origin = np.array( [ 0, ( canshape.width()*cf - imgshape[1] ) / 2. ] )
        else:
            cf = imgshape[0]/canshape.height()
            origin = np.array( [ ( canshape.width()*cf - imgshape[1] ) / 2., 0 ] )
        refpos = self.canvas2.mapFromGlobal(QtGui.QCursor.pos())
        refpos = np.array([refpos.x(),refpos.y()])*cf - origin
        refpos = np.append(refpos,self.sl.value())

        # find the closest cell to the cursor
        bodymask = self.df['rowtype']=='body'
        cellmask = self.df['rowtype']=='cell'
        tpmask = self.df['tidx'] == (self.tp.value()-self.hatch.value())
        
        idx, cell = closer_cell( refpos, self.df[ cellmask & tpmask ], self.fMetaList[self.tp.value()] )

        if any( [ event.key() == cn for cn in cellsname ] ):
            self.df.ix[ idx, 'cname' ] = QtGui.QKeySequence(event.key()).toString().lower() + '.'

        elif any( [ event.key() == cp for cp in cellspos ] ):

            if event.key() == cellspos[0]:      self.df.ix[ idx, 'cname' ] += 'a'
            elif event.key() == cellspos[1]:    self.df.ix[ idx, 'cname' ] += 'p'
                
        elif event.key() == QtCore.Qt.Key_Backspace:
            
            self.df.ix[ idx, 'cname' ] = self.df.ix[ idx, 'cname' ][:-1]

        elif event.key() == QtCore.Qt.Key_I:
            
            sideLmask = self.df['cside']=='L'
            self.df.ix[ tpmask & cellmask & sideLmask,'cside'] = 'R'
            sideRmask = self.df['cside']=='R'
            self.df.ix[ tpmask & cellmask & sideRmask,'cside'] = 'L'

        elif event.key() == QtCore.Qt.Key_U:
            
            if self.df.ix[idx,'cside']=='L':    self.df.ix[ idx,'cside'] = 'R'
            elif self.df.ix[idx,'cside']=='R':  self.df.ix[ idx,'cside'] = 'L'
        
        self.df = self.df.sort(['tidx','rowtype','cside','cXpos']).reset_index(drop=True)
        self.updateTable()
        self.updateAllCanvas()
        
    def onMouseClickOnCanvas1(self, event):
        
        # print(event.button,event.xdata,event.ydata)
        
        tp = self.tp.value()
        sl = self.sl.value()
        
        rowmask = self.df['rowtype']=='body'
        outline = self.df.ix[ rowmask, 'outline' ].values
        spline = self.df.ix[ rowmask, 'spline' ].values

        x = event.xdata
        y = event.ydata

        # left button: add a point to the outline
        if event.button == 1:   
            
            if not np.all(np.isnan(outline[tp])):
                idx = find_closest_point( [ x, y ], outline[tp] )
            else:
                idx = 0                
            outline[tp] = np.insert( outline[tp], idx+1, [ x, y, sl ], axis=0 )
            
            # if the first line is still full of nan, remove it
            if np.all(np.isnan(outline[tp][0])):
                outline[tp] = np.delete(outline[tp],0,axis=0)
                                
        # right button: remove the closest point from the outline
        elif event.button == 3 and len(outline[tp]) > 0:
            
            idx = find_closest_point([x,y],outline[tp])
            if outline[tp].shape[0]!=1:
                outline[tp] = np.delete( outline[tp], idx, axis=0 )
            else:
                outline[tp] = np.array([[np.nan,np.nan,np.nan]])
                
        # update the dataframe with the new outline and spline
        spline[tp] = interp_spline( outline[tp] )
        self.df.ix[ rowmask, 'outline' ] = outline
        self.df.ix[ rowmask, 'spline' ] = spline
        self.updateCanvas1()
        self.setFocus()

    def onMouseClickOnCanvas2(self, event):
        
        refpos = np.array( [ event.xdata, event.ydata, self.sl.value() ] )  
        # print(refpos)
        bodymask = self.df['rowtype']=='body'
        cellmask = self.df['rowtype']=='cell'
        tpmask = self.df['tidx'] == (self.tp.value()-self.hatch.value())
        # print( 'mouse button pressed:', event.button, 'at pos:', refpos )
        if all( refpos[:-1] ) > 0:
            if event.button == 1:
    
                # create an empty cell: the only entries are tidx, times, xyzpos, side
                newcell = create_cell( refpos, self.df[ bodymask & tpmask ], self.side  )
                self.df = pd.concat( [ self.df, newcell ] )
                
            elif event.button == 3:

                if any( self.df[tpmask].rowtype == 'cell' ):
                    idx, cell = closer_cell( refpos, self.df[ cellmask & tpmask ], self.fMetaList[self.tp.value()] )
                    self.df = self.df.drop([idx])
            
            elif event.button == 2:
    
                if self.side == 'L': self.side = 'R'
                elif self.side == 'R': self.side = 'L'
                self.sideLbl.setText(self.lbltxt % self.side)

        self.df = self.df.sort(['tidx','rowtype','cside','cXpos']).reset_index(drop=True)
        self.updateCanvas2()
        self.updateTable()
        self.setFocus()
                
    #-----------------------------------------------------------------------------------------------
    # UTILS
    #-----------------------------------------------------------------------------------------------

    def updateRadioBtn(self):
        if self.channel == '488nm':
            self._488nmBtn.setChecked(True)
        elif self.channel == '561nm':
            self._561nmBtn.setChecked(True)
        elif self.channel == 'CoolLED':
            self.CoolLEDBtn.setChecked(True)
        self.setFocus()

    def setBCslidersMinMax(self):
        self.sld1.setMaximum(np.max(self.stacks[self.channel]))
        self.sld1.setMinimum(np.min(self.stacks[self.channel]))
        self.sld2.setMaximum(np.max(self.stacks[self.channel]))
        self.sld2.setMinimum(np.min(self.stacks[self.channel]))

    def resetBC(self):
        self.sld1.setValue(np.min(self.stacks[self.channel]))
        self.sld2.setValue(np.max(self.stacks[self.channel]))
        
    def updateCanvas1(self):
        
        rowmask = self.df['rowtype']=='body'
        tidxmask = self.df['tidx']==(self.tp.value()-self.hatch.value())
        
        # extract and rescale the outline and spline
        outline = np.copy( self.df.ix[ rowmask & tidxmask, 'outline' ].values[0] )
        spline = np.copy( self.df.ix[ rowmask & tidxmask, 'spline' ].values[0] )

        # plot the image
        self.ax1.cla()
        imgplot = self.ax1.imshow(self.stacks[self.channel][self.sl.value()], cmap = 'gray')
        
        # remove the white borders and plot outline and spline
        self.ax1.autoscale(False)
        self.ax1.axis('Off')
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        # change brightness and contrast
        self.sld1.setValue(np.min([self.sld1.value(),self.sld2.value()]))
        self.sld2.setValue(np.max([self.sld1.value(),self.sld2.value()]))
        imgplot.set_clim(self.sld1.value(), self.sld2.value())        

        # print(outline, spline)
        self.ax1.plot( outline[:,0], outline[:,1], 'o', color='red', ms=6, mew=1, alpha=.5, lw = 1 )
        self.ax1.plot( spline[:,0], spline[:,1], '-', color='yellow', lw = 1 )
        
        # redraw the canvas
        self.canvas1.draw()
        self.setFocus()

    def updateCanvas2(self):
        
        # plot the image
        self.ax2.cla()

        if not os.path.isfile(self.straightFile):
            self.fig2.clf()
            self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
            self.ax2 = self.fig2.add_subplot(111)
            self.canvas2.draw()
            return
            
        imgplot = self.ax2.imshow(self.stacksStraight[self.channel][self.sl.value()], cmap = 'gray')
        imgplot.set_clim(self.sld1.value(), self.sld2.value())

        # remove the white borders
        self.ax2.autoscale(False)
        self.ax2.axis('Off')
        self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        
        # cell text on the image
        tpmask = self.df['tidx'] == (self.tp.value() - self.hatch.value())
        cellmask = self.df['rowtype'] == 'cell'

        for idx, cell in self.df[tpmask & cellmask].iterrows():

            if cell.cZpos == self.sl.value():
                clabel = str(cell.cname) + ' ' + str(cell.cside)
                self.ax2.text( cell.cXpos, cell.cYpos + 10, clabel, color='red', size='medium', alpha=.8,
                        rotation=90)
                self.ax2.plot( cell.cXpos, cell.cYpos, 'x', color='red', alpha = .8 )


        # redraw the canvas
        self.canvas2.draw()
        self.setFocus()

    def updateCanvas3(self):
        # print('updating canvas 3')
        
        tidxNow = self.tp.value() - self.hatch.value()

        tidxNowMask = self.df.tidx == tidxNow
        cellMask = self.df.rowtype == 'cell'
        sideMask = self.df.cside == self.side

        # find the previous timepoint with labeled cells
        if np.sum( sideMask & (self.df.tidx < tidxNow) ) == 0:
            self.fig3.clf()
            self.fig3.subplots_adjust(left=0., right=1., top=1., bottom=0.)
            self.ax3 = self.fig3.add_subplot(111)
            self.canvas3.draw()
            return   
        else:
            tidxPrev = np.max( self.df.ix[ sideMask & ( self.df.tidx<tidxNow ), 'tidx' ] )
            tidxPrevMask = self.df['tidx'] == tidxPrev

        # load images
        prevstacksStraight = {}
        for key in self.fList.keys():
            prevstacksStraight[key] = loadstack(self.path+self.worm.split('_')[0]+'_straighten\\straight%.3d_%s.tif'%(tidxPrev+self.hatch.value()+1,key))

        # plot the image
        self.ax3.cla()

        imgplot = self.ax3.imshow(prevstacksStraight[self.channel][self.sl.value()], cmap = 'gray')
        imgplot.set_clim(self.sld1.value(), self.sld2.value())

        # remove the white borders
        self.ax3.autoscale(False)
        self.ax3.axis('Off')
        self.fig3.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        
        # cell text on the image
        for idx, cell in self.df[tidxPrevMask & cellMask].iterrows():

            if cell.cZpos == self.sl.value():
                clabel = str(cell.cname) + ' ' + str(cell.cside)
                self.ax3.text( cell.cXpos, cell.cYpos + 10, clabel, color='red', size='small', alpha=.8,
                        rotation=90)
                self.ax3.plot( cell.cXpos, cell.cYpos, 'x', color='red', alpha = .8 )

        # redraw the canvas
        self.canvas3.draw()
        self.setFocus()

    def closestCell(self,cell,clist):
        dist = np.abs( clist.relPos - cell.relPos )
        return clist[ dist == np.min(dist) ]
  
    def changeSpaceTime(self, whatToChange, increment):

        whatToChange.setValue( whatToChange.value() + increment )
        
    def rescaledImage(self, img):
        
        Nbig = img.shape[0]
        Nsmall = img.shape[0]/self.scaleFactor
        return img.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1)

    def updateTable(self):
        self.cellTbl.clear()

        tidxNow = self.tp.value() - self.hatch.value()

        tidxNowMask = self.df.tidx == tidxNow
        cellMask = self.df.rowtype == 'cell'

        cellsNow = self.df[cellMask&tidxNowMask]
        cellsPrev = pd.DataFrame({})
        if len(cellsNow) > 0:
	        sideMask = self.df.cside == cellsNow.cside.values[0]

	        # find the previous timepoint with labeled cells
	        if np.sum( sideMask & (self.df.tidx < tidxNow) ) == 0:
	        	cellsPrev = pd.DataFrame({})
	        else:
	            tidxPrev = np.max( self.df.ix[ sideMask & ( self.df.tidx<tidxNow ), 'tidx' ] )
	            tidxPrevMask = self.df['tidx'] == tidxPrev
	            cellsPrev = self.df[cellMask&tidxPrevMask]
	      
	        horHeaders = ['tidx','times','cell name','cell side','cellXpos','cellYpos','cellZpos','-',
	        				'tidx','times','cell name','cell side','cellXpos','cellYpos','cellZpos']
	        self.cellTbl.setColumnCount(len(horHeaders))
	        self.cellTbl.setRowCount(np.max([len(cellsNow),len(cellsPrev)]))        
	        self.cellTbl.setHorizontalHeaderLabels(horHeaders)
	        self.cellTbl.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
	        self.cellTbl.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)

	        row = 0
	        for idx, cell in cellsNow.iterrows():
	            self.cellTbl.setItem(row,0,QtGui.QTableWidgetItem(str(int(cell.tidx)),.0001))
	            self.cellTbl.setItem(row,1,QtGui.QTableWidgetItem(str('%.2f'%cell.times)))
	            self.cellTbl.setItem(row,2,QtGui.QTableWidgetItem(str(cell.cname)))
	            self.cellTbl.setItem(row,3,QtGui.QTableWidgetItem(cell.cside))
	            self.cellTbl.setItem(row,4,QtGui.QTableWidgetItem(str(int(cell.cXpos))))
	            self.cellTbl.setItem(row,5,QtGui.QTableWidgetItem(str(int(cell.cYpos))))
	            self.cellTbl.setItem(row,6,QtGui.QTableWidgetItem(str(int(cell.cZpos))))
	            row += 1

	        row = 0
	        for idx, cell in cellsPrev.iterrows():
	            self.cellTbl.setItem(row,8,QtGui.QTableWidgetItem(str(int(cell.tidx))))
	            self.cellTbl.setItem(row,9,QtGui.QTableWidgetItem(str('%.2f'%cell.times)))
	            self.cellTbl.setItem(row,10,QtGui.QTableWidgetItem(str(cell.cname)))
	            self.cellTbl.setItem(row,11,QtGui.QTableWidgetItem(cell.cside))
	            self.cellTbl.setItem(row,12,QtGui.QTableWidgetItem(str(int(cell.cXpos))))
	            self.cellTbl.setItem(row,13,QtGui.QTableWidgetItem(str(int(cell.cYpos))))
	            self.cellTbl.setItem(row,14,QtGui.QTableWidgetItem(str(int(cell.cZpos))))
	            row += 1
        
        self.setFocus()

if __name__ == '__main__':
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication(sys.argv)
    
    gui = GUI()
    app.setStyle("plastique")
    # app.installEventFilter(gui)
    sys.exit(app.exec_())
    


