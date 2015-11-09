#!/usr/local/bin/python2.7
# Licence: GPLv2.0

import sys,os

import pygtk
pygtk.require("2.0")

import matplotlib
matplotlib.use('GTK')

import gtk
import gtk.glade

import matplotlib.cm as cm

#from matplotlib.figure import *
#from matplotlib.axes import *
#from matplotlib.backends.backend_gtk import *
#from pylab import *

from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas

from exp_setup import *

class app:
    target = ""
    def __init__(self):
        self.init_app()

    def init_app(self):
        builder = gtk.Builder()
        builder.add_from_file(PYXS_PATH+"/pyXS-gui-2.glade")
        builder.connect_signals(self)
        self.window1 = builder.get_object("window1")
        self.window1.show()     

        self.fileList=builder.get_object("liststore1")
        self.fileTv=builder.get_object("tvFileList")
        self.fig2d=None
        self.fig1d=None
        self.fn2d=None
        self.ax1=None
        self.ax2=None

        self.btnMaskOn=builder.get_object("checkbuttonMaskOn")
        self.showMask=False
        self.btnMaskOn.set_active(self.showMask)

        self.btnRingOn=builder.get_object("checkbuttonRingsOn")
        self.showRing=False
        self.btnRingOn.set_active(self.showRing)

        self.SAXScLo=0
        self.SAXScHi=200
        self.WAXScLo=10
        self.WAXScHi=8000
        self.etSAXScLo=builder.get_object("entrySAXScolorLo")
        self.etSAXScHi=builder.get_object("entrySAXScolorHi")
        self.etWAXScLo=builder.get_object("entryWAXScolorLo")
        self.etWAXScHi=builder.get_object("entryWAXScolorHi")
        self.etSAXScLo.set_text("%d" % self.SAXScLo)
        self.etSAXScHi.set_text("%d" % self.SAXScHi)
        self.etWAXScLo.set_text("%d" % self.WAXScLo)
        self.etWAXScHi.set_text("%d" % self.WAXScHi)

        self.etOutFn=builder.get_object("entryOutputFile")

        self.refreshFileList(self)

    def reloadExpSetup(self,widget):
        if os.path.exists(FLAT_FILE_S): os.remove(FLAT_FILE_S)
        if os.path.exists(FLAT_FILE_W): os.remove(FLAT_FILE_W)
        if os.path.exists(DARK_FILE_S): os.remove(DARK_FILE_S)
        if os.path.exists(DARK_FILE_W): os.remove(DARK_FILE_W)
        execfile("exp_setup.py")
        return

    def update2DView(self,widget):

        if self.fn2d==None or self.fn2d=='': return

        self.showMask=self.btnMaskOn.get_active()
        self.showRing=self.btnRingOn.get_active()
        self.SAXScLo=int(self.etSAXScLo.get_text())
        self.SAXScHi=int(self.etSAXScHi.get_text())
        self.WAXScLo=int(self.etWAXScLo.get_text())
        self.WAXScHi=int(self.etWAXScHi.get_text())

        self.plot2D()
        return

    def refreshFileList(self,widget):
        dirList=os.listdir(".")
        dirList.sort()
        dataList=[]
        self.fileList.clear()
        # only keep the file name that has both _SAXS and _WAXS endings 
        for fn in dirList:
            fn1 = fn.rpartition("_SAXS")
            if fn1[0]!='' and fn1[2]=='' and fn1[0].count("dark")==0:
                if fn1[0]+"_WAXS" in dirList:
                    dataList.append(fn1[0])
                    cur=self.fileList.append()
                    self.fileList.set(cur,0,False,1,False,2,fn1[0])
        
        # now set toggle bottons sensitive
        return

    def on_FileSel_changed(self,widget):
        sel=self.fileTv.get_cursor()[0]
        self.fn2d=self.fileList[sel][2]
        self.update2DView(None)
        return

    def sample_file_toggled(self,widget,Nrow):
        v=not widget.get_active()
        r_iter=self.fileList.iter_nth_child(None,int(Nrow))
        self.fileList.set_value(r_iter, 0, v)
        if self.fileList.get_value(r_iter, 1)==v and v==True:
            self.fileList.set_value(r_iter, 1, False)
        return

    def buffer_file_toggled(self,widget,Nrow):
        v=not widget.get_active()
        r_iter=self.fileList.iter_nth_child(None,int(Nrow))
        self.fileList.set_value(r_iter, 1, v)
        if self.fileList.get_value(r_iter, 0)==v and v==True:
            self.fileList.set_value(r_iter, 0, False)
        return

    def process_data(self,widget):
        outFn=self.etOutFn.get_text()
        sample_files=[]
        buffer_files=[]
        r_iter=self.fileList.get_iter_first()
        while (r_iter!=None):
            if self.fileList.get_value(r_iter, 0)==True:
                sample_files.append(self.fileList.get_value(r_iter, 2))
            if self.fileList.get_value(r_iter, 1)==True:
                buffer_files.append(self.fileList.get_value(r_iter, 2))
            r_iter=self.fileList.iter_next(r_iter)
        #print sample_files, buffer_files, outFn

        d1 = proc_SWAXS(sample_files,
                        buffer_files,
                        sdark,wdark,
                        qmax=0.16,qmin=0.115,reft=2500,
                        conc=0,
                        save1d=True,
			#saxsflat=sflat,
                        waxsflat=wflat,
                        fix_scale=-36.95
                        )
        d1.save(outFn)

        analyze(d1,qstart=0.02,qend=0.09,fix_qe=True,qcutoff=0.9,dmax=100)
        plt.show()

        return

    def quit(self,para):
        gtk.main_quit()

    def plot2D(self):

        if self.fn2d==None: return
        saxs_file = self.fn2d + '_SAXS'
        waxs_file = self.fn2d + '_WAXS'
    
        if self.fig2d==None:
            self.fig2d=plt.figure(figsize=(10,6))
            plt.show()
        #plt.set_current_figure(self.fig2d)
        plt.ioff()
        #plt.clf()

        if os.path.isfile(saxs_file):
            if self.ax1==None: 
                plt.subplot(121)
                #plt.subplots_adjust(bottom=0.1,top=0.9)
                self.ax1 = plt.gca()
            else:
                plt.sca(self.ax1)
            plt.cla()
            self.dsaxs = Data2d(saxs_file)
            self.dsaxs.set_exp_para(es)
            self.dsaxs.subtract(sdark.d2data)
            self.pax1 = Axes2dplot(self.ax1,self.dsaxs)
            if self.showMask:
                self.pax1.plot(mask=sdark.mask)
            else:
                self.pax1.plot()
    
            #self.pax1.set_color_scale(cm.gray)
            self.pax1.set_color_scale(cm.gist_yarg)
            cm.gist_yarg.set_under(color='b')
            cm.gist_yarg.set_over(color='r')
            self.pax1.img.set_clim(self.SAXScLo,self.SAXScHi)

            if (self.showRing):
                self.pax1.add_dec("P 0. 0. k+")
                self.pax1.add_dec("Q 0.1076 72 k--")
                self.pax1.add_dec("Q 0.2152 72 k--")
                self.pax1.add_dec("Q 0.3228 72 k--")
                self.pax1.add_dec("Q 0.4304 72 k--")

        if os.path.isfile(waxs_file):
            if self.ax2==None: 
                plt.subplot(122)
                #plt.subplots_adjust(bottom=0.1,top=0.9)
                self.ax2 = plt.gca()
            else:
                plt.sca(self.ax2)
            plt.cla()
            self.dwaxs = Data2d(waxs_file)
            self.dwaxs.set_exp_para(ew)
            self.dwaxs.subtract(wdark.d2data)
            self.pax2 = Axes2dplot(self.ax2,self.dwaxs)
        
            if self.showMask:
                self.pax2.plot(mask=wdark.mask)
            else:
                self.pax2.plot()
            
            #self.pax2.set_color_scale(cm.gray)
            self.pax2.set_color_scale(cm.gist_yarg)
            self.pax2.img.set_clim(self.WAXScLo,self.WAXScHi)
            
            if (self.showRing):
                self.pax2.add_dec("P 0. 0. k+")
                self.pax2.add_dec("Q 0.1076 72 k--")
                self.pax2.add_dec("Q 0.2152 72 k--")
                self.pax2.add_dec("Q 0.3228 72 k--")
                self.pax2.add_dec("Q 0.4304 72 k--")
                self.pax2.add_dec("Q 0.5380 72 k--")
                self.pax2.add_dec("Q 0.6456 72 k--")
                self.pax2.add_dec("Q 0.7532 72 k--")
                self.pax2.add_dec("Q 0.8608 72 k--")
                self.pax2.add_dec("Q 0.9684 72 k--")
                self.pax2.add_dec("Q 1.0760 72 k--")
                self.pax2.add_dec("Q 1.37 72 k--")

        plt.draw()
        plt.ion()

        return     

app = app()
gtk.main()
