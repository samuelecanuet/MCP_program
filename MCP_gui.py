import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Rectangle, Wedge, Polygon
from matplotlib.widgets import Button, Slider, RadioButtons
import csv
from ROOT import *
import subprocess
import argparse
import ctypes
import seaborn as sns
from matplotlib.backend_bases import MouseButton

custom_params = {
        "xtick.direction" : "out",
        "ytick.direction" : "out",
        "lines.markeredgecolor" : "k",
        "lines.markeredgewidth" : 0.5,
        "lines.linewidth" : 1,
        "lines.markersize" : 5,
        "figure.figsize" : (16,9),
        "font.family" : "serif",
        "ytick.labelsize" : 15,
        "xtick.labelsize" : 15,
        "axes.labelsize" : 20,
        "axes.titlesize" : 20,
        "legend.fontsize" : 15,
        "text.usetex" : True,
        # 'figure.subplot.left': 0.20, 
        # 'figure.subplot.bottom': 0.15, 
        # 'figure.subplot.right': 0.95, 
        # 'figure.subplot.top': 0.90
        }
#sns.set_theme(style = "ticks", rc=custom_params)

class MovePoint(object):
    def __init__(self, ax, graf, rec=None):
        self.ratio = 25
        self.ax = ax
        self.figcanvas = self.ax.figure.canvas
        self.graf = graf
        self.moved = None
        self.pointx = None
        self.pointy = None
        self.pressed = False
        self.start = False
        self.delete = False
        self.find = False
        self.initx = graf.get_center()[0]
        self.inity = graf.get_center()[1]
        self.label = graf.get_label()


        self.find = False
        self.graf.set_center((self.initx/self.ratio, self.inity/self.ratio))
        self.graf.set_radius(graf.get_radius()/self.ratio)
        self.graph.set_alpha(0.8)


        self.figcanvas.mpl_connect('button_press_event', self.mouse_press)
        self.figcanvas.mpl_connect('button_release_event', self.mouse_release)
        self.figcanvas.mpl_connect('motion_notify_event', self.mouse_move)

    def mouse_release(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return

        if self.pressed: 
            self.pressed = False
            self.start = False
            self.pointx = None
            self.pointy = None
            if self.graf.get_facecolor()==(1.0, 0.0, 0.0, 1):
                self.graf.set_facecolor('green')
            return
        
    def mouse_press(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if self.start: return
        if self.delete: return
        
        self.pointx = event.xdata
        self.pointy = event.ydata
        if self.graf.contains(event)[0]:
            if event.button is MouseButton.RIGHT:
                self.graf.remove() 
                self.delete = True
            self.pressed = True
        

    def mouse_move(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if not self.pressed: return
        self.start = True
        self.graf.set_center((event.xdata, event.ydata))

    


class MoveRec(object):
    def __init__(self, ax, rec):
        self.ratio = 20
        self.ax = ax
        self.figcanvas = self.ax.figure.canvas
        self.moved = None
        if str(rec)[:9] == "Rectangle":
            self.pointx = round(rec.get_xy()[0]-rec.get_height(), 5)
            self.pointy = round(rec.get_xy()[1]-rec.get_height(), 5)
        self.pressed = False
        self.start = False
        self.delete = False
        self.selected = False
        self.valid = False
        self.rec = rec
        self.label = rec.get_label()
        self.x = 0
        self.y = 0
                       
        self.figcanvas.mpl_connect('button_press_event', self.mouse_press)


    def mouse_press(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if self.start: return
        if self.rec.contains(event)[0] and event.inaxes == axs[0,1]:
            self.pressed = True
            if self.rec.get_facecolor() == (1.0, 0.0, 0.0, 1.0):
                self.rec.set_facecolor('orange')
                self.selected = True
                self.figcanvas.draw()
                return 
            elif self.rec.get_facecolor() == (1.0, 0.6470588235294118, 0.0, 1.0):
                self.rec.set_facecolor('green')
                self.valid = True
                self.select = False
                self.figcanvas.draw()
                return
            else:
                self.rec.set_facecolor('red')
                self.figcanvas.draw()
                self.selected = False
                self.valid = False
                return
            
    def GetCoordinates(self):
        return [self.pointx, self.pointy, self.x, self.y]
    
        
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

gROOT.ProcessLine("gROOT->SetBatch(kTRUE)")  

def gauss(x, amplitude, mean, sigma):
    return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))

def DisplayTGraph(Graph, ax, color='black', label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, marker=".", size=10):
    Graph = Graph.GetListOfPrimitives()[0]
    x, y = [], []
    for i in range(Graph.GetN()):
        x_val = ctypes.c_double()
        y_val = ctypes.c_double()
        Graph.GetPoint(i, x_val, y_val)
        x.append(float(x_val.value))
        y.append(float(y_val.value))

    
    if xlog   != None: ax.set_xscale('log')
    if ylog   != None: ax.set_yscale('log')

    sc = ax.scatter(x, y, color="black", marker = ".")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_title(title)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return [sc, ax]

def DisplayTH1D(Hist, ax, color=None, label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, rebin=None, normalized=None):
        if rebin   != None: Hist.Rebin(rebin)
        if normalized == True: integral = Hist.Integral()
        else : integral = 1.
        nbins_x = Hist.GetNbinsX()
        hist_data = np.zeros(nbins_x)
        bin_centers_x = np.zeros(nbins_x)

        for i in range(1, nbins_x + 1):
            hist_data[i - 1] = Hist.GetBinContent(i)/integral
            bin_centers_x[i - 1] = Hist.GetXaxis().GetBinCenter(i)

        

        if color    == None: color = "black"
        if label    == None: label = Hist.GetTitle()               
        if title    == None: title = Hist.GetTitle()
        if xlabel   == None: xlabel = Hist.GetXaxis().GetTitle()
        if ylabel   == None: ylabel = Hist.GetYaxis().GetTitle()
        if normalized==True: ylabel = "Normalized" + ylabel
        if xlim     == None: xlim = ( bin_centers_x.min(), bin_centers_x.max() )
        if ylim     == None and hist_data.max()*1.1 > ax.get_ylim()[1] : ylim = ( 0, hist_data.max()*1.1 )
        if xtick    != None: ax.set_xticks(np.linspace(xlim[0], xlim[1], xtick))
        if ytick    != None: ax.set_yticks(np.linspace(ylim[0], ylim[1], ytick))
        if xlog     != None: ax.set_xscale('log')
        if ylog     != None: 
            ax.set_yscale('log')
            ylim = ( 1, hist_data.max()*1.1 )
        
        # ax.bar(bin_centers_x, hist_data, label = label, color=color)
        ax.bar(bin_centers_x, hist_data, label="SRIM", edgecolor='black', color='white', linewidth=0.5, width=bin_centers_x[0]-bin_centers_x[1])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        return [bin_centers_x, hist_data]
    
def DisplayTH2D(Hist, ax, color='plasma', label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, zlog=None, rebinx=None, rebiny=None, vmax = None, vmin = None, view=None):
        if rebinx   != None: Hist.RebinX(rebinx)
        if rebiny   != None: Hist.RebinY(rebiny)

    
        nbins_x = Hist.GetNbinsX()
        nbins_y = Hist.GetNbinsY()

        hist_data = np.zeros((nbins_x, nbins_y))
        bin_centers_x = np.zeros(nbins_x)
        bin_centers_y = np.zeros(nbins_y)

        BD, BG, HG, HD = [0, 0], [0, 0], [0, 0], [0, 0]
        for i in range(1, nbins_x + 1):
            for j in range(1, nbins_y + 1):
                hist_data[i - 1, j - 1] = Hist.GetBinContent(i, j)
                bin_centers_x[i - 1] = Hist.GetXaxis().GetBinCenter(i)
                bin_centers_y[j - 1] = Hist.GetYaxis().GetBinCenter(j)

                if BD[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > BD[1] and Hist.GetXaxis().GetBinCenter(i) > 0. and Hist.GetYaxis().GetBinCenter(j) < 0.:
                    BD = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if BG[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > BG[1] and Hist.GetXaxis().GetBinCenter(i) < 0. and Hist.GetYaxis().GetBinCenter(j) < 0.:
                    BG = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if HG[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > HG[1] and Hist.GetXaxis().GetBinCenter(i) < 0. and Hist.GetYaxis().GetBinCenter(j) > 0.:
                    HG = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if HD[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > HD[1] and Hist.GetXaxis().GetBinCenter(i) > 0. and Hist.GetYaxis().GetBinCenter(j) > 0.:
                    HD = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]

        lim = max([i[1] for i in [BD, BG, HG, HD]])/np.sqrt(2)*1.1


        if label  == None: label = Hist.GetTitle()
        if title  == None: title = Hist.GetTitle()
        if xlabel == None: xlabel = Hist.GetXaxis().GetTitle()
        if ylabel == None: ylabel = Hist.GetYaxis().GetTitle()
        if xlim   == None: xlim = ( bin_centers_x.min(), bin_centers_x.max() )
        if ylim   == None: ylim = (bin_centers_y.min(), bin_centers_y.max())
        if xlim   == 'auto': xlim = (-lim, lim)
        if ylim   == 'auto': ylim = (-lim, lim)
        if xtick  != None: ax.set_xticks(np.linspace(xlim[0], xlim[1], xtick))
        if ytick  != None: ax.set_yticks(np.linspace(ylim[0], ylim[1], ytick))
        if xlog   != None: ax.set_xscale('log')
        if ylog   != None: ax.set_yscale('log')

        cax = ax.imshow(hist_data.T, extent=(bin_centers_x.min(), bin_centers_x.max(), bin_centers_y.min(), bin_centers_y.max()), origin='lower', aspect='auto', cmap=color)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if vmax      != None: cax.set_clim(vmax = vmax)
        if vmin      != None: cax.set_clim(vmin = vmin)

        return [cax, ax]

current=[]
current_patch=[]
sc = None

  
def writer(event):
    with open("coordinate_file.txt", "w") as file:
        writer = csv.writer(file, delimiter='\t')
        for square in liste:
            if square.x and square.y:
                values = square.GetCoordinates()
                writer.writerow([round(values[0], 5), round(values[1], 5), values[2], values[3]])

def reconstruction(event):
    if args.calibration:
        axs[0,2].clear()     
        writer(event)
        print(f" {bcolors.OKGREEN} Running  : fit {bcolors.ENDC}")
        subprocess.run(["fit", input_file, "config"])
        root_file = TFile(input_file)
        Histo = DisplayTH2D(root_file.Get("h_Image_corr"), axs[0,2], xlim='auto', ylim='auto', vmax = 20, title = "MCP Reconstructed")

        x, y = [], []
           
        with open("fit_coordinate.txt" , 'r') as fileout:
            for line in fileout:
                line = line.split(" ")
                x.append(float(line[0]))
                y.append(float(line[1]))
        #axs[0,2].scatter(x, y, color='green', s=10)
        #plt.show()
        return Histo, sc

    else:
        subprocess.run(["fit", input_file])
        root_file = TFile(input_file)
        Histo = DisplayTH2D(root_file.Get("h_Image_corr"), axs[0,2], xlim='auto', ylim='auto', title = "MCP reconstructed")
        plt.show()
        return Histo, None

def printer(nombre, error):
    if error == 0:
        return (nombre, error)
    else:
        str_er = list(str(error))
        if str_er[0] != "0": return (float(str(round(nombre))[0]), round(error))
        for i in range(len(str_er)):
            if str_er[i] != "0" and str_er[i] != ".":
                index = i
                break
        str_er.index(str_er[index])
        return float(str(round(nombre, index-1))), float(str(round(error, index-1)))

def One_fit(event):
        liste2=[]
        liste3=[]
    # try : 
        root_file = TFile(input_file)
        th1d = root_file.Get("h_Image") 
        th1d.GetXaxis().SetRangeUser(axs[0,0].get_xlim()[0], axs[0,0].get_xlim()[1])
        th1d.GetYaxis().SetRangeUser(axs[0,0].get_ylim()[0], axs[0,0].get_ylim()[1])
    
        ### Display and Do Hist-Fiiting
        axs[1,1].clear()
        th_x = th1d.ProjectionX("projectionX")
        DisplayTH1D(th_x, axs[1,1], title="X projection")
        gaussian_x = TF1("gaussian", "gaus", axs[0,0].get_xlim()[0], axs[0,0].get_xlim()[1])
        gaussian_x.SetParameters(th_x.GetMaximum(), 0, 1)
        th_x.Fit("gaussian", "R")
        x = np.linspace(axs[0,0].get_xlim()[0], axs[0,0].get_xlim()[1], 10000)
        axs[1,1].plot(x, gauss(x, gaussian_x.GetParameter(0), gaussian_x.GetParameter(1), gaussian_x.GetParameter(2)), color="red")

        axs[1,0].clear()
        th_y = th1d.ProjectionY("projectionY")
        DisplayTH1D(th_y, axs[1,0], title = "Y projection")
        gaussian_y = TF1("gaussian", "gaus", axs[0,0].get_ylim()[0], axs[0,0].get_ylim()[1])
        gaussian_y.SetParameters(th_y.GetMaximum(), 0, 1)
        th_y.Fit("gaussian", "R")
        y = np.linspace(axs[0,0].get_ylim()[0], axs[0,0].get_ylim()[1], 10000)
        axs[1,0].plot(y, gauss(y, gaussian_y.GetParameter(0), gaussian_y.GetParameter(1), gaussian_y.GetParameter(2)), color="red")

        for square in liste:
            if square.selected:
                square.rec.set_facecolor("green")
                square.valid = True
                square.selected = False
                square.x = gaussian_x.GetParameter(1)
                square.y = gaussian_y.GetParameter(1)
                break


        liste2.append(Circle( (gaussian_x.GetParameter(1), gaussian_y.GetParameter(1)), 0.02, facecolor="green", alpha=0.5))
        liste3.append(MoveRec(axs[0,0], liste2[-1]))
        axs[0,0].add_patch(liste2[-1])  




    # except AttributeError:
    #      print(f" {bcolors.WARNING} Please, reconstruct the image before fitting. {bcolors.ENDC}")


def gaussian2D(x, par):
    return par[0] * TMath.Gaus(x[0], par[1], par[2]) * TMath.Gaus(x[1], par[3], par[4])

def Two_fit(event):
    # try :
        root_file = TFile(input_file)
        th2d = root_file.Get("h_Image_corr")
        th2d.GetXaxis().SetRangeUser(axs[0, 1].get_xlim()[0], axs[0, 1].get_xlim()[1])
        th2d.GetYaxis().SetRangeUser(axs[0, 1].get_ylim()[0], axs[0, 1].get_ylim()[1])
        gaussian = TF2("gaussian_function", gaussian2D, 
                                axs[0, 1].get_xlim()[0], axs[0, 1].get_xlim()[1], 
                                axs[0, 1].get_ylim()[0], axs[0, 1].get_ylim()[1], 5)
        gaussian.SetParameters(1, 0, 1, 0, 1)
        th2d.Fit(gaussian, "R")

        ### Display Results
        fontsize = 15
        axs[1,2].clear()
        axs[1,2].axis("off")
        axs[1,2].text(0.0, 0.95, "Fit Paramters : ", fontsize = fontsize)
        axs[1,2].text(0.05, 0.7, "Mean = {} ± {} mm".format(printer(gaussian.GetParameter(1), gaussian.GetParError(1))[0], printer(gaussian.GetParameter(1), gaussian.GetParError(1))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.6, "Std = {} ± {} mm".format(printer(gaussian.GetParameter(2), gaussian.GetParError(2))[0], printer(gaussian.GetParameter(2), gaussian.GetParError(2))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.3, "Mean = {} ± {} mm".format(printer(gaussian.GetParameter(3), gaussian.GetParError(3))[0], printer(gaussian.GetParameter(3), gaussian.GetParError(3))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.2, "Std = {} ± {} mm".format(printer(gaussian.GetParameter(4), gaussian.GetParError(4))[0], printer(gaussian.GetParameter(4), gaussian.GetParError(4))[1]), fontsize = fontsize)
        


    # except AttributeError:
    #      print(f" {bcolors.WARNING} Please, reconstruct the image before fitting. {bcolors.ENDC}")

def update_max(val):
    HIST[0].set_clim(vmax = val)
    if len(axs[0,2].get_images()) > 0: axs[0,2].get_images()[0].set_clim(vmax = val)

def update_min(val):
    HIST[0].set_clim(vmin = val)
    if len(axs[0,2].get_images()) > 0: axs[0,2].get_images()[0].set_clim(vmin = val)

if __name__ == '__main__': 

    ### ARGUMENTS 
    parser = argparse.ArgumentParser(description="MCP reconstruction script",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--filename", action="store", help="Input filename", default=False)
    parser.add_argument("-c", "--calibration", action="store_true", help="Calibration mode")
    args = parser.parse_args()
    config = vars(args)       

    ### INPUT MCP DATA  
    fig, ax = plt.subplots()
    
    if args.filename:
        input_file = args.filename

        if "root" in input_file:
            try : root_file = TFile(input_file); print(f"\n {bcolors.OKGREEN} Openning : {input_file} {bcolors.ENDC}")
            except OSError : print(f"\n {bcolors.FAIL} Wrong file name {bcolors.ENDC}\n"), exit(0)
            
            try : HIST = DisplayTH2D(root_file.Get("h_Image"), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : h_Image {bcolors.ENDC}"); name = True
            except AttributeError : name = None
           
            while name == None:
                name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}")
                try : HIST = DisplayTH2D(root_file.Get(name), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : {name} {bcolors.ENDC}")
                except AttributeError : name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}"); name = None
            


        elif "fast" in input_file:
            try : print(f"\n {bcolors.OKGREEN} Running  : fast2root {bcolors.ENDC}"); print(f" {bcolors.OKGREEN} Openning : {input_file} {bcolors.ENDC}"); resultat = subprocess.run(["fast2root", input_file[:-5]], stdout=subprocess.DEVNULL)
            except FileNotFoundError : print(f"\n {bcolors.FAIL} Wrong Executable Name File (Fast Reader Name group2tree_SL_MCP) {bcolors.ENDC}\n"), exit(0)

            input_file = input_file[:-5]+".root"
            try : root_file = TFile(input_file); print(f" {bcolors.OKGREEN} Openning : {input_file} {bcolors.ENDC}")
            except OSError : print(f"\n {bcolors.FAIL} Wrong file name {bcolors.ENDC}\n"), exit(0)
            
            try : HIST = DisplayTH2D(root_file.Get("h_Image"), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : h_Image {bcolors.ENDC}"); name = True
            except AttributeError : name = None
           
            while name == None:
                name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}")
                try : HIST = DisplayTH2D(root_file.Get(name), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : {name} {bcolors.ENDC}")
                except AttributeError : name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}"); name = None
            


        else:
            print(f"\n{bcolors.FAIL} Error : Please give a root or fast file as argument {bcolors.ENDC}\n")
            exit(0)
    else:
        print(f"\n{bcolors.FAIL} Error : Please give a root or fast file as argument {bcolors.ENDC}\n")
        exit(0)
    plt.close()    


    ### MODE SELECTION
    if args.calibration:
        ### CANVAS
        list_point = []
        global current_point
        for i in range(16*16):
            list_point.append([[None, None],[None, None]])


        fig, axs = plt.subplots(2, 3, figsize = (16, 8), gridspec_kw=dict(height_ratios=(2, 1)))
        axs[0,0].set_title("MCP data")
        axs[0,1].set_title("MCP sample")
        axs[0,1].set_xlabel("X (mm)")
        axs[0,1].set_xlabel("Y (mm)")
        axs[0,1].set_xlim(-10, 10)
        axs[0,1].set_ylim(-10, 10)
        axs[0,2].set_title("MCP reconstructed")
        axs[0,2].set_xlabel("X (mm)")
        axs[0,0].set_xlabel("Y (mm)")
        axs[1,0].axis("off")
        axs[1,1].axis("off")
        axs[1,2].remove()

        HIST = DisplayTH2D(root_file.Get("h_Image"), axs[0,0], title = "MCP data", xlim='auto' , ylim='auto')

        ### CONSTRUCT MCP SKETCH
        background = Rectangle((-10, -10), 20, 20, fc='black', alpha = 0.6)
        axs[0,1].add_patch(background)

        pitch = 2
        size = 1.2
        radius = 0.15
        MCP = 7.5

        liste = []
        liste1=[]
        for x in range(0,8):
            for y in range(0,8):

                allpoint=0
                x_corr = x*pitch-4*pitch+(pitch-size)/2
                y_corr = y*pitch-4*pitch+(pitch-size)/2
        
                liste1.append(Rectangle( (round(x_corr, 5), round(y_corr, 5)), size, size, facecolor="red"))
                liste.append(MoveRec(axs[0,1], liste1[-1]))
                axs[0,1].add_patch(liste1[-1])    

        wedge = Wedge((0,0), 15, 0, 360, width=7.5, color="black", alpha=0.5)
        axs[0,1].add_patch(wedge)

        ### CONSTRUCT peaks
        # radius = 0.02
        # liste2 = []
        # liste3 = []
        # with open("coordinate_file_back.txt", "r") as file:
        #     for line in file:
        #         line = line.split("\t")

        #         liste2.append(Circle( (float(line[0]), float(line[1])), radius, facecolor="red"))
        #         liste3.append(MoveRec(axs[0,0], liste2[-1]))
        #         axs[0,0].add_patch(liste2[-1])      


        ### WRITER BUTTON
        button_pos = plt.axes([0.72, 0.14, 0.15, 0.05])
        button_w = Button(button_pos, 'Write coordinates file')
        button_w.on_clicked(writer)

        ### RECONSTRUCT BUTTON
        button_reconstruction_pos = plt.axes([0.72, 0.08, 0.15, 0.05])
        button_reconstruction = Button(button_reconstruction_pos, 'Reconstruction')
        button_reconstruction.on_clicked(reconstruction)

        ### VMAX SLIDER
        axfreq_max = fig.add_axes([0.74, 0.24, 0.15, 0.01])
        freq_slider_max = Slider(
        ax=axfreq_max,
        label='Max Value',
        valmin=1,
        valmax=500,
        valinit=10)   
        freq_slider_max.on_changed(update_max)

        ### VMIN SLIDER
        axfreq_min = fig.add_axes([0.74, 0.29, 0.15, 0.01])
        freq_slider_min = Slider(
        ax=axfreq_min,
        label='Min Value',
        valmin=1,
        valmax=500,
        valinit=1)   
        freq_slider_min.on_changed(update_min)

        ###HIDER
        def hider(label):
            if label == 'hide':
                HIST[0].set_alpha(0)
            else:
                HIST[0].set_alpha(1)
            fig.canvas.draw()

        radio = RadioButtons(fig.add_axes([0.875, 0.08, 0.05, 0.11]), ('show', 'hide'))
        radio.on_clicked(hider)

        ### 1D fit BUTTON
        button_fit1d_pos = plt.axes([0.665, 0.08, 0.05, 0.11])
        button_fit1d = Button(button_fit1d_pos, '1D fit')
        button_fit1d.on_clicked(One_fit)
        

        plt.subplots_adjust(bottom=0.03, top=0.95, left=0.05, right=0.95)
        plt.show()

    else:
        fig, axs = plt.subplots(2, 3, figsize = (14, 9))
        axs[0,0].set_title("MCP data")
        HIST = DisplayTH2D(root_file.Get("h_Image"), axs[0,0], title = "MCP data", xlim='auto' , ylim='auto')
        axs[0,1].set_title("MCP reconstructed")
        axs[0,1].set_xlabel("X (mm)")
        axs[0,1].set_xlabel("Y (mm)")
        axs[0,2].set_title("MCP fit")
        axs[1,0].remove()
        axs[1,2].axis("off")

        ### RECONSTRUCT BUTTON
        button_reconstruction_pos = plt.axes([0.135, 0.25, 0.2, 0.075])
        button_reconstruction = Button(button_reconstruction_pos, 'Reconstruction')
        button_reconstruction.on_clicked(reconstruction)

        ### 1D fit BUTTON
        button_fit1d_pos = plt.axes([0.135, 0.15, 0.09, 0.075])
        button_fit1d = Button(button_fit1d_pos, '1D fit')
        button_fit1d.on_clicked(One_fit)

        ### 2D fit BUTTON
        button_fit2d_pos = plt.axes([0.245, 0.15, 0.09, 0.075])
        button_fit2d = Button(button_fit2d_pos, '2D fit')
        button_fit2d.on_clicked(Two_fit)

        plt.subplots_adjust(hspace=0.3, wspace=0.3)
        plt.show()
