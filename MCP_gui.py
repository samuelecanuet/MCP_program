import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Rectangle, Wedge
from matplotlib.widgets import Button, Slider
import csv
from ROOT import *
import subprocess
import argparse
import ctypes
import seaborn as sns

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
sns.set_theme(style = "ticks", rc=custom_params)

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
    ax.set_title(title)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return [sc, ax]

def DisplayTH1D(Hist, ax, color=None, label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, rebin=None, normalized=None):
        if rebin   != None: Hist.Rebin(rebin)
        if normalized == True: integral = Hist.Integral()
        else : integral = 1.

        print(Hist)
        
        nbins_x = Hist.GetNbinsX()
        print(nbins_x)

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
def mouse_event(event):
    try:
        if event.inaxes.get_title() == "MCP sample":
            for i, patches in enumerate(event.inaxes.patches):
                if patches.contains(event)[0] and patches.get_label() != "":
                    
                    if patches.get_facecolor() == (1.0, 0.6470588235294118, 0.0, 1.0) or patches.get_facecolor() == (0.0, 0.5019607843137255, 0.0, 1.0):
                        patches.set_facecolor("red")
                        list_point[int(patches.get_label())][1][2].remove()
                        event.canvas.draw()
                        list_point[int(patches.get_label())][1] = [None, None, None]
                        
                    else: 
                        patches.set_facecolor("orange")
                        event.canvas.draw()
                        current.append(int(patches.get_label()))  
                        current_patch.append(i)  
                            
                    
        if event.inaxes.get_title() == "MCP data":   
            list_point[current[-1]][1] = [event.xdata, event.ydata, event.inaxes.scatter(event.xdata, event.ydata, color="white", marker = '+', s=50)]
            event.canvas.figure.axes[1].patches[current_patch[-1]].set_facecolor("green")
            event.canvas.draw()
    except AttributeError:
        None

def on_move(event):
    global sc
    if event.inaxes is axs[0,0]:
        x, y = event.xdata, event.ydata
        del axs[1,0].get_children()[0]
        axs[1,0].set_xlim(x - 0.15, x + 0.15)
        axs[1,0].set_ylim(y - 0.15, y + 0.15)
        
        if sc: sc.remove()
        sc = axs[1,0].scatter(x, y, marker = "+", s = 125, color="red")
        
        event.canvas.draw()
   
def writer(event):
    with open("coordinate_file.txt", "w") as file:
        writer = csv.writer(file, delimiter='\t')
        for values in list_point:
            if values[1][0] == None or values[1][1] == None:
                continue
            writer.writerow([values[0][0], values[0][1], values[1][0], values[1][1]])

def reconstruction(event):
    if args.calibration:
        writer(event)
        print(f" {bcolors.OKGREEN} Running  : fit {bcolors.ENDC}")
        subprocess.run(["fit", input_file, "config"])
        root_file = TFile(input_file)
        Histo = DisplayTH2D(root_file.Get("h_Image_corr"), axs[0,2], xlim='auto', ylim='auto', vmax = 20, title = "MCP Reconstructed")
        sc = DisplayTGraph(root_file.Get("Residus"), axs[1,2])
        plt.show()
        return Histo, sc

    else:
        subprocess.run(["fit", input_file])
        root_file = TFile(input_file)
        Histo = DisplayTH2D(root_file.Get("h_Image_corr"), axs[0,1], xlim='auto', ylim='auto', title = "MCP reconstructed")
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
    try : 
        root_file = TFile(input_file)
        th1d = root_file.Get("h_Image_corr") 
        th1d.GetXaxis().SetRangeUser(axs[0,1].get_xlim()[0], axs[0,1].get_xlim()[1])
        th1d.GetYaxis().SetRangeUser(axs[0,1].get_ylim()[0], axs[0,1].get_ylim()[1])
        
        ### Display and Do Hist-Fiiting
        axs[1,1].clear()
        th_x = th1d.ProjectionX("projectionX")
        DisplayTH1D(th_x, axs[1,1], title="X projection")
        gaussian_x = TF1("gaussian", "gaus", axs[0,1].get_xlim()[0], axs[0,1].get_xlim()[1])
        gaussian_x.SetParameters(th_x.GetMaximum(), 0, 1)
        th_x.Fit("gaussian", "R")
        x = np.linspace(axs[0,1].get_xlim()[0], axs[0,1].get_xlim()[1], 10000)
        axs[1,1].plot(x, gauss(x, gaussian_x.GetParameter(0), gaussian_x.GetParameter(1), gaussian_x.GetParameter(2)), color="red")

        axs[0,2].clear()
        th_y = th1d.ProjectionY("projectionY")
        DisplayTH1D(th_y, axs[0,2], title = "Y projection")
        gaussian_y = TF1("gaussian", "gaus", axs[0,1].get_ylim()[0], axs[0,1].get_ylim()[1])
        gaussian_y.SetParameters(th_y.GetMaximum(), 0, 1)
        th_y.Fit("gaussian", "R")
        y = np.linspace(axs[0,1].get_ylim()[0], axs[0,1].get_ylim()[1], 10000)
        axs[0,2].plot(y, gauss(y, gaussian_y.GetParameter(0), gaussian_y.GetParameter(1), gaussian_y.GetParameter(2)), color="red")

        ### Display Results
        fontsize = 15
        axs[1,2].clear()
        axs[1,2].axis("off")
        axs[1,2].text(0.0, 0.95, "Fit Paramters : ", fontsize = fontsize)
        axs[1,2].text(0.05, 0.8, "X Projection : ", fontsize = fontsize)
        axs[1,2].text(0.05, 0.7, "Mean = {} +/- {} mm".format(printer(gaussian_x.GetParameter(1), gaussian_x.GetParError(1))[0], printer(gaussian_x.GetParameter(1), gaussian_x.GetParError(1))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.6, "Std = {} +/- {} mm".format(printer(gaussian_x.GetParameter(2), gaussian_x.GetParError(2))[0], printer(gaussian_x.GetParameter(2), gaussian_x.GetParError(2))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.4, "Y Projection : ", fontsize = fontsize)
        axs[1,2].text(0.05, 0.3, "Mean = {} +/- {} mm".format(printer(gaussian_y.GetParameter(1), gaussian_y.GetParError(1))[0], printer(gaussian_y.GetParameter(1), gaussian_y.GetParError(1))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.2, "Std = {} +/- {} mm".format(printer(gaussian_y.GetParameter(2), gaussian_y.GetParError(2))[0], printer(gaussian_y.GetParameter(2), gaussian_y.GetParError(2))[1]), fontsize = fontsize)
        

    except AttributeError:
         print(f" {bcolors.WARNING} Please, reconstruct the image before fitting. {bcolors.ENDC}")


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
        axs[1,2].text(0.05, 0.7, "Mean = {} +/- {} mm".format(printer(gaussian.GetParameter(1), gaussian.GetParError(1))[0], printer(gaussian.GetParameter(1), gaussian.GetParError(1))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.6, "Std = {} +/- {} mm".format(printer(gaussian.GetParameter(2), gaussian.GetParError(2))[0], printer(gaussian.GetParameter(2), gaussian.GetParError(2))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.3, "Mean = {} +/- {} mm".format(printer(gaussian.GetParameter(3), gaussian.GetParError(3))[0], printer(gaussian.GetParameter(3), gaussian.GetParError(3))[1]), fontsize = fontsize)
        axs[1,2].text(0.05, 0.2, "Std = {} +/- {} mm".format(printer(gaussian.GetParameter(4), gaussian.GetParError(4))[0], printer(gaussian.GetParameter(4), gaussian.GetParError(4))[1]), fontsize = fontsize)
        


    # except AttributeError:
    #      print(f" {bcolors.WARNING} Please, reconstruct the image before fitting. {bcolors.ENDC}")

def update_max(val):
    HIST[0].set_clim(vmax = val)
    HIST_z[0].set_clim(vmax = val)
    if len(axs[0,2].get_images()) > 0: axs[0,2].get_images()[0].set_clim(vmax = val)

def update_min(val):
    HIST[0].set_clim(vmin = val)
    HIST_z[0].set_clim(vmin = val)
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


        fig, axs = plt.subplots(2, 3, figsize = (15, 9), gridspec_kw=dict(height_ratios=(2, 1)))
        axs[0,0].set_title("MCP data")
        axs[0,1].set_title("MCP sample")
        axs[0,1].set_xlabel("X (mm)")
        axs[0,1].set_xlabel("Y (mm)")
        axs[0,2].set_title("MCP reconstructed")
        axs[0,2].set_xlabel("X (mm)")
        axs[0,2].set_xlabel("Y (mm)")
        axs[1,1].remove()
        axs[1,2].set_xlabel("Distance from (0,0)")
        axs[1,2].set_ylabel("Distance from real point")

        fig.suptitle('MCP Selection Points')

        fig.canvas.mpl_connect('button_press_event', mouse_event)

        HIST = DisplayTH2D(root_file.Get("h_Image"), axs[0,0], title = "MCP data", xlim='auto' , ylim='auto')
        HIST_z = DisplayTH2D(root_file.Get("h_Image"), axs[1,0], title = "MCP data", xlim='auto' , ylim='auto')

        ### CONSTRUCT MCP SKETCH
        pitch = 2
        size = 1.2
        radius = 0.2
        background = Rectangle((-10, -10), 20, 20, fc='black', alpha = 0.6)
        axs[0,1].add_patch(background)

        for x in range(0,8):
            for y in range(0,8):
                x_corr = x*pitch-4*pitch+(pitch-size)/2
                y_corr = y*pitch-4*pitch+(pitch-size)/2
                rec = Rectangle((x_corr, y_corr), size, size, fc='blue')
                c1 = Circle((x_corr, y_corr), radius, fc='red', label=2*x+8*4*y)
                list_point[2*x+8*4*y][0] = [x_corr, y_corr]
                c2 = Circle((x_corr+size, y_corr+size), radius, fc='red', label=2*x+8*4*y+1)
                list_point[2*x+8*4*y+1][0] = [round(x_corr+size, 5), round(y_corr+size, 5)]
                c3 = Circle((x_corr+size, y_corr), radius, fc='red', label=2*x+8*4*y+2*8)
                list_point[2*x+8*4*y+2*8][0] = [round(x_corr+size, 5), y_corr]
                c4 = Circle((x_corr, y_corr+size), radius, fc='red', label=2*x+8*4*y+2*8+1)
                list_point[2*x+8*4*y+2*8+1][0] = [x_corr, round(y_corr+size, 5)]
                
                axs[0,1].add_patch(rec)
                axs[0,1].add_patch(c1)
                axs[0,1].add_patch(c2)
                axs[0,1].add_patch(c3)
                axs[0,1].add_patch(c4)

        wedge = Wedge((0,0), 15, 0, 360, width=7.5, color="black", alpha=0.5)
        axs[0,1].add_patch(wedge)
        axs[0,1].set_xlim(-10, 10)
        axs[0,1].set_ylim(-10, 10)

        ### WRITER BUTTON
        button_pos = plt.axes([0.4, 0.05, 0.2, 0.075])
        button_w = Button(button_pos, 'Write coordinates file')
        button_w.on_clicked(writer)

        ### RECONSTRUCT BUTTON
        button_reconstruction_pos = plt.axes([0.68, 0.05, 0.2, 0.075])
        button_reconstruction = Button(button_reconstruction_pos, 'Reconstruction')
        button_reconstruction.on_clicked(reconstruction)

        ### VMAX SLIDER
        axfreq_max = fig.add_axes([0.1, 0.05, 0.2, 0.01])
        freq_slider_max = Slider(
        ax=axfreq_max,
        label='Max Value',
        valmin=1,
        valmax=150,
        valinit=10)   
        freq_slider_max.on_changed(update_max)

        ### VMIN SLIDER
        axfreq_min = fig.add_axes([0.1, 0.1, 0.2, 0.01])
        freq_slider_min = Slider(
        ax=axfreq_min,
        label='Min Value',
        valmin=1,
        valmax=100,
        valinit=1)   
        freq_slider_min.on_changed(update_min)

        ### LENS
        fig.canvas.mpl_connect('motion_notify_event', on_move)
        axs[1,0].set_xticks([])
        axs[1,0].set_yticks([])
        axs[1,0].set_title("")
        

        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.85)
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
