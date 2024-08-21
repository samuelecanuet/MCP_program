import numpy as np
from matplotlib import pyplot as plt
from ROOT import *
import seaborn as sns
import subprocess
import matplotlib.colors as mcolors
import argparse


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

def DisplayTH1D(Hist, ax, color=None, label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, rebin=None, normalized=None, scaler=None):
        if rebin   != None: Hist.Rebin(rebin)
        if normalized == True: integral = Hist.Integral()
        else : integral = 1.
        if scaler != None: integral = scaler
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
        ax.step(bin_centers_x, hist_data, label=label)#", edgecolor='black', color='white', linewidth=0.5, width=bin_centers_x[0]-bin_centers_x[1])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        return [bin_centers_x, hist_data]
    

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

def ParcerParameter(B=None, V=None, Anode=0, Grid=0, runs=[]):
    print("\n########## Experimental Parameter #########")
    if B == False:
        print(f"## B = VARIABLE \t\t\t")
    else:
        print(f"## B = {B} \t\t\t\t")
    if V == False:
        print(f"## MCP : VARIABLE \t\t\t")
    else:
        print(f"## MCP : {V}V \t\t\t\t ")
    print(f"## Anode : {Anode}V \t\t\t")
    print(f"## Grid : {Grid}V \t\t\t\t")
    print(f"## Runs : {str(runs)} \t")
    print("###########################################\n")

root = "../../../../../../../mnt/hgfs/shared-2/data_2023-12_test_MCP/data/"

### ARGUMENTS 
parser = argparse.ArgumentParser(description="MCP reconstruction script",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-V", "--voltage", action="store", help="Voltage Value", default=False)
parser.add_argument("-B", "--field", action="store", help="Field Value")
parser.add_argument("-type", "--type", action="store", help="Field Value")
args = parser.parse_args()

######################################################
############ Q VALUE VOLTAGE ON BG / B ## ############
######################################################


dic_files_BG = {}
dic_files_BG["0T"] = [
    "MCP3_1700V_100V_BG_0017_0001",
    "MCP3_1800V_100V_BG_0016_0001",
    "MCP3_1900V_100V_BG_0017_0001",
    "MCP3_2000V_100V_BG_0012_0001",
    "MCP3_2200V_100V_BG_0013_0001",
    "MCP3_2500V_100V_BG_0014_0001",
    "MCP3_2600V_100V_BG_0018_0001",
    "MCP3_2400V_100V_BG_0015_0001",
    "MCP3_3000V_100V_Grid40V_BG_0019_0001"
]
dic_files_BG["2T"] = [
    "MCP3_2500V_100V_Grid45V_BG_2T_0025_0001",
    "MCP3_2800V_100V_Grid45V_BG_2T_0029_0001",
    "MCP3_3000V_100V_Grid45V_BG_2T_0030_0001"
]
dic_files_BG["3T"] = [
    "MCP3_3000V_100V_Grid45V_BG_3T_0032_0001"
]

dic_files_BEAM = {}
dic_files_BEAM["0T"] = [
    "MCP3_1800V_100V_Grid40V_WithBeam_0004_0001",
    "MCP3_2000V_100V_Grid45V_WithBeam_0T_0020_0001"
]
dic_files_BEAM["1T"] = [
    "MCP3_1800V_100V_Grid40V_WithBeam_1T_0007_0001",
    "MCP3_2000V_100V_Grid45V_WithBeam_1T_0021_0001",
    "MCP3_2500V_100V_Grid40V_WithBeam_1T_0009_0001",
]
dic_files_BEAM["2T"] = [
    "MCP3_2500V_100V_Grid45V_WithBeam_2T_0024_0001",
    "MCP3_2800V_100V_Grid45V_WithBeam_2T_0027_0001",
]
dic_files_BEAM["3T"] = [
    "MCP3_3000V_100V_Grid45V_WithBeam_3T_0033_0001"
]
dic_files_BEAM["4T"] = [
    "MCP3_3000V_100V_Grid45V_WithBeam_4T_0034_0001",
    "MCP3_3250V_100V_Grid45V_WithBeam_4T_0036_0001",
]
dic_files = {}
dic_files["BG"] = dic_files_BG
dic_files["BEAM"] = dic_files_BEAM

list_file = []
while args.type != "BEAM" and args.type != "BG" and args.type != "COMPARE":
    args.type = input("Specify the type of data : 'BG' or 'BEAM' (or COMPARE): ")

if args.type == "COMPARE":
    types = ["BEAM", "BG"]
else:
    types = [args.type]

for type in types:
    if args.field:
        LABEL = "Voltage (V)"
        list_file_field = dic_files[type][args.field]
        if args.voltage:
            for file in list_file_field:
                if args.voltage in file:
                    list_file.append(file)
        else:
            list_file = list_file_field
    else:
        if args.voltage:
            LABEL = "B (T)"     
            for key, dic in dic_files[type].items():
                for file in dic:
                    if args.voltage in file:
                        list_file.append(file)
        else:
            print("Give a Voltage or a Field Value")
            exit(0)

if len(list_file) == 0:
    print("No file has these conditions")
    exit(0)
maximum = 1
max_list = []
mean_list = []
std_list = []
max_bin_list = []
voltage_list = []
max_nullbin_list=[]
runs = []
fig, ax = plt.subplots(figsize=(8,5))

print("")
for name, color in zip(list_file, list(mcolors.TABLEAU_COLORS)):
    name = "Save-03-12-2023-193926_Grid_on_MCPBeamPos_HV2000V"
    input_file = root + name
    #runs.append(int(name[len(name)-9:len(name)-5]))
    print(f" {bcolors.OKGREEN} Running  : fast2root {bcolors.ENDC} for {name}")
    subprocess.run(["fast2root", input_file], stdout=subprocess.DEVNULL)
    root_file = TFile(input_file+".root")
    
    histo = root_file.Get("h_MCP_SIG_q1")
    if "3000" in name:
        histo.GetXaxis().SetRangeUser(5000, 80000)
    max_list.append(histo.GetMaximum())
    mean_list.append(histo.GetMean())
    std_list.append(histo.GetRMS())
    max_bin_list.append(histo.GetMaximumBin())
    if histo.GetMaximum() > maximum:
         maximum = histo.GetMaximum()
    for i in range(histo.GetNbinsX(), 0, -1):
        if histo.GetBinContent(i) > histo.GetMaximum()*1e-5:
            max_nullbin_list.append(histo.GetBinCenter(i))
            break
    
print("")
ParcerParameter(B=args.field, V=args.voltage, runs=runs)

##renormalization
for name, maxi, mean, max_bin in zip(list_file, max_list, mean_list, max_bin_list):
    input_file = root + name
    root_file = TFile(input_file+".root")
    if args.type == "COMPARE":
        if "BG" in name:
            label = "Background"
        else:
            label = "Beam"
    else:
        if args.voltage:
            try:
                label = name[name.index("T_")-1:name.index("T_")+1]
            except ValueError:
                label = "0T"
        else:
            label = name[name.index("V_")-4:name.index("V_")+1]
            
            
        voltage_list.append(int(label[:-1]))
    histo, data = DisplayTH1D(root_file.Get("h_MCP_SIG_q1"), ax, rebin = 1, label = label, title="", scaler = maxi/maximum)
    

#HITOGRAMS
ax.set_title(f"MCP Histogram Q value on {args.type}")
ax.set_xlabel("Channels")
ax.set_ylabel("Counts")
ax.set_yticks([])
ax.set_xlim(0, max(max_nullbin_list))
ax.legend()

if (len(list_file) != 1 and args.type != "COMPARE"):
    fig, ax1 = plt.subplots(figsize=(8,5))
    #SCATTER
    ax1.errorbar(voltage_list, mean_list, yerr=std_list, label="mean", fmt="o", capsize=5, color="black")
    ax1.set_title(f"MCP Mean Q value on {args.type}")
    ax1.set_xlabel(LABEL)
    ax1.set_ylabel("Counts")
    ax1.set_yticks([])
plt.show()
