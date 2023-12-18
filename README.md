# MCP Reconstruction Script

This Python script is designed for the reconstruction of MCP (Micro Channel Plate) data. It includes functionalities for calibration, data visualization, and fitting.

## Prerequisites

- **GSL** (https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz)
- **fasterac** (https://faster.in2p3.fr/index.php/download/category/2-software?download=20:fasterac-2-17-0-tar-gz)
- **ROOT** and **PyROOT**
- run commands :

```bash
./configure
```

```bash
source ~/.cshrc
```
## Usage

To run the script, use the following command:

```bash
python MCP_gui.py -f <input_filename> [-c]
```

 
- `-f` or `--filename`: Specify the input file (either .root or .fast format). 
- `-c` or `--calibration`: Enable calibration mode.
## Script Overview

### 1. Calibration Mode
- Execute the script with the -c argument.
- Identify the red square of the MCP grid sample and link them to the corresponding area on the image from the dataset.
- Record the coordinate correspondences in a file and/or gain match the image without deformation.

### 2. Beam Characterization Mode
- Run the script without the -c argument, using a pre-generated coordinate file.
- Reconstruct the image without deformation.
- Perform a fitting process and access the fitting parameters.

## Fitting Details
- The reconstruction is performed with the following formulas : 
```math x = -\frac{
\log(\frac{HD*BD}{HG*BG})}
{
\log(\frac{HD*BD*HG*BG}{TOT^4})}
\quad\quad\quad\quad\quad\quad\quad\quad
 y = -\frac{\log\left(\frac{HD \times BG}{HG \times BD}\right)}{\log\left(\frac{HD \times BD \times HG \times BG}{TOT^4}\right)}```

- Gain match is performed using two 3rd order polynomial fits: x<sub>real</sub> = g(x<sub>data</sub>, y<sub>data</sub>) and y<sub>real</sub> = f(x<sub>data</sub>, y<sub>data</sub>). You can modify the degree of the function in the "fit.C" file. 

### Data File Example
Background at 0T for Gain Match : *MCP3_2600V_100V_BG_0015_0001.fast*


## Version 
### v1.0 : 
- GUI for calibration and reconstruction
- Manual calibration with corner point
- Polynomial fit
- Gaussian 1D

### v2.0 : 
- GUI for calibration and log reconstruction
- Manual calibration with square area
- Polynomial fit
- Gaussian 1D for beam