# MCP Reconstruction Script

This Python script is designed for the reconstruction of MCP (Micro Channel Plate) data. It includes functionalities for calibration, data visualization, and fitting.

## Prerequisites

- **GSL** (https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz)
- **fasterac** (https://faster.in2p3.fr/index.php/download/category/2-software?download=20:fasterac-2-17-0-tar-gz)
- **ROOT**
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
python script_name.py -f <input_filename> [-c]
```

 
- `-f` or `--filename`: Specify the input file (either .root or .fast format). 
- `-c` or `--calibration`: Enable calibration mode.
## Script Overview

### 1. Calibration Mode
- Execute the script with the -c argument.
- Identify the red corners of the MCP grid sample and link them to the corresponding points on the image from the dataset.
- Record the coordinate correspondences in a file (sample/data) or reconstruct the image without deformation.

### 2. Beam Characterization Mode
- Run the script without the -c argument, using a pre-generated coordinate file.
- Reconstruct the image without deformation.
- Perform a fitting process and access the fitting parameters.

## Fitting Details
- Reconstruction is performed using two polynomial fits: x<sub>real</sub> = g(x<sub>data</sub>, y<sub>data</sub>) and y<sub>real</sub> = f(x<sub>data</sub>, y<sub>data</sub>). You can modify the degree of the function in the "fit.C" file. 
- The reconstruction is constrained within the phase space (x<sub>data</sub>,y<sub>data</sub>) based on the points selection in Calibration Mode.


## Version 
### v1.0 : 
- GUI for calibration and reconstruction
- Manual calibration with corner point
- Polynomial fit
- Gaussian 1D