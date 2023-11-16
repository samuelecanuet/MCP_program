# MCP Reconstruction Script

This Python script is designed for the reconstruction of MCP (Micro Channel Plate) data. It includes functionalities for calibration, data visualization, and fitting.

## Prerequisites

Make sure you have the required Python libraries installed:

```bash
pip install -r requirement.txt
```

```bash
make
```

Additionally, you need to have ROOT installed.
## Usage

To run the script, use the following command:

```bash
python script_name.py -f <input_filename> [-c]
```

 
- `-f` or `--filename`: Specify the input file (either .root or .fast format). 
- `-c` or `--calibration`: Enable calibration mode.
## Script Overview
### Libraries Used 
- `numpy` 
- `matplotlib` 
- `seaborn` 
- `ROOT` 
- `subprocess` 
- `argparse` 
- `ctypes`


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
- The reconstruction is constrained within the phase space (xdata,ydatax_{\text{data}}, y_{\text{data}}xdata​,ydata​) based on the selected points in Calibration Mode.


## Version 
### v1.0 : 
- GUI for calibration and reconstruction
- Manual calibration with corner point
- Polynomial fit
- Gaussian 1D and 2D fit

## Disclaimer

**Note:** The use of this software/project comes with absolutely no guarantees or warranties. The authors and contributors are not responsible for any consequences or damages that may arise from the use of this software. Use it at your own risk.

This project is provided "as is," and it is your responsibility to assess its fitness for your intended purpose. The authors make no claims regarding its reliability, completeness, or suitability for any particular use. Contributions and suggestions are welcome but are not guaranteed to be incorporated into the project.