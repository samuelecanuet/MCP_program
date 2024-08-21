import matplotlib.pyplot as plt
import numpy as np


with open ('fit_coordinate.txt', 'r') as f:
    for line in f:
        x, y = line.split(" ")
        print(x, y)
        plt.scatter(float(x), float(y), color='red')

plt.show()