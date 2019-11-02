# xps
Python tool for xps peak fitting

### Auxiliary:
- **gauss(x, C, mu, sigma)**, Gaussian function with coefficient C

### s_orb-class:
- **s_orbital(self, x, C, mu, sig)**, s-orbital constructed from one Gaussian

### p_orb-class:
- **dsep(self)**, get doublet separation
- **split_p_orbital(self, x, C, mu1, mu2, sig)**, split p-orbital (p1/3, p2/3) constructed from two Gaussians

### d_orb-class:
- **dsep(self)**, get doublet separation
- **split_d_orbital(self, x, C, mu1, mu2, sig)**, split d-orbital (d2/5, d3/5) constructed from two Gaussians


### XPS-class:
- **__init__(self, peaks, background, x_data, y_data)**, setup for guess of inital peak position and size, data to fit, and background.
- **analyticPeaks(self, x, *args)**, function used to curve fit the peaks. Adds peaks to fit according to information in self.peaks.
- **peakfit(self)**, does the peak fitting. Returns popt and pcov so the user may make their own plots if they wish.
- **plotPeakfit(self)**, plots the peak fitting data.
- **shirley_passive(self)**, xps passive shirley background.
- **linear_background(self, x, a, b)**, xps linear background.


## Usage
An example of how to use the program can be found in example.py. This code can be copied, and with a few small changes, used to curve
fit other XPS-data. Firstly, change the variables *path* and *filename* to specify the desired peak to be fitted. The file must be of
the type .txt and contain only x- and y-values separated by whitespace.

```
path = "XPS_data/"
filename = "Rb3d.txt"
```

Change the *peaks* variable according to what peaks are expected to exist. This variable is a dictionary where the keys are the strings
's', 'p' and 'd' describing the orbital type. The values are lists of lists with the initial guesses for height and position for the
peaks, also containing the doublet separation (in eV) in the case of p and d orbitals. In example.py we have looked for two separated 
d-orbitals:

```
peaks = {'d': [[111, 1.0e+5, 1.48],[114, 1.0e+4, 1.48]]}
```

Resulting peak fit:
<br>
![Peakfit](https://github.com/jacobdben/XPyS/blob/master/xpsfitted.png)
