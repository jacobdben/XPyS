# xps
Python tool for xps peak fitting

### Auxiliary:
- **gauss(x, C, mu, sigma)**, Gaussian function with coefficient C
- **s_orbital(x, C, mu, sig)**, s-orbital constructed from one Gaussian
- **split_p_orbital(x, C, mu1, mu2, sig)**, split p-orbital (p1/3, p2/3) constructed from two Gaussians
- **split_d_orbital(x, C, mu1, mu2, sig)**, split d-orbital (d2/5, d3/5) constructed from two Gaussians
- **linear_background(x, a, b)**, xps linear background

### XPS-class:
- **__init__(self, peaks, background, x_data, y_data)**, setup for guess of inital peak position and size, data to fit, and background.
- **analyticPeaks(self, x, *args)**, function used to curve fit the peaks. Adds peaks to fit according to information in self.peaks.
- **peakfit(self)**, does the peak fitting. Returns popt and pcov so the user may make their own plots if they wish.
- **plotPeakfit(self)**, plots the peak fitting data.
