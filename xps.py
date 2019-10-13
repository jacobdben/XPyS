from scipy.optimize import curve_fit
from numpy import trapz
import numpy as np
import matplotlib.pyplot as plt

def gauss(x, C, mu, sigma):
    return abs(C)/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

class s_orb:
    def s_orbital(self, x, C, mu, sig):
        return gauss(x, C, mu, sig)

class p_orb:
    def __init__(self, doublet_sep):
        self.dsp = doublet_sep

    def dsep(self):
        return self.dsp

    def split_p_orbital(self, x, C, mu1, sig):
        return gauss(x, C, mu1, sig) + gauss(x, 1 / 2 * C, mu1 + self.dsp, sig)

class d_orb:
    def __init__(self, doublet_sep):
        self.dsp = doublet_sep

    def dsep(self):
        return self.dsp

    def split_d_orbital(self, x, C, mu1, sig):
        return gauss(x, C, mu1, sig) + gauss(x, 2 / 3 * C, mu1 + self.dsp, sig)

class XPS:
    def __init__(self, peaks, background_type, x_data, y_data):
        """
        :param peaks: dict containing type of orbital with list of approx position and height of each peak.
        Has the form   {"s": [ [peak1pos, peak1height], [peak2pos, peak2height], ...],
                        "p": [ [peak1pos, peak1height, peak1sep], [peak2pos, peak2height, peak1sep], ...,
                        "d": [ [peak1pos, peak1height, peak1sep], [peak2pos, peak2height, peak1sep], ...}
        :param background: String containing type of background, e.g "linear" or "shirley"
        """
        self.peaks = peaks                          # Information about peaks
        self.bg_type = background_type              # String for choosing background type
        self.x = x_data                             # x-values to be curve-fitted
        self.y = y_data                             # y-values to be curve-fitted
        self.background = np.zeros(x_data.shape[0]) # Background for plotting
        self.s_orbitals = []                        # List of s-orb objects
        self.p_orbitals = []                        # List of p-orb objects
        self.d_orbitals = []                        # List of d-orb objects
        self.p0 = []                                # Initial parameter values for curve-fitting

        # Reverse x-axis if given x-data is low to high
        if self.x[0] < self.x[-1]:
            self.x = self.x[-1:0:-1]

        # For each s-orbital, get initial parameter values and create an s-orb object
        if 's' in peaks:
            for s in peaks.get('s'):
                self.p0.append(s[1])
                self.p0.append(s[0])
                self.p0.append(1.0)
                self.s_orbitals.append(s_orb())

        # For each p-orbital, get initial parameter values and create a p-orb object
        if 'p' in peaks:
            for p in peaks.get('p'):
                self.p0.append(p[1])
                self.p0.append(p[0])
                self.p0.append(1.0)
                self.p_orbitals.append(p_orb(p[2]))

        # For each d-orbital, get initial parameter values and create a d-orb object
        if 'd' in peaks:
            for d in peaks.get('d'):
                self.p0.append(d[1])
                self.p0.append(d[0])
                self.p0.append(1.0)
                self.d_orbitals.append(d_orb(d[2]))

        if self.bg_type == 'linear':
            self.p0.append(0.0)
            self.p0.append(0.0)

        if self.bg_type == 'active_shirley':
            self.p0.append(1.0)


    def analyticPeaks(self, x, *args):
        """
        The function to be fitted. Is a sum of Gaussians and the background according to how many orbitals were given
        and what type of orbitals they are.
        :param x: x-values
        :param args: Parameters to be determinied by curve-fitting (given in the list p0)
        :return: A sum of Gaussians and the background
        """

        solution = 0
        argindex = 0

        if 's' in self.peaks:
            for s in range(0, len(self.peaks.get('s'))):
                solution += self.s_orbitals[s].s_orbital(x, args[argindex], args[argindex + 1], args[argindex + 2])
                argindex += 3
        if 'p' in self.peaks:
            for p in range(0, len(self.peaks.get('p'))):
                solution += self.p_orbitals[p].split_p_orbital(x, args[argindex], args[argindex + 1], args[argindex + 2])
                argindex += 3
        if 'd' in self.peaks:
            for d in range(0, len(self.peaks.get('d'))):
                solution += self.d_orbitals[d].split_d_orbital(x, args[argindex], args[argindex + 1], args[argindex + 2])
                argindex += 3

        if self.bg_type == "linear":
            solution += self.linear_background(x, args[argindex], args[argindex + 1])
            argindex += 2

        if self.bg_type == "active_shirley":
            solution += self.shirley_active(args[argindex])
            argindex += 1

        if self.bg_type == "passive_shirley":
            solution += self.shirley_passive()

        return solution

    def peakfit(self):
        """
        Calls scipy.optimize.curve_fit
        :return: popt - Parameters for the function to curve-fit, pcov - Covariance
        """

        popt, pcov = curve_fit(self.analyticPeaks, self.x, self.y, p0=self.p0)
        return popt, pcov

    def plotPeakfit(self):
        """
        Plots the curve-fitted peak with the individual orbital Gaussians and the background
        :return:
        """

        popt, pcov = self.peakfit()

        plt.plot(self.x, self.y, 'x', label='Data')
        plt.plot(self.x, self.analyticPeaks(self.x, *popt), label='Curve fit')

        argindex = 0
        if 's' in self.peaks:
            for s in range(0, len(self.peaks.get('s'))):
                plt.plot(self.x, gauss(self.x, popt[argindex], popt[argindex + 1], popt[argindex + 2]), color='C3',
                         linestyle='--')
                argindex += 3
        if 'p' in self.peaks:
            for p in range(0, len(self.peaks.get('p'))):
                plt.plot(self.x, gauss(self.x, popt[argindex], popt[argindex + 1], popt[argindex + 2]), color='C4',
                         linestyle='--')
                plt.plot(self.x, gauss(self.x, 1 / 2 * popt[argindex], popt[argindex + 1] + self.p_orbitals[p].dsep(), popt[argindex + 2]),
                         color='C4', linestyle='--')
                argindex += 3
        if 'd' in self.peaks:
            for d in range(0, len(self.peaks.get('d'))):
                plt.plot(self.x, gauss(self.x, popt[argindex], popt[argindex + 1], popt[argindex + 2]), color='C5',
                         linestyle='--')
                plt.plot(self.x, gauss(self.x, 2 / 3 * popt[argindex], popt[argindex + 1] + self.d_orbitals[d].dsep(), popt[argindex + 2]),
                         color='C5', linestyle='--')
                argindex += 3

        plt.plot(self.x, self.background, 'k:')
        plt.grid()
        plt.legend()
        plt.show()

    # Doesn't work
    def shirley_active(self, k):
        B = np.zeros(self.x.shape[0])

        for n in range(0, 10):

            for i in range(0, B.shape[0]):
                B[i] = k * trapz(self.y[i:] - self.y[-1] - B[i:], x=self.x[i:])

        self.background = B + self.y[-1]
        return B + self.y[-1]

    def shirley_passive(self):
        """
        Gives a passive Shirley background
        :return: Passive Shirley background
        """
        B = np.zeros(self.x.shape[0])

        for n in range(0, 10):

            k = (self.y[0] - self.y[-1])/trapz(self.y - self.y[-1], x=self.x)
            for i in range(0, B.shape[0]):
                B[i] = k * trapz(self.y[i:] - self.y[-1] - B[i:], x=self.x[i:])

        self.background = B + self.y[-1]
        return B + self.y[-1]

    def linear_background(self, x, a, b):
        self.background = a * x + b
        return a * x + b

