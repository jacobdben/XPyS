from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def gauss(x, C, mu, sigma):
    return C/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def s_orbital(x, C, mu, sig):
    return gauss(x, C, mu, sig)

def split_p_orbital(x, C, mu1, mu2, sig):
    return gauss(x, C, mu1, sig) + gauss(x, 1/2*C, mu2, sig)

def split_d_orbital(x, C, mu1, mu2, sig):
    return gauss(x, C, mu1, sig) + gauss(x, 2/3*C, mu2, sig)

def linear_background(x, a, b):
    return a*x + b


class XPS:
    def __init__(self, peaks, background, x_data, y_data):
        """
        :param peaks: dict containing type of orbital with list of approx position and height of each peak.
        Has the form {"s" : [ [peak1pos, peak1height], [peak2pos, peak2height], ...], "p" : ..., "d": ...}
        :param background: String containing type of background, e.g "linear" or "shirley"
        """
        self.peaks = peaks
        self.background = background
        self.x = x_data
        self.y = y_data

        self.p0 = []

        if 's' in peaks:
            for s in peaks.get('s'):
                self.p0.append(s[1])
                self.p0.append(s[0])
                self.p0.append(1.0)

        if 'p' in peaks:
            for p in peaks.get('p'):
                self.p0.append(p[1])
                self.p0.append(p[0])
                self.p0.append(p[0] + 1.0)
                self.p0.append(1.0)

        if 'd' in peaks:
            for d in peaks.get('d'):
                self.p0.append(d[1])
                self.p0.append(d[0])
                self.p0.append(d[0] + 1.0)
                self.p0.append(1.0)

        if background == 'linear':
            self.p0.append(1.0)
            self.p0.append(1.0)


    def analyticPeaks(self, x, *args):
        solution = 0
        argindex = 0

        if 's' in self.peaks:
            for s in range(0, len(self.peaks.get('s'))):
                solution += s_orbital(x, args[argindex], args[argindex + 1], args[argindex + 2])
                argindex += 3
        if 'p' in self.peaks:
            for p in range(0, len(self.peaks.get('p'))):
                solution += split_p_orbital(x, args[argindex], args[argindex + 1], args[argindex + 2], args[argindex + 3])
                argindex += 4
        if 'd' in self.peaks:
            for d in range(0, len(self.peaks.get('d'))):
                solution += split_d_orbital(x, args[argindex], args[argindex + 1], args[argindex + 2], args[argindex + 3])
                argindex += 4

        if self.background == "linear":
            solution += linear_background(x, args[argindex], args[argindex + 1])
            argindex += 1

        return solution

    def peakfit(self):

        popt, pcov = curve_fit(self.analyticPeaks, self.x, self.y, p0=self.p0)



        return popt, pcov

    def plotPeakfit(self):

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
                plt.plot(self.x, gauss(self.x, popt[argindex], popt[argindex + 1], popt[argindex + 3]), color='C4',
                         linestyle='--')
                plt.plot(self.x, gauss(self.x, 1 / 2 * popt[argindex], popt[argindex + 2], popt[argindex + 3]),
                         color='C4', linestyle='--')
                argindex += 4
        if 'd' in self.peaks:
            for d in range(0, len(self.peaks.get('d'))):
                plt.plot(self.x, gauss(self.x, popt[argindex], popt[argindex + 1], popt[argindex + 3]), color='C5',
                         linestyle='--')
                plt.plot(self.x, gauss(self.x, 2 / 3 * popt[argindex], popt[argindex + 2], popt[argindex + 3]),
                         color='C5', linestyle='--')
                argindex += 4

        plt.grid()
        plt.legend()
        plt.show()
