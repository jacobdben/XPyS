from xps import *

path = "XPS_data/"

# Change this according to what file you want to analyse
filename = "Rb3d.txt"

########### Read file (do not change) ############
file = open(path + filename, 'r')

x_data = np.array([])
y_data = np.array([])

if file.mode == 'r':
    lines = file.readlines()
    for line in lines:
        if line == '':
            continue
        x_val, y_val = line.split()
        x_data = np.append(x_data, [float(x_val)])
        y_data = np.append(y_data, [float(y_val)])

plt.plot(x_data, y_data, 'x', label='Data')
plt.show()
##################################################

# Here you set what orbitals to fit, and give guesses for position, height and doublet separation
# Warning: the initial guesses for position must be good (so try tweaking the position if the curve fit fails)
peaks = {'d': [[112, 1.0e+5, 1.39],[114, 1.0e+0, 1.39]]}

# Ignore this
window = XPS(peaks, 'passive_shirley', x_data, y_data)
window.peakfit()
window.plotPeakfit()