'''Methods for calculating and graphing results of M^2 measurements.'''

import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat

def remove_wings(marginals, kernel_length, comparison_length, error, symmetric=False):
    '''Remove wings from 1D data by a rolling average comparison.'''
    y0_index = np.argmax(marginals)

    # Right limit
    max_index = len(marginals) - y0_index - kernel_length
    right_averages = np.zeros(max_index)

    right_limit = len(marginals)
    for i in range(0, max_index):
        rolling_average = np.mean(marginals[y0_index + i:y0_index + kernel_length + i])
        right_averages[i] = rolling_average

        if i <= comparison_length:
            continue

        deviation = np.mean(np.abs(rolling_average - right_averages[i - comparison_length:i]))
        if deviation < error:
            right_limit = y0_index - comparison_length + i
            break

    # Left limit
    min_index = -(y0_index - kernel_length)
    left_averages = np.zeros(np.abs(min_index))

    left_limit = 0
    for i in range(0, min_index, -1):
        rolling_average = np.mean(marginals[y0_index - kernel_length + i:y0_index + i])
        left_averages[np.abs(i)] = rolling_average

        if np.abs(i) <= comparison_length:
            continue

        deviation = np.mean(np.abs(rolling_average - left_averages[np.abs(i) - comparison_length:np.abs(i)]))
        if deviation < error:
            left_limit = y0_index + comparison_length + i
            break

    right_noise = np.mean(marginals[right_limit:])
    left_noise = np.mean(marginals[:left_limit])
    baseline_noise = np.min((left_noise, right_noise))
    if np.isnan(baseline_noise):
        baseline_noise = 0
    if symmetric:
        max_travel = np.max((right_limit - y0_index, left_limit - y0_index))
        left_limit = y0_index - max_travel
        right_limit = y0_index + max_travel
    marginals = marginals[left_limit:right_limit] - baseline_noise


    return marginals, left_limit, right_limit

def get_Wy2(irradiance, y, dx, dy,
            kernel_length=10, comparison_length=3, error=1E-1, symmetric=False):
    '''Obtain waist squared in y-direction (Wy2) through D4σ method.'''
    # Normalize irradiance
    irradiance = (irradiance - np.min(irradiance)) / (np.max(irradiance) - np.min(irradiance))
    marginals = trapezoid(irradiance, dx=dx)

    # Denoising
    original_marginals = marginals.copy()
    marginals, left_limit, right_limit = remove_wings(marginals, kernel_length, comparison_length, error, symmetric)
    y = y[left_limit:right_limit]

    # Total power
    P = trapezoid(marginals, dx=dy)
    # Normalize marginals
    marginals = marginals / P

    # Center of mass
    y0 = trapezoid(y * marginals, dx=dy)
    # Standard deviation squared
    sigma2_y = trapezoid((y - y0)**2 * marginals, dx=dy)
    # Waist squared
    Wy2 = 4 * sigma2_y

    return Wy2, y0, y, marginals, original_marginals

def get_Wx2(irradiance, x, dx, dy,
            kernel_length=10, comparison_length=3, error=1E-1, symmetric=False):
    '''Obtain waist squared in x-direction (Wy2) through D4σ method.'''
    return get_Wy2(irradiance.T, x, dy, dx, kernel_length, comparison_length, error, symmetric)

def W2_fit(z, W2s, wl):
    '''Fit beam parameters to waists as a function of propagation direction.'''
    W2_fit = lambda z, z0, W0, M2: W0**2 + M2**2 * (wl / (np.pi * W0))**2 * (z - z0)**2
    # Initial guesses:
    z00 = z[len(z)//2]
    W00 = np.min(W2s)
    M20 = 1
    popt, pcov = curve_fit(W2_fit, z, W2s, p0=(z00, W00, M20))
    z0, W0, M2 = popt
    sigma_z0, sigma_W0, sigma_M2 = np.sqrt(np.diagonal(pcov))

    z0 = ufloat(z0, sigma_z0)
    W0 = ufloat(W0, sigma_W0)
    M2 = ufloat(M2, sigma_M2)

    W2_func = lambda z: W2_fit(z, *popt)

    return z0, W0, M2, W2_func


'''Plotting methods.'''
def save_1D_irradiance(y0, y, marginals_y, x0, x, marginals_x, title, save=False):
    fig, ax = plt.subplots()

    ax.scatter(y, marginals_y, s=5, color='#b92f4f', label="Marginals(y)")
    ax.axvline(y0, color='#b92f4f', ls="--", alpha=0.3)
    ax.scatter(x, marginals_x, s=5, color='#0097d8', label="Marginals(x)")
    plt.axvline(x0, color='#0097d8', ls="--", alpha=0.3)
    ax.set(title=title, xlabel="Position $y$ ó $x$ [mm]", ylabel="Normalized 1D Irradiance [u.a]")
    ax.legend()

    if save:
        plt.savefig(f"Figures/Analysis/{title}.png", dpi=200)

def plot_W2_fit(z, W2s, W2_func, z0, W0, M2, direction: str, color, fig, ax):
    z0, W0, M2 = ufloat_to_str(z0), ufloat_to_str(W0), ufloat_to_str(M2)
    scatter_label = f"$W^2_{direction}(z)$"
    plot_label = f"$z_{{0, {direction}}} = {z0}$ mm\n" \
                 f"$W_{{0, {direction}}} = {W0}$ mm\n" \
                 f"$M^2_{direction} = {M2}$"
    ax.scatter(z, W2s, label=scatter_label, s=5, color=color)
    ax.plot(z, W2_func(z), label=plot_label, color=color)
    ax.grid(True, ls="--", color="#ccc")
    ax.legend(fontsize=6)

def ufloat_to_str(X):
    Xn, Xs = X.n, X.s
    Xn, Xs = E_to_10pow(f"{float(Xn):.2E}"), E_to_10pow(f"{float(Xs):.2E}")
    string = f"({Xn}) \pm ({Xs})"
    string = string.replace("E-0", "E-").replace("+0", "")
    string = string.replace("E0)", ")")
    string = string.replace("E", "\\mathrm{E}")
    return string

def E_to_10pow(string):
    mantissa, exponent = string.split("E")
    exponent = int(exponent)
    if exponent == 0:
        return mantissa
    if exponent == 1:
        return mantissa + "\\times 10"
    return mantissa + f"\\times 10^{{{exponent}}}"


if __name__ == '__main__':
    pass
