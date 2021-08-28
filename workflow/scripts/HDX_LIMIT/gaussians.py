from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import numpy as np


def Gauss_1(xdata, H, A, mean, std):
    return H + A * np.exp(-(xdata - mean) ** 2 / (2 * std ** 2))


def Gauss_2(xdata, H1, A1, mean1, std1, A2, mean2, std2, ):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2))


def Gauss_3(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2))


def Gauss_4(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2))


def Gauss_5(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4, A5, mean5, std5):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2)) + A5 * np.exp(-(xdata - mean5) ** 2 / (2 * std5 ** 2))


def Gauss_6(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4, A5, mean5, std5,
            A6, mean6, std6):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2)) + A5 * np.exp(-(xdata - mean5) ** 2 / (2 * std5 ** 2)) + A6 * np.exp(
        -(xdata - mean6) ** 2 / (2 * std6 ** 2))


def Gauss_7(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4, A5, mean5, std5,
            A6, mean6, std6, A7, mean7, std7):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2)) + A5 * np.exp(-(xdata - mean5) ** 2 / (2 * std5 ** 2)) + A6 * np.exp(
        -(xdata - mean6) ** 2 / (2 * std6 ** 2)) + A7 * np.exp(-(xdata - mean7) ** 2 / (2 * std7 ** 2))


def Gauss_8(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4, A5, mean5, std5,
            A6, mean6, std6, A7, mean7, std7, A8, mean8, std8):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2)) + A5 * np.exp(-(xdata - mean5) ** 2 / (2 * std5 ** 2)) + A6 * np.exp(
        -(xdata - mean6) ** 2 / (2 * std6 ** 2)) + A7 * np.exp(-(xdata - mean7) ** 2 / (2 * std7 ** 2)) + A8 * np.exp(
        -(xdata - mean8) ** 2 / (2 * std8 ** 2))


def Gauss_9(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4, A5, mean5, std5,
            A6, mean6, std6, A7, mean7, std7, A8, mean8, std8, A9, mean9, std9):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2)) + A5 * np.exp(-(xdata - mean5) ** 2 / (2 * std5 ** 2)) + A6 * np.exp(
        -(xdata - mean6) ** 2 / (2 * std6 ** 2)) + A7 * np.exp(-(xdata - mean7) ** 2 / (2 * std7 ** 2)) + A8 * np.exp(
        -(xdata - mean8) ** 2 / (2 * std8 ** 2)) + A9 * np.exp(-(xdata - mean9) ** 2 / (2 * std9 ** 2))


def Gauss_10(xdata, H1, A1, mean1, std1, A2, mean2, std2, A3, mean3, std3, A4, mean4, std4, A5, mean5, std5,
             A6, mean6, std6, A7, mean7, std7, A8, mean8, std8, A9, mean9, std9, A10, mean10, std10):
    return H1 + A1 * np.exp(-(xdata - mean1) ** 2 / (2 * std1 ** 2)) + A2 * np.exp(
        -(xdata - mean2) ** 2 / (2 * std2 ** 2)) + A3 * np.exp(-(xdata - mean3) ** 2 / (2 * std3 ** 2)) + A4 * np.exp(
        -(xdata - mean4) ** 2 / (2 * std4 ** 2)) + A5 * np.exp(-(xdata - mean5) ** 2 / (2 * std5 ** 2)) + A6 * np.exp(
        -(xdata - mean6) ** 2 / (2 * std6 ** 2)) + A7 * np.exp(-(xdata - mean7) ** 2 / (2 * std7 ** 2)) + A8 * np.exp(
        -(xdata - mean8) ** 2 / (2 * std8 ** 2)) + A9 * np.exp(-(xdata - mean9) ** 2 / (2 * std9 ** 2)) + A10 * np.exp(
        -(xdata - mean10) ** 2 / (2 * std10 ** 2))


def get_popt(f, prominence, width):
    std_guess = 3
    std_min = 1
    std_max = 15
    tol_min = 0.99
    tol_max = 1.01

    means = find_peaks(f, prominence=prominence, width=width)[0]

    if len(means) == 0:
        return []
    elif len(means) == 1:
        popt, pcov = curve_fit(Gauss_1, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], 1],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], 20]))
    elif len(means) == 2:
        popt, pcov = curve_fit(Gauss_2, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], 1, tol_min * f[means[1]],
                                        tol_min * means[1], 1],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1],
                                        20]))
    elif len(means) == 3:
        popt, pcov = curve_fit(Gauss_3, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], 1, tol_min * f[means[1]],
                                        tol_min * means[1], 1,
                                        tol_min * f[means[2]], tol_min * means[2], 1],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], 20]))

    elif len(means) == 4:
        popt, pcov = curve_fit(Gauss_4, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3,
                                                                  f[means[3]], means[3], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], 1, tol_min * f[means[1]],
                                        tol_min * means[1], 1,
                                        tol_min * f[means[2]], tol_min * means[2], 1, tol_min * f[means[3]],
                                        tol_min * means[3], 1],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], 20]))
    elif len(means) == 5:
        popt, pcov = curve_fit(Gauss_5, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3,
                                                                  f[means[3]], means[3], 3, f[means[4]], means[4], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], 1, tol_min * f[means[1]],
                                        tol_min * means[1], 1,
                                        tol_min * f[means[2]], tol_min * means[2], 1, tol_min * f[means[3]],
                                        tol_min * means[3], 1,
                                        tol_min * f[means[4]], tol_min * means[4], 1],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], std_max,
                                        tol_max * f[means[4]], tol_max * means[4], std_max]))
    elif len(means) == 6:
        popt, pcov = curve_fit(Gauss_6, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3,
                                                                  f[means[3]], means[3], 3, f[means[4]], means[4], 3,
                                                                  f[means[5]], means[5], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], std_min, tol_min * f[means[1]],
                                        tol_min * means[1], std_min,
                                        tol_min * f[means[2]], tol_min * means[2], std_min, tol_min * f[means[3]],
                                        tol_min * means[3], std_min,
                                        tol_min * f[means[4]], tol_min * means[4], std_min, tol_min * f[means[5]],
                                        tol_min * means[5], std_min],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], std_max,
                                        tol_max * f[means[4]], tol_max * means[4], std_max, tol_max * f[means[5]],
                                        tol_max * means[5], std_max]))
    elif len(means) == 7:
        popt, pcov = curve_fit(Gauss_7, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3,
                                                                  f[means[3]], means[3], 3, f[means[4]], means[4], 3,
                                                                  f[means[5]], means[5], 3, f[means[6]], means[6], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], std_min, tol_min * f[means[1]],
                                        tol_min * means[1], std_min,
                                        tol_min * f[means[2]], tol_min * means[2], std_min, tol_min * f[means[3]],
                                        tol_min * means[3], std_min,
                                        tol_min * f[means[4]], tol_min * means[4], std_min, tol_min * f[means[5]],
                                        tol_min * means[5], std_min,
                                        tol_min * f[means[6]], tol_min * means[6], std_min],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], std_max,
                                        tol_max * f[means[4]], tol_max * means[4], std_max, tol_max * f[means[5]],
                                        tol_max * means[5], std_max,
                                        tol_max * f[means[6]], tol_max * means[6], std_max]))
    elif len(means) == 8:
        popt, pcov = curve_fit(Gauss_8, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3,
                                                                  f[means[3]], means[3], 3, f[means[4]], means[4], 3,
                                                                  f[means[5]], means[5], 3, f[means[6]], means[6], 3,
                                                                  f[means[7]], means[7], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], std_min, tol_min * f[means[1]],
                                        tol_min * means[1], std_min,
                                        tol_min * f[means[2]], tol_min * means[2], std_min, tol_min * f[means[3]],
                                        tol_min * means[3], std_min,
                                        tol_min * f[means[4]], tol_min * means[4], std_min, tol_min * f[means[5]],
                                        tol_min * means[5], std_min,
                                        tol_min * f[means[6]], tol_min * means[6], std_min, tol_min * f[means[7]],
                                        tol_min * means[7], std_min],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], std_max,
                                        tol_max * f[means[4]], tol_max * means[4], std_max, tol_max * f[means[5]],
                                        tol_max * means[5], std_max,
                                        tol_max * f[means[6]], tol_max * means[6], std_max, tol_max * f[means[7]],
                                        tol_max * means[7], std_max]))
    elif len(means) == 9:
        popt, pcov = curve_fit(Gauss_9, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                  f[means[1]], means[1], 3,
                                                                  f[means[2]], means[2], 3,
                                                                  f[means[3]], means[3], 3, f[means[4]], means[4], 3,
                                                                  f[means[5]], means[5], 3, f[means[6]], means[6], 3,
                                                                  f[means[7]], means[7], 3, f[means[8]], means[8], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], std_min, tol_min * f[means[1]],
                                        tol_min * means[1], std_min,
                                        tol_min * f[means[2]], tol_min * means[2], std_min, tol_min * f[means[3]],
                                        tol_min * means[3], std_min,
                                        tol_min * f[means[4]], tol_min * means[4], std_min, tol_min * f[means[5]],
                                        tol_min * means[5], std_min,
                                        tol_min * f[means[6]], tol_min * means[6], std_min, tol_min * f[means[7]],
                                        tol_min * means[7], std_min,
                                        tol_min * f[means[8]], tol_min * means[8], std_min],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], std_max,
                                        tol_max * f[means[4]], tol_max * means[4], std_max, tol_max * f[means[5]],
                                        tol_max * means[5], std_max,
                                        tol_max * f[means[6]], tol_max * means[6], std_max, tol_max * f[means[7]],
                                        tol_max * means[7], std_max,
                                        tol_max * f[means[8]], tol_max * means[8], std_max]))
    else:
        popt, pcov = curve_fit(Gauss_10, np.arange(len(f)), f, p0=[0, f[means[0]], means[0], 3,
                                                                   f[means[1]], means[1], 3,
                                                                   f[means[2]], means[2], 3,
                                                                   f[means[3]], means[3], 3, f[means[4]], means[4], 3,
                                                                   f[means[5]], means[5], 3, f[means[6]], means[6], 3,
                                                                   f[means[7]], means[7], 3, f[means[8]], means[8], 3,
                                                                   f[means[9]], means[9], 3],
                               bounds=([0, tol_min * f[means[0]], tol_min * means[0], std_min, tol_min * f[means[1]],
                                        tol_min * means[1], std_min,
                                        tol_min * f[means[2]], tol_min * means[2], std_min, tol_min * f[means[3]],
                                        tol_min * means[3], std_min,
                                        tol_min * f[means[4]], tol_min * means[4], std_min, tol_min * f[means[5]],
                                        tol_min * means[5], std_min,
                                        tol_min * f[means[6]], tol_min * means[6], std_min, tol_min * f[means[7]],
                                        tol_min * means[7], std_min,
                                        tol_min * f[means[8]], tol_min * means[8], std_min, tol_min * f[means[9]],
                                        tol_min * means[9], std_min],
                                       [0.0005, tol_max * f[means[0]], tol_max * means[0], std_max,
                                        tol_max * f[means[1]],
                                        tol_max * means[1], std_max,
                                        tol_max * f[means[2]], tol_max * means[2], std_max, tol_max * f[means[3]],
                                        tol_max * means[3], std_max,
                                        tol_max * f[means[4]], tol_max * means[4], std_max, tol_max * f[means[5]],
                                        tol_max * means[5], std_max,
                                        tol_max * f[means[6]], tol_max * means[6], std_max, tol_max * f[means[7]],
                                        tol_max * means[7], std_max,
                                        tol_max * f[means[8]], tol_max * means[8], std_max, tol_max * f[means[9]],
                                        tol_max * means[9], std_max]))
    return popt


def fit_gaussians(factor, mz_data, prominence=0.001, width=0.5):
    factor_max = max(factor)
    mz_data_max = max(mz_data)

    popts = get_popt(factor / factor_max, prominence=prominence, width=width)

    integrated_mz_gaussian = []
    mz_data_gaussian = []

    if popts == []:
        return []
    else:
        params = []
        for idx, _ in enumerate(range(len(popts) // 3)):
            params.append(np.concatenate([np.reshape(popts[0], 1), np.array(popts[3 * idx + 1:3 * idx + 4])]))
        params.sort(key=lambda x: x[1], reverse=True)
        for idx, param in enumerate(params):
            # Generated filtered IC from factor
            mask = Gauss_1(np.arange(len(factor)), *param)
            factor_filtered = (mask * (factor / factor_max) / max(mask * (factor / factor_max))) * factor[
                np.argmax(mask)]
            if sum(factor_filtered) > 5:
                integrated_mz_gaussian.append(factor_filtered)

                # Generate filtered mz data
                mask = np.repeat(mask, 7)
                mz_data_filtered = (mask * (mz_data / mz_data_max) / max(mask * (mz_data / mz_data_max))) * mz_data[
                    np.argmax(mask)]
                mz_data_gaussian.append(mz_data_filtered)


        return integrated_mz_gaussian, mz_data_gaussian
