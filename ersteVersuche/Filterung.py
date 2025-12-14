import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
from DiscreteFourierTransform import myfft, gaussianLowPassFilter
import time
import math

def convolve2d(image:np.ndarray, kernel: np.ndarray,verbose: bool=True)->np.ndarray:
    kh, kw = kernel.shape
    pad_h = kh // 2 
    pad_w = kw // 2

    # Bild mit Nullen auffüllen (Padding)
    padded = np.pad(image.astype(np.float32), ((pad_h, pad_h), (pad_w, pad_w)), mode="constant",constant_values=0)

    output = np.zeros_like(image,dtype=np.float32)

    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            region = padded[y:y+kh, x:x+kw]
            output[y, x] = np.sum(region * kernel)


    disp = np.abs(output)
    disp= 255* (disp/disp.max())
    disp = disp.astype(np.uint8)

    if verbose:
        cv2.imshow("convolve2D filtered image",disp)
        cv2.waitKey(0)  
        cv2.destroyAllWindows()

    return output

import numpy as np

def convolve2d_3x3_vectorized(image: np.ndarray,
                              kernel: np.ndarray,
                              verbose:bool = True) -> np.ndarray:
    """
    Vektorisierte 2D-Faltung für beliebige 3x3-Kernel.

    Parameter
    ----------
    image : np.ndarray (H x W)
        Eingabebild (Graustufen), vorzugsweise float32 oder float64

    kernel : np.ndarray (3 x 3)
        Faltungskernel (z.B. Laplace, Sobel, Schärfen, etc.)

    padding : str
        Padding-Modus für np.pad (z.B. 'constant', 'edge', 'reflect')

    Returns
    -------
    output : np.ndarray (H x W)
        Gefaltetes Bild (float), ohne Anzeige-Normalisierung
    """

    # --- Sicherheitschecks ---
    if kernel.shape != (3, 3):
        raise ValueError("Kernel muss die Form (3, 3) haben")

    # Sicherstellen, dass mit float gerechnet wird
    image = image.astype(np.float32)

    # --- Padding ---
    # 1 Pixel Rand auf jeder Seite, damit alle Nachbarn existieren
    padded = np.pad(image, 1,mode="constant",constant_values=0)

    # --- Kernel-Koeffizienten ---
    # Explizit benannt für Lesbarkeit
    k00, k01, k02 = kernel[0, 0], kernel[0, 1], kernel[0, 2]
    k10, k11, k12 = kernel[1, 0], kernel[1, 1], kernel[1, 2]
    k20, k21, k22 = kernel[2, 0], kernel[2, 1], kernel[2, 2]

    # --- Vektorisierte Faltung ---
    # Jede Zeile greift auf einen verschobenen View des gepaddeten Bildes zu
    output = (
        k00 * padded[0:-2, 0:-2] +  # oben links
        k01 * padded[0:-2, 1:-1] +  # oben
        k02 * padded[0:-2, 2:  ] +  # oben rechts
        k10 * padded[1:-1, 0:-2] +  # links
        k11 * padded[1:-1, 1:-1] +  # Zentrum
        k12 * padded[1:-1, 2:  ] +  # rechts
        k20 * padded[2:  , 0:-2] +  # unten links
        k21 * padded[2:  , 1:-1] +  # unten
        k22 * padded[2:  , 2:  ]    # unten rechts
    )

    disp = np.abs(output)
    disp= 255* (disp/disp.max())
    disp = disp.astype(np.uint8)

    if verbose:
        cv2.imshow("convolve2d_3x3_vectorized image",disp)
        cv2.waitKey(0)  
        cv2.destroyAllWindows()

    return output


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    x_t = cv2.imread("lena.bmp",cv2.IMREAD_UNCHANGED)
    sigmaInFreqDomain = 50
    sigmaInSpacialDomain = 1/(2*np.pi*sigmaInFreqDomain)
    x_t_gausFiltered,x_t_gausFiltered_norm = gaussianLowPassFilter(x_t,sigma=sigmaInFreqDomain,verbose=True)
    x_t_gausFiltered2 = cv2.GaussianBlur(x_t,ksize=(0,0),sigmaX=sigmaInSpacialDomain,sigmaY=sigmaInSpacialDomain,borderType=cv2.BORDER_WRAP)
    # vergleiche x_t_gausFiltered und x_t_gausFiltered2
    difference = np.abs(x_t_gausFiltered - x_t_gausFiltered2)  
    print(f"Max difference between custom Gaussian filter and cv2.GaussianBlur: {np.max(difference)}")

    laplace_kernel = np.array([[0, 1, 0],
                               [1,-4, 1],
                               [0, 1, 0]],dtype=np.float32)

    Tstart1 = time.monotonic_ns()
    X_t_filtered = convolve2d(x_t_gausFiltered2,laplace_kernel,verbose=True)
    Tend1_start2 = time.monotonic_ns()
    X_t_filtered = convolve2d_3x3_vectorized(x_t_gausFiltered2,laplace_kernel,verbose=True)
    Tend2_start3 = time.monotonic_ns()
    X_t_filtered = cv2.filter2D(x_t_gausFiltered2, -1, laplace_kernel)
    Tend3 = time.monotonic_ns()
    print(f"convolve2d time: {(Tend1_start2 - Tstart1)/1000000:.4f} miliseconds")
    print(f"convolve2d_3x3_vectorized time: {(Tend2_start3 - Tend1_start2)/1000000:.4f} miliseconds")
    print(f"cv2.filter2D time: {(Tend3 - Tend2_start3)/1000000:.4f} miliseconds")

    cv2.imshow("cv2filter",X_t_filtered.astype(np.uint8))
    cv2.waitKey(0)  
    cv2.destroyAllWindows()



if __name__ =="__main__":
    main()