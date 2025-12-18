import cv2
import os
import numpy as np
from DiscreteFourierTransform import myfft,myinvfft, dispFreqDomainAmplitude

def add_periodicNoise(X_f:np.ndarray,noise_freqs:tuple,noiseAbs:float = 1e6)->np.ndarray:
    X_noise = X_f.copy()
    rows,cols = X_noise.shape
    for u0,v0 in noise_freqs:
        Z = X_noise[rows//2 +u0, cols//2 +v0]
        if np.abs(Z) < 1e-12:
            delta = noiseAbs + 0j
        else:
            delta = noiseAbs*Z/np.abs(Z)
        X_noise[rows//2 +u0, cols//2 +v0] = Z + delta
        # --- konjugierte Gegenfrequenz ---
        X_noise[rows//2 - u0, cols//2 - v0] = np.conj(X_noise[rows//2 + u0, cols//2 + v0])

    return X_noise


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    x = np.asarray(cv2.imread("lena.bmp",cv2.IMREAD_UNCHANGED),dtype=np.uint8)
    X_f,X_mag,X_phase = myfft(x)

    X_noise = add_periodicNoise(X_f,((10,10),(10,-10)),1e7)
    dispFreqDomainAmplitude(X_noise,title="noise")

    xback = myinvfft(X_noise)

    # Differenzanalyse
    diff = np.abs(x.astype(np.float32) - xback)

    print("Max difference:", np.max(diff))
    print("Mean difference:", np.mean(diff))



if __name__ =="__main__":
    main()