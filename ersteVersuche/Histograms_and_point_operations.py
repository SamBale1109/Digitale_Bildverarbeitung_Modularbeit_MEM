import cv2
import numpy as np
import os
import sympy as sp
import matplotlib.pyplot as plt

def create_histogram(img,verbose=True,title="Histogramm der Graustufen"):
    hist, bins = np.histogram(img.flatten(),bins=256, range=(0,256), density=True)
    if verbose:
        plt.figure(figsize=(8,4))
        plt.plot(hist, color='black')
        plt.title(title)
        plt.xlabel("Grauwert")
        plt.ylabel("Häufigkeit")
        plt.grid(True)
        plt.show()
    return hist

def stretch_histogram(unstretched_img,hist,verbose=True):
    mask = hist > 0 # create [True,False] array
    if np.any(mask): # if any True exists
        L = np.argmax(mask) # idx of first True occurance
        H = len(mask) - 1 - np.argmax((mask)[::-1]) # idx of last True occurance
        print(f"L: {L} | H: {H}")
    else: 
        L = 0
        H = 256
    stretched_img = (unstretched_img.astype(np.float32)-L)/(H-L)
    stretched_img = np.clip(stretched_img,0,1) # clip values to 0-1
    stretched_img = np.uint8(stretched_img*255) # scale back to 0-255
    stretched_hist = create_histogram(stretched_img,verbose=verbose,title= "stretched Histogram")
    return stretched_img, stretched_hist


# Image Histogram as a probability mass function:
# Cumulative Distribution Function (CDF) sagt, wie viele Pixel haben Grauwert ≤ i
# Für eine diskrete Zufallsvariable X:
# CDF(x)=P(X≤x) → also die Wahrscheinlichkeit, dass der Grauwert(X) kleiner oder gleich x ist.
# Histogram soll waagerechte gerade annähern (alle Graustufenwerte gleich wahrscheinlich)
def hist_equalize(img, hist, verbose=True):
    # sicherstellen, dass Graustufenbild
    img = img.astype(np.uint8)

    # Kumulative Verteilung (CDF)
    cdf = hist.cumsum()

    # Nullwerte der CDF entfernen (damit untere Grenzen nicht stören)
    cdf_nonzero = cdf[cdf > 0]

    # Normalisierte CDF (0–255)
    cdf_norm = (cdf - cdf_nonzero[0]) / (cdf_nonzero[-1] - cdf_nonzero[0]) * 255
    cdf_norm = cdf_norm.astype(np.uint8)

    # Lookup: Für jeden Pixel neuen Wert zuweisen
    equalized_img = cdf_norm[img]
    equalized_hist = create_histogram(equalized_img,verbose=True,title= "equalized Histogram")

    return equalized_img, equalized_hist


os.chdir(os.path.dirname(os.path.abspath(__file__)))

unstretched_img = cv2.imread("unequalized.jpg",cv2.IMREAD_UNCHANGED)
hist = create_histogram(unstretched_img)

stretched_img, stretched_histogram = stretch_histogram(unstretched_img,hist)
equalized_img, equalized_hist = hist_equalize(unstretched_img,hist)

cv2.imshow("stretched_img",stretched_img)
cv2.waitKey(0)  # wartet auf Tastendruck
cv2.imshow("equalized_img",equalized_img)
cv2.waitKey(0)
cv2.destroyAllWindows()