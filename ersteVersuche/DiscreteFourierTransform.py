import numpy as np, matplotlib.pyplot as plt
import cv2
import os



def dispImg(img:np.ndarray,title:str = "dispImg Output"):
    temp_img = img.copy()
    temp_img = 255*(temp_img/temp_img.max())
    cv2.imshow(title, temp_img.astype(dtype=np.uint8))
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    

def dispFreqDomainAmplitude(X_f:np.ndarray,title:str = "Magnitude F(u,v)"):
    X_mag = np.abs(X_f)
    X_mag = np.log1p(X_mag)
    X_mag = 255 * X_mag / np.max(X_mag)
    X_mag = X_mag.astype(np.uint8)
    cv2.imshow(title, X_mag)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def myfft(img:np.ndarray, verbose:bool=True):
    """
    Berechnet die 2D-Fouriertransformation eines Graustufenbildes und stellt
    Amplituden- und Phasenspektrum bereit.

    Die Funktion führt folgende Schritte aus:
    1. 2D-FFT des Eingabebildes
    2. Zentrierung der Nullfrequenz (fftshift)
    3. Berechnung und logarithmische Skalierung des Amplitudenspektrums
    4. Berechnung des Phasenspektrums

    Parameter
    ----------
    img : np.ndarray (H x W)
        Eingabebild im Ortsraum (Graustufen).
        Der Datentyp sollte numerisch sein (z. B. uint8 oder float).

    verbose : bool, optional
        Wenn True, wird das Amplitudenspektrum mit OpenCV angezeigt.
        Default: True.

    Returns
    -------
    X_f : np.ndarray (H x W), complex
        Komplexes Fourier-Spektrum des Bildes mit zentrierter Nullfrequenz.

    X_mag : np.ndarray (H x W), uint8
        Logarithmisch skaliertes Amplitudenspektrum zur Visualisierung,
        normalisiert auf den Wertebereich [0, 255].

    X_phase : np.ndarray (H x W), float
        Phasenspektrum in Grad (−180° bis 180°).

    Notes
    -----
    - Die logarithmische Skalierung (log1p) dient der besseren Darstellung
      des Amplitudenspektrums, da Fourier-Amplituden typischerweise einen
      sehr großen Dynamikbereich besitzen.
    - Das zurückgegebene Amplitudenspektrum (X_mag) ist ausschließlich
      für Visualisierungszwecke geeignet und nicht für numerische
      Weiterverarbeitung.
    - Die FFT geht implizit von periodischen Randbedingungen aus.

    """
    X_f = np.fft.fft2(img)
    X_f = np.fft.fftshift(X_f) # verschiebt die 0-Frequenz des Spektrums in die Bildmitte

    X_mag = np.log1p(np.abs(X_f)) # np.abs zur Darstellung der Amplitude, np.log zur Skalierung, log1p zur vermeidung von log(0)
    X_mag = 255 * X_mag / np.max(X_mag)
    X_mag = X_mag.astype(np.uint8)

    X_phase = np.angle(X_f, deg=True)

    if verbose:
        cv2.imshow("Amplitude F(u,v)", X_mag)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    return X_f, X_mag, X_phase


def myinvfft(X: np.ndarray, verbose: bool = True) -> np.ndarray:
    """
    Berechnet die inverse 2D-Fouriertransformation eines zentrierten
    Spektrums und rekonstruiert das Bild im Ortsraum.

    Parameter
    ----------
    X : np.ndarray (H x W), complex
        Fourier-Spektrum mit zentrierter Nullfrequenz (fftshifted).

    verbose : bool
        Wenn True, wird das rekonstruierte Bild zu Anzeigezwecken dargestellt.

    Returns
    -------
    x : np.ndarray (H x W), float
        Rekonstruiertes Bild im Ortsraum (float, nicht normalisiert).
    """

    X_backshifted = np.fft.ifftshift(X)
    x = np.fft.ifft2(X_backshifted)

    # numerisch kleiner Imaginärteil → verwerfen
    x = np.real(x)

    if verbose:
        disp = x - x.min()
        disp = disp / disp.max() * 255
        disp = disp.astype(np.uint8)

        cv2.imshow("x_reconstructed", disp)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    return x

def lowPassFilter(img: np.ndarray,filterShape: tuple,verbose=True) -> np.ndarray:

    X_f,X_mag, X_phase = myfft(img,verbose=False)
    
    rows,cols = img.shape
    lowPass_mask = np.zeros([rows,cols],dtype=np.uint8)
    lowPass_mask[(rows//2-filterShape[0]//2):(rows//2+filterShape[0]//2),(cols//2-filterShape[1]//2):(cols//2+filterShape[1]//2)] = 1
     
    # multiply X_mag with lowpass filter:
    X_f_filtered = X_f*lowPass_mask
    X_t_filtered = np.fft.ifft2(np.fft.ifftshift(X_f_filtered))
    X_t_filtered = np.abs(X_t_filtered)
    X_t_filtered = (255 * X_t_filtered / np.max(X_t_filtered)).astype(np.uint8)
    if verbose:
        cv2.imshow("lowpass mask",lowPass_mask*255)
        cv2.waitKey(0) 
        cv2.imshow("X_t_filtered",np.abs(X_t_filtered))
        cv2.waitKey(0)  
        cv2.destroyAllWindows()
    
    return X_t_filtered

def gaussianLowPassFilter(img: np.ndarray, sigma: float,verbose=True) -> np.ndarray:
    X_f,X_mag, X_phase = myfft(img,verbose=False)

    rows, cols = img.shape
    crow, ccol = rows // 2, cols // 2

    # Koordinatenraster (Distanzmatrix D(u,v))
    u = np.arange(rows)
    v = np.arange(cols)
    V, U = np.meshgrid(v, u)
    
    # H(u,v)=e^(D(u,v)^2 / 2*σ^2)​ ; D(u,v) = D2 = quadratische Entferung zum Mittelpunkt(crow|ccol)
    D2 = (U - crow)**2 + (V - ccol)**2

    H = np.exp(-(D2) / (2 * (sigma**2)))

    # Anwenden auf komplexes Spektrum
    X_f_filtered = X_f * H

    # inverse FFT
    X_back = np.fft.ifft2(np.fft.ifftshift(X_f_filtered))
    X_back = np.abs(X_back)

    # normalisieren auf 0–255
    X_back_norm = 255 * X_back / np.max(X_back)
    # --- verbose 3D Visualization ---
    if verbose:
        # Gaussian Filter im Ortsraum ist ebenfalls ein 2D Gaussian
        h_spatial = np.fft.ifft2(np.fft.ifftshift(H))
        h_spatial = np.abs(h_spatial)
        h_spatial = np.fft.fftshift(h_spatial)
        h_spatial = h_spatial/h_spatial.max()

        fig = plt.figure(figsize=(15, 4))

        # --- Plot 1: 3D Frequency Domain Filter ---
        ax1 = fig.add_subplot(131, projection="3d")
        ax1.plot_surface(U, V, H, cmap='viridis')
        ax1.set_title("Gaussian Lowpass (Frequency Domain)")
        ax1.set_xlabel("u")
        ax1.set_ylabel("v")
        ax1.set_zlabel("H(u,v)")

        # --- Plot 2: 3D Spatial Domain Filter ---
        ax2 = fig.add_subplot(132, projection="3d")
        ax2.plot_surface(U, V, h_spatial, cmap='magma')
        ax2.set_title("Gaussian Kernel (Spatial Domain)")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_zlabel("h(x,y)")

        # --- Plot 3: Filtered Image ---
        ax3 = fig.add_subplot(133)
        ax3.imshow(X_back_norm.astype(np.uint8), cmap='gray')
        ax3.set_title("Filtered Image")
        ax3.axis("off")

        plt.tight_layout()
        plt.show()

    return X_back,X_back_norm

def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    x_t = cv2.imread("lena.bmp",cv2.IMREAD_UNCHANGED)
    # X_f,X_mag, X_phase = myfft(x_t)

    # X_t_filtered = lowPassFilter(x_t,(100,100))
    X_t_filtered,X_t_filtered_normiert = gaussianLowPassFilter(x_t,sigma=50)

if __name__ =="__main__":
    main()