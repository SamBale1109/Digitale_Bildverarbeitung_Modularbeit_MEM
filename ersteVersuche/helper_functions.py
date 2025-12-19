
import numpy as np
import matplotlib.pyplot as plt


def plotPDF(z,pdf,mu,sigma):
    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(z, pdf, label=r"$p(z)$ Gaussian PDF")

    # Markante Punkte
    plt.axvline(mu, linestyle="--", label=r"Mean $\mu$")
    plt.axvline(mu - sigma, linestyle=":", label=r"$\mu - \sigma$")
    plt.axvline(mu + sigma, linestyle=":", label=r"$\mu + \sigma$")
    plt.axvline(mu - 2*sigma, linestyle=":", alpha=0.7, label=r"$\mu \pm 2\sigma$")
    plt.axvline(mu + 2*sigma, linestyle=":", alpha=0.7)

    # Achsenbeschriftung
    plt.xlabel("Intensity value $z$")
    plt.ylabel("Probability density $p(z)$")

    # Titel
    plt.title("Gaussian Probability Density Function")

    # Gitter & Legende
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()

    plt.tight_layout()
    plt.show()

def dispImg(img:np.ndarray,title:str = "dispImg:"):
        real_img = np.real(img)
        real_img = 255*(real_img/real_img.max())
        plt.imshow(real_img, cmap="gray")
        plt.title(title)
        plt.axis("off")
        plt.show()

def dispFreqDomainAmplitude(X_f:np.ndarray,title:str = "Magnitude F(u,v)"):
    X = X_f.copy()
    X = np.fft.ifftshift(X)
    X_mag = np.abs(X)
    X_mag = np.log1p(X_mag)
    X_mag = 255 * X_mag / np.max(X_mag)
    X_mag = X_mag.astype(np.uint8)
    plt.imshow(X_mag, cmap="gray")
    plt.title(title)
    plt.axis("off")
    plt.show()