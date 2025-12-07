import cv2
import numpy as np
import os
import sympy as sp
import matplotlib.pyplot as plt

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def create_histogram(img,verbose=True):
    hist, bins = np.histogram(img.flatten(),bins=256, range=(0,256), density=True)
    if verbose:
        plt.figure(figsize=(8,4))
        plt.plot(hist, color='black')
        plt.title("Histogramm der Graustufen")
        plt.xlabel("Grauwert")
        plt.ylabel("Häufigkeit")
        plt.grid(True)
        plt.show()
    return hist

def quantize_bad(img,anzGrauwerte):
    nrows,ncols = np.shape(img)
    f,q,a,b = sp.symbols("f q a b")
    fun = (f-q)**2
    E = sp.integrate(fun,(f,a,b))
    dE = sp.diff(E,q)
    q_opt = sp.solve(sp.simplify(dE), q)[0]
    print("q_opt =", q_opt)
    step = 256//anzGrauwerte
    Intervalle = [(A*step,(A+1)*step) for A in range(anzGrauwerte)]
    print(f"Intervalle: {Intervalle}")
    q_values = [int(q_opt.subs({a: I[0], b: I[1]})) for I in Intervalle]
    # mappe alle werte in Interval[0] auf q_values[0]
    qant_img = np.zeros([nrows,ncols],dtype=np.uint8)
    for x in range(nrows):
        for y in range(ncols):
            for interval_idx,I in enumerate(Intervalle):
                if I[0] <= img[x, y] < I[1]:
                    qant_img[x, y] = q_values[interval_idx]
                    break
    
    return qant_img

def quantize(img,anzGrauwerte):
    f,q,a,b = sp.symbols("f q a b")
    fun = (f-q)**2
    E = sp.integrate(fun,(f,a,b))
    dE = sp.diff(E,q)
    q_opt = sp.solve(sp.simplify(dE), q)[0]
    print("q_opt =", q_opt)

    step = 256//anzGrauwerte
    Intervalle = [(A*step,(A+1)*step) for A in range(anzGrauwerte)]
    print(f"Intervalle: {Intervalle}")

    q_values = np.array([int(q_opt.subs({a: I[0], b: I[1]})) for I in Intervalle],dtype=np.uint8)

    bounds = np.array([A for (A,B) in Intervalle] + [256])
    # np.digitize ordnet jedes Element eines Arrays einem Intervallindex zu. -> array aus indizes
    # durch fancy indexing wird für jeden idx in idx_array der entsrechende q_value aus q_values genommen und array der selben struktur wie idx-array erzeugt
    idx_array = np.digitize(img,bins=bounds[:-1])
    quant_img = q_values[idx_array-1]

    return quant_img



lena_bmp = cv2.imread("lena.bmp",cv2.IMREAD_UNCHANGED)
hist_lena = create_histogram(lena_bmp)

quant_lena = quantize(lena_bmp,16)
hist_lena_quant = create_histogram(quant_lena)

cv2.imshow("quant_lena",quant_lena)
cv2.waitKey(0)  # wartet auf Tastendruck
cv2.destroyAllWindows()
