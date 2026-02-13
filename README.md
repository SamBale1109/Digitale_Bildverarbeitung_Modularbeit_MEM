## Table of Contents

- [Wiener-Filter in der digitalen Bildverarbeitung](#wiener-filter-in-der-digitalen-bildverarbeitung)
  - [Statistische Annahmen](#statistische-annahmen)
  - [Bildmodell (Degradation + Rauschen)](#bildmodell-degradation--rauschen)
  - [Optimierungsprinzip (MMSE und Orthogonalitätsprinzip)](#optimierungsprinzip-mmse-und-orthogonalitätsprinzip)
  - [Herleitung](#herleitung)
  - [Verhalten des Wiener-Filters](#verhalten-des-wiener-filters)
  - [Ermittlung des Rausch-Signal-Verhältnis](#ermittlung-des-rausch-signal-verhältnis)
    - [Schätzung der Rauschstatistik](#schätzung-der-rauschstatistik)
  - [Degradation des Bildes](#degradation-des-bildes)
    - [Lineare Bewegungsunschärfe (engl. Motion Blur)](#lineare-bewegungsunschärfe-engl-motion-blur)
    - [Gaußsches Rauschen](#gaußsches-rauschen)
- [Bildrekonstruktion](#bildrekonstruktion)
  - [Bekannte Bilddaten (theoretischer Idealfall)](#bekannte-bilddaten-theoretischer-idealfall)
  - [Unbekannte Bilddaten (Realfall)](#unbekannte-bilddaten-realfall)
  - [Artefakte](#artefakte)
  - [Fazit](#fazit)
  - [Anhang: MATLAB-Implementierung des Wiener-Filters](#anhang-matlab-implementierung-des-wiener-filters)

\renewcommand{\thepage}{\arabic{page}}
\setcounter{page}{1}

# Wiener-Filter in der digitalen Bildverarbeitung
Der Wiener-Filter, auch als Minimum Mean Square Error (MMSEMMSE) -Filter bezeichnet, ist ein stochastisches Verfahren zur Bildrestaurierung. Ziel ist es, aus einem verschlechterten Bild eine statistisch optimale Schätzung $\hat{f}$\nomenclature{$\hat{f}(x,y)$}{Restauriertes Bild} des ursprünglichen Bildes $f$$f(x, y)$ zu bestimmen. 
Im Gegensatz zu reinen Inversionsverfahren berücksichtigt der Wiener-Filter die statistischen Eigenschaften von Bild und Rauschen explizit. Eine frequenzabhängige Gewichtung auf Basis der Leistungsdichtespektren dämpft stark verrauschte Frequenzanteile und bevorzugt informationshaltige Komponenten, wodurch im statistischen Mittel ein optimales Restaurierungsergebnis erzielt wird.
## Statistische Annahmen

Der Wiener-Filter basiert auf dem Ansatz, Bilder und Rauschen nicht als deterministische Größen, sondern als Zufallsvariablen zu behandeln. Dafür werden folgende zentrale Annahmen getroffen:
\begin{itemize}
    \item Unkorreliertheit von Bild und Rauschen: Es wird vorausgesetzt, dass das Rauschen und das Bildsignal statistisch unabhängig voneinander sind.
    \item Zero Mean: Die Herleitung geht davon aus, dass entweder das Bildsignal oder das Rauschen (oder beide) einen Mittelwert von Null aufweisen
    \item Linearität: Die Schätzung des optimalen Bildes wird als eine lineare Funktion der Intensitätswerte des gestörten Bildes betrachtet. Auch der zugrunde liegende Verschlechterungsprozess (die Degradation H$H(u,v)$) wird als ein lineares System modelliert.\cite{Gonzalez.2018}
\end{itemize}

## Bildmodell (Degradation + Rauschen)
Zunächst wird das lineare ortsinvariante Modell für ein verschmiertes und verrauschtes Bild aufgestellt:

$$
g(x,y) = h(x,y) * f(x,y) + n(x,y),$g(x,y)$$n(x,y)$ \tag{1}
$$

wobei $f(x,y)$ das Originalbild, $h(x,y)$ $h(x,y)$die Degradationsfunktion, $n(x,y)$ additives Rauschen und $g(x,y)$ das degradierte Bild ist. Das Symbol $*$ bezeichnet die zweidimensionale Faltung.
Durch Anwendung der diskreten zweidimensionalen Fourier-Transformation wird in den Frequenzraum transformiert:

$$
G(u,v) = H(u,v)F(u,v) + N(u,v).$G(u,v)$$F(u,v)$$N(u,v)$ \tag{2}
$$

Die Restauration des Bildes erfolgt durch einen linearen Filter $W(u,v)$ $W(u,v)$im Frequenzraum. Das restaurierte Bild $\widehat{F}(u,v)$ ergibt sich dabei zu:

$$
\widehat{F}(u,v) = W(u,v)G(u,v),\nomenclature{$\widehat{F}(u,v)$}{Restauriertes Bild} \tag{3}
$$



## Optimierungsprinzip (MMSE und Orthogonalitätsprinzip)

Die Bestimmung des Wiener-Filters basiert auf dem Orthogonalitätsprinzip. Es besagt, dass eine Schätzung $\hat{f}$ des ursprünglichen Bildes $f$ genau dann optimal im Sinne des minimalen mittleren quadratischen Fehlers ist, wenn der verbleibende Schätzfehler orthogonal, also statistisch unkorreliert, zum beobachteten degradierten Bild $g$ ist.
Für Zufallsfelder in der Bildverarbeitung bedeutet Orthogonalität, dass der Erwartungswert des Produkts aus Schätzfehler und Beobachtung verschwindet:

$$
E\{(\hat{f}-f)g\} = 0 \nomenclature{$E\{\}$}{Erwartungswertoperator} \tag{4}
$$

Das degradierte Bild $g$ spannt den Raum der verfügbaren Bildinformationen auf, während die Schätzung $\hat{f}$ als Funktion dieser Beobachtungen innerhalb dieses Raumes liegt. Die Orthogonalitätsbedingung fordert, dass der minimale Schätzfehler senkrecht zu diesem Informationsraum steht und somit keine verwertbaren Bildinformationen mehr enthält. Der verbleibende Fehler ist dann reiner Zufall, der in keinerlei Beziehung mehr zu den verfügbaren Daten steht und somit mathematisch nicht weiter reduziert werden kann.\cite{Kailath.1981}

## Herleitung
In dieser Arbeit wird der Wiener-Filter nicht direkt über die Orthogonalitätsbedingung,
sondern über den äquivalenten Ansatz der expliziten Minimierung des mittleren
quadratischen Fehlers hergeleitet. Beide Vorgehensweisen führen auf dieselbe
Optimallösung. Der quadratische Fehler $e^2$ $e^2$wird definiert als:


$$
e^2 = (f(x,y)-\widehat{f}(x,y))^2 \tag{5}
$$

Durch Einsetzen von $\widehat{f}(x,y)=w(x,y)*g(x,y)$ und Transformation in den Frequenzraum ergibt sich der quadratische Fehler für eine einzelne Frequenz (u,v) (vereinfacht ohne Indizes) zu:

$$
|E_f|^2 = |F - WG|^2
= (F-WG)(F-WG)^* \tag{6}
$$

wobei $(\cdot)^*$ die komplexe Konjugation bezeichnet.\\
Setzt man das Modell~ ein, erhält man:

$$
E_f = F - W(HF + N) = F(1 - WH) - WN \tag{7}
$$


$$
|E_f|^2 = [F(1 - WH) - WN] \cdot [F(1 - WH) - WN]^* \tag{8}
$$

Das ergibt ausmultipliziert:
%\begin{enumerate}
%\item $|F(1 - WH)|^2$ (Quadrat des ersten Teils)
%\item $-F(1 - WH)(WN)^*$ (Kreuzterm 1)
%\item $-(F(1 - WH))^*(WN)$ (Kreuzterm 2)
%\item $+|WN|^2$ (Quadrat des zweiten Teils)
%\end{enumerate}

$$
\begin{aligned}
|E_f|^2
&= |F(1 - WH)|^2
- F(1 - WH)(WN)^* \quad \text{(Kreuzterm 1)}
- (F(1 - WH))^*(WN) \quad \text{(Kreuzterm 2)}
+ |WN|^2
\end{aligned} \tag{9}
$$

Zur Berechnung des Mean Square Error (MSE)MSE wird der Erwartungswert-Operator angewandt:

$$
MSE = E\{|E_f|^2\} \tag{10}
$$


$$
MSE = E\{|F(1 - WH)|^2-F(1 - WH)(WN)^*-(F(1 - WH))^*(WN)+|WN|^2\} \notag \tag{11}
$$

Hier kommen die statistischen Annahmen zum Tragen:
Die angenommene Unkorreliertheit des Rauschens und des Bildsignals hat zur Folge, dass der Erwartungswert des Produkts von Bild und Rauschen Null ergibt und die Kreuzterme 1 und 2 entfallen:

$$
E\{F \cdot N^*\} = 0 \quad \text{und} \quad E\{F^* \cdot N\} = 0 \tag{12}
$$


$$
MSE=E\{|F(1 - WH)|^2\}+E\{|WN|^2\} \tag{13}
$$

Da $W$ und $H$ feste Filterfunktionen (keine Zufallsvariablen) sind, können konstante Terme vor den Erwartungswert gezogen werden:

$$
MSE = |1 - WH|^2 \cdot E\{|F|^2\} + |W|^2 \cdot E\{|N|^2\} \tag{14}
$$

Diese Erwartungswerte der quadrierten Beträge entsprechen den Leistungsdichtespektren:
\begin{itemize}
\item $S_f(u,v) = \mathbb{E}\{|F(u,v)|^2\}$: Leistungsspektrum des Originalbildes.
\item $S_\eta(u,v) = \mathbb{E}\{|N(u,v)|^2\}$: Leistungsspektrum des Rauschens.
\end{itemize}
Durch Ausmultiplizieren mit ($|1 - WH|^2 = (1 - WH)(1 - W^* H^*)$) folgt:

$$
MSE = (1 - WH - W^* H^* + |W|^2 |H|^2)\, S_f + |W|^2 S_\eta \tag{15}
$$

Da der Filterkoeffizient $W(u,v)$ im Frequenzbereich komplex ist, handelt es sich bei der MSE-Funktion um eine reellwertige Funktion einer komplexen Variablen. Eine direkte Ableitung nach $W$ ist daher nicht zulässig. Stattdessen wird der Wirtinger-Kalkül verwendet, bei dem $W$ und $W^*$ als formal unabhängige Variablen betrachtet werden. 
Für reellwertige Kostenfunktionen ist das Minimum erreicht, wenn die Ableitung nach dem komplex konjugierten Parameter verschwindet\cite{Brandwood.1983}. Aus diesem Grund wird das Minimum des Fehlers gefunden indem der Ausdruck nach $W^*$ abgeleitet wird und das Ergebnis gleich Null gesetzt wird:

$$
\frac{\partial MSE}{\partial W^*}=- H^* S_f+ W |H|^2 S_f+ W S_\eta= 0 \tag{16}
$$

Die resultierende Gleichung kann nach $W$ aufgelöst werden und ergibt den Wiener-Filter im Frequenzbereich\cite{Wiener.1949}:

$$
W ( |H|^2 S_f + S_\eta ) = H^* S_f \notag \tag{17}
$$


$$
W(u,v)=\frac{H^*(u,v) S_f(u,v)}{|H(u,v)|^2 S_f(u,v) + S_\eta(u,v)} \tag{18}
$$

Da die Leistungsspektren von Bild und Rauschen meist unbekannt sind wird die Gleichung durch Kürzen mit $S_f(u,v)$ in eine vereinfachte Form mit dem Rausch-Signal-Verhältnis $\mathrm{NSR}(u,v)$(engl. Noise-to-Signal Ratio)NSR überführt:

$$
NSR(u,v) = \frac{Leistungsdichtespektrum\ des\ Rauschens}{Leistungsdichtespektrum\ des\ Originalbildes}=\frac{S_\eta(u,v)}{S_f(u,v)} \tag{19}
$$


$$
W(u,v) = \left[\frac{H^*(u,v)}{|H(u,v)|^2 + S_\eta(u,v)/S_f(u,v)}\right] \tag{20}
$$


$$
W(u,v) = \left[\frac{H^*(u,v)}{|H(u,v)|^2 + NSR(u,v)}\right] \tag{21}
$$


$$
W(u,v)=\frac{1}{H(u,v)}\left[ \frac{1}{1+\frac{NSR(u,v)}{|H(u,v)|^2}}\right] \tag{22}
$$



## Verhalten des Wiener-Filters
Aus der umgestellten Wiener-Filter-Gleichung~ geht hervor, dass der Wiener-Filter im Wesentlichen einem Inversfilter mit zusätzlichem, frequenzabhängigen Gewichtungsterm zur Steuerung der Rauschunterdrückung entspricht. Für ein hohes NSR nimmt der Gewichtungsterm einen kleinen Wert an, wodurch die entsprechende Frequenzkomponente stark abgeschwächt wird. Ist das NSR hingegen klein, d.\,h.\ der Rauschanteil ist gering, geht der Gewichtungsterm gegen Eins. In diesem Fall entspricht der Wiener-Filter näherungsweise dem Inversfilter.
Der im Gewichtungsterm enthaltene Faktor $|H(u,v)|^2$ stellt sicher, dass Frequenzanteile, die durch den Degradationsprozess stark gedämpft oder vollständig unterdrückt wurden, nicht verstärkt werden. Dadurch verhindert der Wiener-Filter eine unkontrollierte Verstärkung von Rauschen in Frequenzbereichen, in denen keine verlässliche Bildinformation vorhanden ist. Abbildung~ zeigt die Übertragungsfunktionen des Wiener-Filters für unterschiedliche konstante NSR-Werte. Sie verdeutlicht das verhalten eines Inversfilters  und dessen Ähnlichkeit mit einem Wiener Filter kleiner NSRs hin zu einer glättenden Filtercharakteristik bei größeren NSRs.
\begin{figure}[ht]
\centering
\includesvg[width=0.9\textwidth]{Bilder/Vergleich_Wiener_Invers_H}
\caption{Vergleich der Frequenzgänge von Inversfilter, Wiener-Filter und Degradationsfilter}

\end{figure}
\FloatBarrier

\clearpage
## Ermittlung des Rausch-Signal-Verhältnis
Da die Leistungsdichtespektren $S_f(u,v)$ $S_f(u,v)$und $S_\eta(u,v)$$S_\eta(u,v)$  für ein unbekanntes Originalbild in der Praxis nicht vorliegen, wird das eigentlich frequenzabhängige Rausch-Signal-Verhältnis häufig durch den skalaren Regularisierungsparameter K approximiert. Damit werden bei der Rekonstruktion alle Frequenzen, unabhängig des tatsächlichen Anteils des Rauschens, gleich behandelt. K kann entweder experimentell durch visuelle Beurteilung des rekonstruierten Bildes oder wie in dieser Arbeit durch eine Schätzung aus dem Quotienten der Varianzen von Rauschen und Signal ermittelt werden.
Dabei wird die Varianz als Maß der mittleren Leistung herangezogen, da sie im Gegensatz zum Mittelwert ausschließlich die Intensitätsschwankungen um den Gleichanteil und damit die strukturell relevanten AC-Anteile (Wechselanteile) wie Kanten und Texturen des Bildes beschreibt.

### Schätzung der Rauschstatistik
Zur Schätzung der Varianz des Rauschens $\sigma_n^2$ wird im degradierten Bild eine möglichst strukturlose, monochrome Bildregion ausgewählt, deren Pixel im Originalbild einen konstanten Grauwert aufweisen würden. Die vorhandenen Intensitätsschwankungen in diesem Ausschnitt können daher dem Rauschen zugeschrieben werden. Aus dem normierten Histogramm des Bildausschnitts wird so die Art und Varianz des Rauschens ermittelt.
% \begin{wrapfigure}{r}{0.4\textwidth}
%     \centering
%     \includegraphics[width=0.3\textwidth]{Bilder/BildMarkierterAusschnitt.png} 
%     \caption{Bildausschnitt aus $G(u,v)$}
%     
%     % \includegraphics[width=0.4\textwidth]{Bilder/HistogramRauschen2.png} 
%     % \caption{Normiertes Histogramm des Bildausschnitts. Die Intensitätsverteilung zeigt eine annähernd gaußförmige Verteilung, was auf additives gaußsches Rauschen hinweist}
%     % 
% \end{wrapfigure}
Für eine 8 Bit Darstellung mit $L=256$ möglichen Graustufen und den Intensitätswerten $z_i$ $z_i$
liefert das normierte Histogramm des Bildausschnitts die Wahrscheinlichkeitsverteilung $p_n(z_i)$,$p_n(z_i)$ aus der sich der Mittelwert $\hat{z}$ und die Varianz $\sigma_n^2$ des Rauschens berechnen lassen:

$$
\hat{z} = \sum_{i=0}^{L-1} z_i \, p_n(z_i). \tag{23}
$$


$$
\sigma_n^2 = \sum_{i=0}^{L-1} (z_i - \hat{z})^2 \, p_n(z_i). \tag{24}
$$


$$
\sigma_f^2 \approx \sigma_g^2 - \sigma_n^2. \tag{25}
$$


$$
\mathrm{K} \approx  \frac{\sigma_n^2}{\sigma_f^2}. \tag{26}
$$

wobei $\sigma_g^2$ $\sigma^2$die Varianz des degradierten Bildes und $\sigma_f^2$ die Varianz des Originalbildes bezeichnet\cite{Gonzalez.2018}.
\FloatBarrier

\clearpage
## Degradation des Bildes
Die mathematische Darstellung von Bildstörungen ist Voraussetzung für die Bildrestaurierung mithilfe des Wiener Filters. Physikalische Einflüsse wie Bewegungen zwischen Kamera und Bild während der Belichtung, aber auch Sensorrauschen müssen in ein Modell überführt werden. 
### Lineare Bewegungsunschärfe (engl. Motion Blur)
 Lineare Bewegungsunschärfe resultiert aus der relativen Bewegung zwischen Sensor und Objekt. Die Berechnung der Degradationsfunktion $H(u,v)$ basiert auf dem in
\cite{Gonzalez.2018} beschriebenen Modell. Unter der Annahme, dass sich der Verschluss (engl. Shutter) augenblicklich öffnet und schließt und der optische Abbildungsprozess ideal ist, können die durch die Bewegung verursachten Effekte isoliert werden. Das verschlechterte Bild $g(x,y)$ wird durch Integration der zeitabhängigen Verschiebungskomponenten $x_0(t)$ und $y_0(t)$ modelliert und ergibt sich mit der Belichtungszeit TT zu:
 

$$
g(x,y) = \int_{0}^{T} f\!\bigl[x - x_0(t),\, y - y_0(t)\bigr] \, dt \tag{27}
$$

% Um die Degradationsfunktion zu isolieren wird die kontinuierliche Fourier-Transformation auf das Integral angewendet:
% 
$$
% G(u,v) = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}
% \left[ \int_{0}^{T} f\!\bigl([x - x_0(t),\, y - y_0(t)\bigr]\, dt \right]
% e^{-j 2\pi (u x + v y)} \, dx \, dy
% \tag{28}
$$

% Durch Vertauschen der Intergrationsreihenfolge folgt:
% 
$$
% G(u,v) =  \int_{0}^{T}F(u,v)e^{-j 2\pi \left[ u x_0(t) + v y_0(t) \right]} \, dt
% \tag{29}
$$

% Zieht man nun $F(u,v)$ vor das zeitintegral kann die Degradationsfunktion dargestellt werden:
% 
$$
% H(u,v) = \frac{G(u,v)}{F(u,v)}=\int_{0}^{T} e^{-j 2\pi \left[ u x_0(t) + v y_0(t) \right]} \, dt
% 
% \tag{30}
$$

Angenommen die Lineare Bewegung in y-Richtung ist Null $y_0(t)=0$ und die Variable $x_0(t)$ ist bekannt, folgt für die Transferfunktion:

$$
H(u,v) = \int_{0}^{T}e^{-j 2\pi  u x_0(t) } \, dt \tag{31}
$$

Nehmen wir zusätzlich die Bewegung in x-Richtung als linear an, so bewegt sich das Bild mit der konstanten Geschwindigkeit $a/T$, sodass $x_0(t)=\frac{at}{T}$. Hier beschreibt a adie Verschiebungsdistanz während der Belichtung. Durch einsetzten in das Integral und Integration folgt:

$$
H(u,v) = \frac{T}{\pi u a } \,\sin \,(\pi u a ) \,e^{-j \pi u a } \tag{32}
$$

Um eine diskrete Übertragungsfunktion der Größe M Mx NN zu erhalten, muss die Gleichung über die gesamte Größe abgetastet werden. Durch Multiplikation mit dem Originalbild im Frequenzraum ergibt sich ein degradiertes Frequenzspektrum.
\begin{figure}[htbp]
    \centering
    \begin{subfigure}[t]{0.32\textwidth}
        \includegraphics[width=\linewidth]{Stefan Bilder/Originalbild.png}
        \caption{Originalbild}
    \end{subfigure}
    \hspace{0.02\textwidth}
    \begin{subfigure}[t]{0.32\textwidth}
        \includegraphics[width=\linewidth]{Stefan Bilder/Motion Blur.png}
        \caption{Bewegungsunschärfe}
    \end{subfigure}
    \caption{Originalbild (a) und Ergebnis der Bilddegradation in Form von Bewegungsunschärfe (b) mit $a=25$ und $T=1$.}
    
\end{figure}



### Gaußsches Rauschen

%Die Wahrscheinlichkeitsdichtefunktion PDF(PDF) des gaußschen Rauschens ist definiert als
% 
$$
% p(z) = \frac{1}{\sqrt{2\pi\sigma^2}} 
% \exp\left( -\frac{(z-\hat{z})^2}{2\sigma^2} \right),
% \tag{33}
$$

% wobei $z$ die Intensität eines Pixels, $\hat{z}$ \nomenclature{$\hat{z}$}{Mittelwert}den Mittelwert und $\sigma$ $\sigma$ die Standardabweichung beschreibt.  
Ein häufig verwendetes Rauschmodell in der digitalen Bildverarbeitung ist das gaußsche Rauschen. Es beschreibt zufällige Intensitätsschwankungen, die unter anderem durch Sensorrauschen oder elektronische Störeinflüsse entstehen können. In dieser Arbeit wird gaußsches Rauschen mit einem Mittelwert von $\mu = 0$ erzeugt, indem mehrere gleichverteilte Zufallsvariablen aufsummiert und anschließend normiert werden. Dieses Vorgehen ist zulässig, da die Summe unabhängiger Zufallsvariablen gegen eine Normalverteilung konvergiert. Die resultierende Rauschmatrix wird mit der gewünschten Standardabweichung skaliert und anschließend gemäß dem Model~ auf das degradierte Bild addiert.
%
$$
%g(x,y) = f(x,y)*h(x,y)+ n(x,y),
% \tag{34}
$$

%wobei $f(x,y)*h(x,y)$ das degradierte Bild in Form von Bewegungsunschärfe und $n(x,y)$ das gaußsche Rauschen beschreibt.  

Es wird additives gaußsches Rauschen mit einer Varianz von $\sigma^2 = 60$ verwendet. Die Standartabweichung des Rauschens beträgt demnach $7,746$ Grauwerte.
Dieses Rauschmodell beschreibt signalunabhängige Störungen, wie sie typischerweise bei Bildsensoren auftreten.  Für den Wiener Filter ist die Unabhänigkeit wie in Kapitel~ bereits beschrieben eine feste Annahme. Die Wahl der Varianz bestimmt die Stärke des Rauschens und erlaubt eine kontrollierte Simulation realistischer Bilddegradation. Das Rauschen wird additiv auf das degradierte Bild aufgebracht:
\begin{figure}[htbp]
    \centering
    \begin{subfigure}[t]{0.32\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/Originalbild.png}
        \caption{Originalbild}
    \end{subfigure}\hfill
    \begin{subfigure}[t]{0.32\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/Motion Blur.png}
        \caption{Bewegungsunschärfe}
    \end{subfigure}\hfill
    \begin{subfigure}[t]{0.32\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/noisyimage.png}
        \caption{Bewegungsunschärfe und Rauschen}
    \end{subfigure}
    \caption{Originalbild und degradierte Bilder nach Anwendung von Bewegungsunschärfe und anschließend additivem gaußschem Rauschen mit Varianz 60.}
    
\end{figure}
\FloatBarrier
# Bildrekonstruktion
Zur Rekonstruktion des Degradierten Bildes $g(x,y)$ wird der Wiener Filter aus Gleichung~ verwendet.
%
$$
%W(u,v)=\frac{H^{*}(u,v)}{|H(u,v)|^{2} + \frac{S_\eta(u,v)}{S_f(u,v)}}
% \tag{35}
$$

Durch Multiplikation mit dem Degradierten Bild und dem $NSR$ ergibt sich für das Rekonstruierte Bild im Frequenzraum $\hat{F}(u,v)$:

$$
\hat{F}(u,v) =
\frac{H^{*}(u,v)}{|H(u,v)|^{2} + NSR}
\, G(u,v) \tag{36}
$$

Die Degradationsfunktion wird in dieser Arbeit als bekannt angenommen. Diese Annahme ist in vielen praktischen Anwendungen zulässig, da bestimmte Bilddegradationen, wie beispielsweise Bewegungsunschärfe, aus den Aufnahmeeigenschaften der Kamera rekonstruiert oder mithilfe integrierter Inertialmesseinheiten bestimmt werden können.\cite{Gonzalez.2018} 
## Bekannte Bilddaten (theoretischer Idealfall)
Beginnend soll der Fall betrachtet werden, bei welchem sowohl das Originalbild als auch das Rauschen bekannt sind. Dieser Fall tritt ein, wenn das Bild wie hier erst künstlich verschlechtert wird. Das $NSR$ kann dabei exakt aus den Leistungsdichtespektren von Originalbild und Rauschen bestimmt werden. Anmerkend zu erwähnen ist, dass das $NSR$ in diesem Fall einer Matrix der Größe des Bildes entspricht. Das Rekonstruierte Bild wird im Frequenzraum berechnet und anschließend Rücktransformiert:
\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.35\textwidth]{Stefan Bilder/Rekonstruiertesbildknown.png}

        \smallskip
    {\small $\mathrm{NSR}\!:\,\frac{S_\eta(u,v)}{S_f(u,v)}$ \quad $\mathrm{MSE}\!:\,0{,}0051$}

    \caption{$\hat{f}_{\text{ideal}}(u,v)$: Ideale Wiener-Rekonstruktion}
    
\end{figure}
\FloatBarrier
Das Rekonstruierte Bild zeigt wesentlich schärfere Kanten. Jedoch ist noch ein deutlicher Unterschied zum Originalbild zu erkennen. Um diese Abweichungen mit einem einzigen Wert beschreiben zu können wird hier der MSE als Gütekriterium genutzt:

$$
\mathrm{MSE} = \frac{1}{MN} \sum_{x=1}^{M} \sum_{y=1}^{N}
\left( f(x,y) - \hat{f}(x,y) \right)^2 \tag{37}
$$

Für diese Berechnung wird auch das Originalbild benötigt. Der MSE wird auch im Folgenden genutzt, um die Qualität der Rekonstruktion zu bewerten und die Schätzung des Rausch Signal Verhältnisses zu validieren.
## Unbekannte Bilddaten (Realfall)
Im Folgenden wird der eigentliche Anwendungsfall des Wiener Filters beschrieben, bei dem sowohl das Originalbild, als auch das Rauschen unbekannt sind. Es gilt eine geeignete Schätzung für das NSR zu finden. Dazu wird wie in Kapitel~ beschrieben zunächst das normierte Histogram des Bildausschnitts ausgewertet(siehe Abbildung ~). Die Intensitätsverteilung zeigt eine annähernd gaußförmige Verteilung aus der sich eine Standardabweichung $\sigma_n$ von 8.5932 berechnet. Dieser Wert liegt nahe an der in ~ beschriebenen tatsächlichen Standardabweichung von 7,746. Analog wird $\sigma_g$ aus dem Histogramm von $g(x,y)$ geschätzt. Aus der Differenz der Varianzen wird $\sigma_f^2$ geschätzt und gemäß der Gleichung ~ der Regularisierungsparameter K gebildet. Mit dem geschätzten Rausch zu Signal Verhältnis von $\mathrm{K}\approx0,026$ wird das Bild rekonstruiert:
\begin{figure}[htbp]
    \centering
    \begin{subfigure}[t]{0.5\textwidth}
        \centering
        \includegraphics[width=1.0\linewidth]{Bilder/HistogramRauschen3.png}
        \captionsetup{width=0.8\textwidth}
        \caption{Normiertes Histogramm des Bildausschnitts}
        
    \end{subfigure}\hfill
    \begin{subfigure}[t]{0.5\textwidth}
        \centering

        \includegraphics[width=0.7\textwidth]{Stefan Bilder/Rekonstruiertesbild.png}
    %     \smallskip
    % {\small $\mathrm{K}\!:\,0{,}026 \quad \mathrm{MSE}\!:\,0{,}0094$}
        \captionsetup{width=0.8\textwidth}
        \caption{$\hat{f}_{\text{real}}(u,v)$: Reale Wiener-Rekonstruktion mit den Werten:\\$\mathrm{K}\!:\,0{,}026 \quad \mathrm{MSE}\!:\,0{,}0094$}
    \end{subfigure}\hfill
    

    
    %
\end{figure}
\FloatBarrier
Zur Validierung des berechneten K wird der Referenzwert durch mehrere benachbarte Werte ergänzt. 
Dazu werden vier zusätzliche Werte logarithmisch um den Referenzwert verteilt. Diese logarithmische Skalierung ermöglicht außerdem eine robuste Bewertung der Sensitivität der Wiener-Rekonstruktion gegenüber Abweichungen des NSR:
\begin{figure}[htbp]
    \centering
    \begin{subfigure}[t]{0.24\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/Rekon_NSR_0.00262_MSE_0.0511.png}
        \caption{$\mathrm{K}\!:\,0{,}0026$ \\ $\mathrm{MSE}\!:\,0{,}0511$}
    \end{subfigure}\hfill
    \begin{subfigure}[t]{0.24\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/Rekon_NSR_0.0122_MSE_0.0122.png}
        \caption{$\mathrm{K}\!:\,0{,}0122$ \\ $\mathrm{MSE}\!:\,0{,}0122$}
    \end{subfigure}\hfill
    \begin{subfigure}[t]{0.24\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/Rekon_NSR_0.0564_MSE_0.01.png}
        \caption{$\mathrm{K}\!:\,0{,}0564$ \\ $\mathrm{MSE}\!:\,0{,}01$}
    \end{subfigure}\hfill
    \begin{subfigure}[t]{0.24\textwidth}
        \centering
        \includegraphics[width=\linewidth]{Stefan Bilder/Rekon_NSR_0.262_MSE_0.023.png}
        \caption{$\mathrm{K}\!:\,0{,}262$ \\ $\mathrm{MSE}\!:\,0{,}032$}
    \end{subfigure}
    \caption{Logarithmische Variation von K: Der Faktor hat Einfluss auf die Charakteristik des Wiener Filters.}
    
\end{figure}
\FloatBarrier
Für kleine Werte von K werden Bilddetails stärker rekonstruiert. Dies führt zu einer höheren Schärfe, verstärkt jedoch gleichzeitig das vorhandene Rauschen. Für große Werte von K dominiert der Gewichtungsterm. Das resultierende Bild ist stärker geglättet, weist weniger Rauschen auf, verliert jedoch an Schärfe und Detailinformation.
Der Parameter K im Wiener-Filter ist immer ein Kompromiss zwischen Schärfe und Rauschunterdrückung. Die Wahl von K ist entscheidend und maßgeblich für die Stabilität des Filter wie bereits in Kapitel~ näher erläutert.
## Artefakte
Bei der Anwendung des Wiener-Filters zur Rekonstruktion von bewegungsunscharfen Bildern können typische Artefakte auftreten. 
Insbesondere sogenanntes Ghosting entsteht durch die Überbetonung bestimmter Frequenzanteile, die durch die Degradationsfunktion der Bewegungsunschärfezuvor stark abgeschwächt wurden. 
Da Kanten alle Frequenzanteile enthalten, werden sie bei der Inversion besonders empfindlich verstärkt und können zu Mehrfachkonturen oder Nachbildern führen. Obwohl der Wiener-Filter die direkte Division durch Null vermeidet, werden Frequenzen mit kleinen Beträgen von $|H(u,v)|$ weiterhin stark gewichtet, sofern der Regularisierungsterm nicht ausreichend groß gewählt ist.

Zusätzlich kann es zu sogenanntem Mottling kommen, bei dem homogene Bildbereiche fleckig oder inhomogen erscheinen, da hochfrequente Rauschanteile übermäßig verstärkt werden.
In den Rekonstruierten Bildern mit variierendem Regularisierungsfaktor K sind diese Artefakte gut zu erkennen. Das Ghosting ist in Form eines weiteren Umrisses um den Kameramann zu sehen. Motlling ist dort eher als normalverteiltes Rauschen und bei hohen Werten von K nur wenig zu sehen. Stärker ausgeprägt ist das Mottling bei der idealen Rekonstruktion mit bekanntem NSR. Der Grund dafür liegt in der frequenzabhängigen Gewichtung des Filters: In Frequenzbereichen, in denen das Bildspektrum kleine Werte annimmt, wird das NSR sehr groß, wodurch hochfrequente Rauschanteile lokal stark verstärkt werden. 
Diese frequenzselektive Überverstärkung führt insbesondere in zuvor homogenen Bildregionen zu inhomogenen, fleckenartigen Strukturen.
Im Gegensatz dazu wirkt ein konstanter NSR-Wert als eine gleichmäßige Dämpfung. Diese wirkt auch auf kritische Frequenzen und bewirkt somit eine visuell stabilere Rekonstruktionen.\cite{Gonzalez.2018}
## Fazit
Der Wiener-Filter stellt ein optimales Rekonstruktionsverfahren im Sinne des mittleren quadratischen Fehlers dar. Durch die explizite Berücksichtigung von Rauscheinflüssen ermöglicht er eine stabilere Rekonstruktion als reine Inversionsverfahren. Der Regularisierungsparameter steuert dabei den Kompromiss zwischen Schärferückgewinnung und Rauschunterdrückung. Kleine Werte von K führen zu schärferen Rekonstruktionen, verstärken jedoch das Rauschen, während große Werte eine stärkere Glättung bewirken. Die richtige Wahl des Parameters ist maßgeblich für die visuelle Qualität des rekonstruierten Bildes. Eine perfekte Rekonstruktion ist jedoch grundsätzlich nicht möglich, da Frequenzanteile mit $H(u,v)=0$ vollständig verloren gehen und nicht wiederhergestellt werden können.

## Anhang: MATLAB-Implementierung des Wiener-Filters
Die in den vorherigen Kapiteln beschriebenen Schritte zur Bilddegradation und Wiener-Rekonstruktion wurden vollständig in MATLAB umgesetzt. Der zugehörige Programmcode ist im Anhang in Form einer \texttt{.mlx}-Datei beigefügt. Standardfunktionen werden lediglich für grundlegende Operationen wie Fourier-Transformationen und Bilddarstellung verwendet, während die eigentliche Signalverarbeitung explizit selbst implementiert ist.
Die Bilddegradation kann wahlweise über eine gaußsche Punktspreizfunktion oder über lineare Bewegungsunschärfe erfolgen. In beiden Fällen wird die jeweilige Degradationsfunktion direkt aus dem zugrunde liegenden analytischen Modell berechnet und nicht über vorgefertigte Filterfunktionen erzeugt. Ebenso werden das additive gaußsche Rauschen, die NSR-Schätzung sowie die Wiener-Filterübertragungsfunktion vollständig manuell berechnet und zusammengesetzt. Die Rekonstruktion erfolgt ausschließlich über elementweise Operationen im Frequenzraum ohne Verwendung von MATLAB-internen Restaurations- oder Regularisierungsfunktionen.
Sowohl der Idealfall mit bekanntem, frequenzabhängigem NSR als auch der praxisnahe Realfall mit konstantem K sind implementiert.
Zur Bewertung der Ergebnisse werden rekonstruierte Bilder für unterschiedliche K-Werte erzeugt und über den mittlere quadratische Fehler mit dem Originalbild verglichen. Ergänzend werden Frequenzgangvergleiche zwischen Inversfilter, Wiener-Filter mit verschiedenen Werten von K und der Degradationsfunktion dargestellt.
