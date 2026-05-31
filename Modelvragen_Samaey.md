# Modelvragen Deel Samaey — Volledige Uitgewerkte Antwoorden

> Gebaseerd op de cursustekst Numerieke Benadering (KU Leuven, 2025–2026).
> Alle bewijzen zijn volledig uitgeschreven. Gebruik dit document als primaire studiereferentie.

---

## V1 — Kubische splinefunctie

### (a) Definitie

**Intuïtie:** Globale hoge-graadsveeltermen oscilleren sterk tussen knooppunten (Runge-fenomeen). Oplossing: verdeel het interval in stukjes en gebruik op elk stuk een veelterm van lage graad, aaneengekoppeld met de juiste vloeiendheid.

**Definitie 5.2.1 — Splinefunctie van graad $k$:**

Zij $a = t_0 < t_1 < \cdots < t_{n-1} < t_n = b$ een strikte knooppuntenrij.

Een **splinefunctie** $s(x)$ van graad $k > 0$ is een functie op $[a,b]$ met:
1. Op elk deelinterval $[t_{i-1}, t_i]$ is $s(x)$ een **veelterm van graad $\leq k$**.
2. De functie $s(x)$ en haar afgeleiden $s', s'', \ldots, s^{(k-1)}$ zijn **continu** op $[a,b]$.

Een **kubische spline** ($k = 3$) bestaat dus uit stukken van graad $\leq 3$, aaneengekoppeld met $C^2$-continuïteit.

> ⚠️ **Examentip:** Ken de definitie letterlijk. Twee eisen: (1) veelterm per stuk, (2) continuïteit t/m $k-1$-de afgeleide.

💡 **Waarom kubisch?** Het menselijk oog kan discontinuïteiten in de **kromming** (de tweede afgeleide) waarnemen. Om visueel vloeiende curves te bekomen, moeten $s$ en $s'$ en $s''$ continu zijn → kubisch ($k = 3$, $C^2$).

---

### (b) Dimensie van de splineruimte

> ⚠️ **Examentip:** Het bewijs is een eenvoudig telbewijs. Leer het van buiten.

**Stelling 5.2.1:** De vectorruimte van splinefuncties van graad $k$ met $n+1$ knooppunten heeft **dimensie $n + k$**.

**Bewijs (telbewijs):**

**Stap 1 — Vrijheidsgraden zonder continuïteit:**

$n$ deelintervallen, op elk stuk een veelterm van graad $\leq k$ → $k+1$ coëfficiënten per stuk:
$$\text{totaal} = n(k+1) \text{ parameters.}$$

**Stap 2 — Continuïteitsvoorwaarden:**

In elk van de $n-1$ **inwendige** knooppunten $t_1, \ldots, t_{n-1}$ worden $k$ voorwaarden opgelegd: continuïteit van $s, s', \ldots, s^{(k-1)}$.

$$\text{aantal beperkingen} = (n-1) \cdot k.$$

**Stap 3 — Vrijheidsgraden na de beperkingen:**

$$\dim = n(k+1) - (n-1)k = nk + n - nk + k = n + k. \quad \square$$

💡 **Intuïtie:** $n$ stukken + $k$ extra vrijheidsgraden voor het globale gedrag (bijv. $k$ randcondities vrij te kiezen voor een kubische spline: $n + 3$).

---

### (c) Basis van afgeknotte machtsfuncties en B-splines

#### Afgeknotte-machtsbasis

**Definitie:**
$$
(x - t_i)_+^k = \begin{cases} 0 & \text{als } x \leq t_i \\ (x - t_i)^k & \text{als } x > t_i \end{cases}
$$

De volledige basis bestaat uit:
$$\{1,\, x,\, x^2,\, \ldots,\, x^k,\, (x-t_1)_+^k,\, \ldots,\, (x-t_{n-1})_+^k\}.$$

Dat zijn $k+1 + (n-1) = n+k$ functies — klopt met de dimensie.

**Tekening (kubisch, $k=3$, $n=2$):**
```
     (x-t₁)³₊
      |
      |         ___/
      |        /
______|________/_____________ x
      t₁
```
Nul links van $t_1$, daarna $= (x-t_1)^3$. Elke functie "activeert" in één knooppunt.

**Nadeel:** Numeriek instabiel (slechte conditionering) — dit is de **"slechte basis"** waarnaar in examenvragen verwezen wordt. De volledige afgeknotte-machtsbasis bevat zowel de veeltermbasis $\{1, x, \ldots, x^k\}$ als de afgeknotte functies $\{(x-t_1)_+^k, \ldots, (x-t_{n-1})_+^k\}$ — het weglaten van de veeltermbasis is een veelgemaakte fout. Gebruik in de praktijk de B-splinebasis.

#### B-splinebasis

**Definitie (via gedeelde differentie):**

De **gewone B-spline** van graad $k$ met knooppunten $t_i < t_{i+1} < \cdots < t_{i+k+1}$:
$$M_{i,k+1}(x) = [t_i, t_{i+1}, \ldots, t_{i+k+1}]_t\, (t - x)_+^k.$$

De **genormaliseerde B-spline:**
$$\boxed{N_{i,k+1}(x) = (t_{i+k+1} - t_i)\, M_{i,k+1}(x).}$$

**Recursiebetrekking** (sleuteleigenschap voor berekening):
$$N_{i,k+1}(x) = \frac{x - t_i}{t_{i+k} - t_i}\, N_{i,k}(x) + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_{i+1}}\, N_{i+1,k}(x).$$

**Basisgeval** ($k=0$): $N_{i,1}(x) = \mathbf{1}_{[t_i, t_{i+1})}(x)$ (karakteristieke functie).

**Tekening (kubische B-splines, $k=3$):**
```
         N_{i,4}
           /\
          /  \
_________/    \_________
    t_i            t_{i+4}
```
Elke B-spline is nul buiten $[t_i, t_{i+k+1}]$ (lokale drager = $k+2$ knooppunten breed).

**Sleuteleigenschappen van genormaliseerde B-splines:**

| Eigenschap | Formule/uitspraak |
|---|---|
| Lokale drager | $N_{i,k+1}(x) = 0$ buiten $[t_i, t_{i+k+1}]$ |
| Positiviteit | $N_{i,k+1}(x) > 0$ voor $x \in (t_i, t_{i+k+1})$ |
| Partitie van eenheid | $\sum_{i=-k}^{n-1} N_{i,k+1}(x) = 1$ op $[t_0, t_n]$ |

---

### (d) B-splinebasis voor benadering van experimentele data

**Context:** Gegeven meetpunten $(x_r, f_r)$, $r = 1, \ldots, N$. We zoeken een spline $s(x)$ van graad $k$ die het **gewogen kleinste-kwadratenfout** minimaliseert:
$$\min_{c} \sum_{r=1}^N w_r\, \bigl(f_r - s(x_r)\bigr)^2.$$

**Voorstelling in B-splinebasis:**
$$s(x) = \sum_{i=-k}^{n-1} c_i\, N_{i,k+1}(x).$$

**Opstellen van het normaalstelsel:**

Substitueer in de kostfunctie:
$$\min_{c} \sum_{r=1}^N w_r\, \Bigl(f_r - \sum_i c_i N_{i,k+1}(x_r)\Bigr)^2.$$

Differenteer naar $c_j$ en stel gelijk aan nul:
$$\sum_r w_r N_{j,k+1}(x_r) \Bigl(\sum_i c_i N_{i,k+1}(x_r) - f_r\Bigr) = 0, \quad \forall j.$$

Dit is het normaalstelsel $Gc = b$ met:
$$G_{j,i} = \sum_{r=1}^N w_r N_{j,k+1}(x_r) N_{i,k+1}(x_r), \qquad b_j = \sum_{r=1}^N w_r N_{j,k+1}(x_r) f_r.$$

**Te kiezen parameters:**
- Graad $k$: doorgaans $k = 3$ (kubisch)
- Aantal en ligging knooppunten $t_0, \ldots, t_n$: te kiezen a.h.v. de data
- Coëfficiënten $c_i$: bepaald door het normaalstelsel

> ⚠️ **Examentip:** De knooppuntenkeuze is de moeilijkste stap. Te weinig → slechtere fit; te veel → overfitting. Noodzakelijke voorwaarde voor een unieke oplossing: **Schönberg-Whitney-voorwaarde** — er moet minstens één datapunt liggen in elk knooppunteninterval $[t_{i}, t_{i+1})$. Zonder dit kan het stelsel $Gc = b$ singulier worden.

---

### (e) Spaarsheidspatroon van het lineaire stelsel

**Sleuteleigenschap (lokale drager):** $N_{i,k+1}(x) = 0$ buiten $[t_i, t_{i+k+1}]$, een interval van $k+1$ deelintervallen. Daardoor is:
$$G_{j,i} = 0 \quad \text{als } |j - i| > k.$$

Het normaalstelsel is een **bandmatrix** met bandbreedte $2k+1$.

**Spaarsheidspatroon voor $k=3$, $n=5$ (7 B-splines):**
```
   [* * * *        ]
   [* * * * *      ]
   [* * * * * *    ]
   [* * * * * * *  ]
   [  * * * * * * *]
   [    * * * * * *]
   [      * * * * *]
```
(Enkel de 7 diagonalen $|i-j| \leq 3$ zijn niet-nul.)

> ⚠️ **Examentip:** Op het examen is enkel het patroon vereist, niet de waarden. Teken een $N \times N$ bandmatrix met breedte $2k+1$.

💡 **Numeriek voordeel:** De bandstructuur laat toe het normaalstelsel in $O(Nk^2)$ bewerkingen op te lossen met Cholesky-factorisatie voor bandmatrices, i.p.v. $O(N^3)$ voor een vol stelsel.

---

## V2 — Tweedimensionele veeltermbenaderering

**Probleemstelling:** Gegeven meetdata $f_{i,j} \approx f(x_i, y_j)$ op een rechthoekig rooster ($1 \leq i \leq M$, $1 \leq j \leq N$), met separabele gewichten $w_{i,j} = w_i^{(1)} w_j^{(2)}$. Zoek een benadering
$$z(x, y) = \sum_{k=0}^m \sum_{l=0}^n c_{k,l}\, x^k y^l.$$

### (a) Orthogonale basis voor het 2D-benaderingsprobleem

**Intuïtie:** Als we een monomiale basis gebruiken, krijgen we een slecht geconditioneerde Grammatrix (bijna-Hilbert). Oplossing: zoek **orthogonale veeltermen** in 2D.

**Aanpak — tensorproduct met separabele gewichten:**

> ⚠️ **Examentip:** De sleutel is de separabiliteit van het gewicht. Zonder die eigenschap werkt dit niet.

Stel $\phi_{k,l}(x,y) = p_k(x)\, q_l(y)$, waarbij:
- $\{p_k\}$ orthogonaal zijn t.o.v. discrete sommatie in $x$: $\sum_{i=1}^M w_i^{(1)} p_k(x_i) p_u(x_i) = 0$ voor $k \neq u$
- $\{q_l\}$ orthogonaal zijn t.o.v. discrete sommatie in $y$: $\sum_{j=1}^N w_j^{(2)} q_l(y_j) q_v(y_j) = 0$ voor $l \neq v$

**Bewijs van 2D-orthogonaliteit:**

$$\sum_{i=1}^M \sum_{j=1}^N w_{i,j}\, \phi_{k,l}(x_i, y_j)\, \phi_{u,v}(x_i, y_j)$$
$$= \sum_{i=1}^M \sum_{j=1}^N w_i^{(1)} w_j^{(2)}\, p_k(x_i) q_l(y_j)\, p_u(x_i) q_v(y_j)$$
$$= \underbrace{\left(\sum_{i=1}^M w_i^{(1)} p_k(x_i) p_u(x_i)\right)}_{=0 \text{ als } k \neq u} \cdot \underbrace{\left(\sum_{j=1}^N w_j^{(2)} q_l(y_j) q_v(y_j)\right)}_{\text{willekeurig}}$$

Dus: $\langle \phi_{k,l}, \phi_{u,v} \rangle = 0$ zolang $k \neq u$ **of** $l \neq v$.

De som is nul zodra minstens één van beide factoren nul is, d.w.z. zodra $(k,l) \neq (u,v)$. $\square$

**Opbouwen van $p_k$ en $q_l$:** via de Forsythe-methode (drietermsrecursiebetrekking):
$$p_k(x) = \lambda_k\left[(x - \alpha_k) p_{k-1}(x) - \beta_k p_{k-2}(x)\right],$$
met discrete recursiecoëfficiënten $\alpha_k$, $\beta_k$ uit de data $\{x_i, w_i^{(1)}\}$.

---

### (b) Coëfficiënten van de beste benadering

**Stelling:** Met orthogonale basis $\{\phi_{k,l} = p_k q_l\}$ herleiden de normaalstelsel-vergelijkingen zich tot de eenvoudige uitdrukking:
$$\boxed{a_{k,l} = \frac{\displaystyle\sum_{i=1}^M \sum_{j=1}^N w_i^{(1)} w_j^{(2)}\, f_{i,j}\, p_k(x_i)\, q_l(y_j)}{\|p_k\|^2 \|q_l\|^2}.}$$

**Bewijs:**

Het normaalstelsel voor de coëfficiënten $a_{k,l}$ (in de orthogonale basis) is:
$$\sum_{k', l'} a_{k',l'}\, \langle \phi_{k',l'}, \phi_{k,l} \rangle = \langle f, \phi_{k,l} \rangle.$$

Door de orthogonaliteit geldt $\langle \phi_{k',l'}, \phi_{k,l} \rangle = 0$ tenzij $(k',l') = (k,l)$. Het stelsel **diagonaliseert volledig**:
$$a_{k,l} \cdot \langle \phi_{k,l}, \phi_{k,l} \rangle = \langle f, \phi_{k,l} \rangle.$$

Uitrekenen:
$$\langle \phi_{k,l}, \phi_{k,l} \rangle = \sum_{i,j} w_{i,j} p_k(x_i)^2 q_l(y_j)^2 = \left(\sum_i w_i^{(1)} p_k^2(x_i)\right)\left(\sum_j w_j^{(2)} q_l^2(y_j)\right) = \|p_k\|^2 \|q_l\|^2,$$

en $\langle f, \phi_{k,l} \rangle = \sum_{i,j} w_{i,j} f_{i,j} p_k(x_i) q_l(y_j)$.

Delen geeft de formule. $\square$

> ⚠️ **Examentip:** De sleutelboodschap is dat het normaalstelsel diagonaal wordt door de orthogonale basis → elke coëfficiënt kan apart berekend worden.

---

### (c) Efficiënte berekening van de coëfficiënten

**Sleutelobservatie:** De formule voor $a_{k,l}$ kan gefactoriseerd worden:
$$a_{k,l} = \frac{\displaystyle\sum_{j=1}^N w_j^{(2)} b_{k,j}\, q_l(y_j)}{\|q_l\|^2}, \qquad \text{waarbij} \quad b_{k,j} = \frac{\displaystyle\sum_{i=1}^M w_i^{(1)} f_{i,j}\, p_k(x_i)}{\|p_k\|^2}.$$

⚙️ **Algoritme: herhaalde eendimensionale benadering**

```
Stap 1 — Benadering langs de x-richting:
  voor j = 1 tot N:
      bereken b[k,j] = <f(:,j), p_k> / ||p_k||²   voor k = 0,...,m
      (dit is een 1D Forsythe-benadering van de j-de datarij)

Stap 2 — Benadering langs de y-richting:
  voor k = 0 tot m:
      bereken a[k,l] = <b[k,:], q_l> / ||q_l||²   voor l = 0,...,n
      (dit is een 1D Forsythe-benadering van de k-de coëfficiëntrij)
```

**Grafische illustratie:**
```
Data f[i,j]:           Stap 1:               Stap 2:
M × N rooster   →  b[k,j] (m+1)×N    →   a[k,l] (m+1)×(n+1)
                   elke kolom j           elke rij k
                   Forsythe in x          Forsythe in y
```

---

### (d) Vergelijking rekencomplexiteit

**Directe methode (monomiale basis, vol normaalstelsel):**

Normaalstelsel: $(m+1)(n+1)$ onbekenden, matrix $(m+1)^2(n+1)^2$ entries.
- Opstellen: $O(MN(m+1)(n+1))$
- Oplossen: $O((m+1)^3(n+1)^3)$ — **kubisch in het aantal basisfuncties**

**Methode via herhaalde 1D Forsythe:**

- Stap 1: $N$ keer Forsythe met $M$ punten, graad $m$ → $O(NM(m+1))$
- Stap 2: $(m+1)$ keer Forsythe met $N$ punten, graad $n$ → $O((m+1)N(n+1))$
- **Totaal: $O(MN(m+n))$** — **lineair** in het aantal data- en basisfuncties

| Methode | Rekenkost | Stabiliteit |
|---|---|---|
| Directe methode | $O((mn)^3)$ | Slecht (Hilbert-matrix) |
| Herhaalde 1D Forsythe | $O(MN(m+n))$ | Uitstekend (diagonaal stelsel) |

> ⚠️ **Examentip:** Het sleutelwoord is "herhaalde eendimensionale benadering". Dankzij separabiliteit reduceert 2D → 1D + 1D.

---

## V3 — Discreet benaderingsprobleem: veeltermbenadering

**Probleemstelling:** Gegeven meetdata $(x_i, f_i)_{i=1}^N$, $a \leq x_1 < \cdots < x_N \leq b$. Zoek een veelterm $y_n(x)$ van graad $n \ll N$.

### (a) Scalair product voor discrete kleinste-kwadraten

Het **discrete gewogen scalair product** op de ruimte van discrete functies is:
$$\boxed{(f, g) = \sum_{i=1}^N w_i\, f(x_i)\, g(x_i), \qquad w_i > 0.}$$

De bijhorende norm: $\|f\|^2 = \sum_{i=1}^N w_i f(x_i)^2$.

**Benaderingsprobleem:** Minimaliseer
$$E(c) = \|f - y_n\|^2 = \sum_{i=1}^N w_i (f_i - y_n(x_i))^2.$$

---

### (b) Normaalstelsel = overgedetermineerd stelsel

**Aanpak 1 — Normaalstelsel:**

Schrijf $y_n(x) = \sum_{k=0}^n c_k \phi_k(x)$. De orthogonaliteitsconditie geeft het **normaalstelsel**:
$$Gc = b, \qquad G_{jk} = (\phi_j, \phi_k) = \sum_{i=1}^N w_i \phi_j(x_i) \phi_k(x_i), \quad b_j = (\phi_j, f) = \sum_{i=1}^N w_i \phi_j(x_i) f_i.$$

**Aanpak 2 — Overgedetermineerd stelsel:**

Definieer $D = \text{diag}(\sqrt{w_1}, \ldots, \sqrt{w_N})$ en de matrix $A \in \mathbb{R}^{N \times (n+1)}$ met $A_{ik} = \phi_k(x_i)$.

**Bewijs van de equivalentie:**

Het overgedetermineerde stelsel is $DAc = Df$ (met $N \gg n+1$ vergelijkingen en $n+1$ onbekenden). De normaalvergelijkingen (i.e. $A^TD^TDAc = A^TD^TDf$) zijn:
$$A^T W A c = A^T W f, \qquad W = D^2 = \text{diag}(w_1, \ldots, w_N).$$

Uitschrijven: $(A^TWA)_{jk} = \sum_i w_i A_{ij} A_{ik} = \sum_i w_i \phi_j(x_i) \phi_k(x_i) = G_{jk}$.

Het zijn **dezelfde stelsels**. $\square$

> ⚠️ **Examentip:** Dit verband is het brugpunt tussen de abstracte normaalstelsel-theorie en de concrete matrixformulering. Ken het.

---

### (c) Detecteren van overfitting

**Overfitting:** Het model past te nauw aan bij de ruizige trainingsdata. Het leert de ruis, niet de structuur van $f$. Herkenbaar doordat het goed presteert op de trainingspunten maar slecht op nieuwe punten.

**Detectiemethode — trainings- vs. testfout:**

Splits de data in twee sets:
- **Trainingsdata:** $\{(x_i, f_i)\}_{i=1}^N$ — voor het opstellen van de benadering
- **Testdata:** $\{(\tilde{x}_i, \tilde{f}_i)\}$ — voor het evalueren op ongeziene data

$$E_{\text{train}} = \sum_{i=1}^N w_i (f_i - y_n(x_i))^2, \qquad E_{\text{test}} = \sum_i \tilde{w}_i (\tilde{f}_i - y_n(\tilde{x}_i))^2.$$

**Gedrag als functie van de graad $n$:**
```
Fout
 |  E_test
 |    \   /
 |     \ /   ← optimale graad
 |      *
 |        \_____ E_train
 |_____________________ n
```

- **Trainingsfout** daalt **monotoon** (meer vrijheidsgraden → lagere fout).
- **Testfout** daalt aanvankelijk, maar **stijgt opnieuw** bij te hoge graad → signaal van overfitting.

De **optimale graad** minimaliseert $E_{\text{test}}$.

---

### (d) Tikhonov regularisatie

**Definitie en motivatie:**

In plaats van de complexiteit (graad) op voorhand te begrenzen, voegen we een **strafterm** toe aan de kostfunctie die grote coëfficiënten ontmoedigt:

$$E_{\text{reg}}(c) = \underbrace{\sum_{i=1}^N w_i (f_i - y_n(x_i))^2}_{\text{datafidelity}} + \lambda \underbrace{\|c\|_2^2}_{\text{regularisatie}}, \qquad \lambda > 0.$$

**Motivatie:**
- Grote coëfficiënten $c_k$ geven aanleiding tot sterk oscillerende veeltermen (overfitting).
- Door $\lambda \|c\|_2^2$ toe te voegen, worden grote coëfficiënten bestraft.
- De parameter $\lambda > 0$ regelt de **balans** tussen trouw aan de data en eenvoud van het model.

| $\lambda$ | Effect |
|---|---|
| $\lambda = 0$ | Geen regularisatie — gewone kleinste-kwadraten |
| $\lambda$ klein | Lichte afschrikking van grote coëfficiënten |
| $\lambda \to \infty$ | Alle $c_k \to 0$ (het model vervalt) |

---

### (e) Wijziging aan het normaalstelsel

**Stelling:** Tikhonovregularisatie geeft slechts een **kleine aanpassing** aan het gewone normaalstelsel: de matrix $A^TWA$ wordt vervangen door $A^TWA + \lambda I$.

**Bewijs (afleiding van het geregulariseerde normaalstelsel):**

Schrijf $E_{\text{reg}}(c) = \|D(Ac - f)\|_2^2 + \lambda \|c\|_2^2$. Uitwerking:
$$E_{\text{reg}}(c) = c^T(A^TWA + \lambda I)c - 2f^TWAc + f^TWf.$$

Differenteer naar $c$ en stel gelijk aan nul:
$$\nabla_c E_{\text{reg}} = 2(A^TWA + \lambda I)c - 2A^TWf = 0.$$

Het **geregulariseerde normaalstelsel:**
$$\boxed{(A^TWA + \lambda I)\hat{c} = A^TWf.}$$

**Vergelijking:**

| | Normaalstelsel | Geregulariseerd |
|---|---|---|
| Matrix | $A^TWA$ | $A^TWA + \lambda I$ |
| Inverteerbaar? | Nee (als $A$ niet vol rang) | Ja (voor $\lambda > 0$) |
| Conditiegetal | Hoog | Lager ($\lambda$ vergroot de minimale eigenwaarde) |

> ⚠️ **Examentip:** De aanpassing is enkel het optellen van $\lambda I$ op de diagonaal van de normaalmatrix. Dit is een kleine maar cruciale wijziging.

💡 **Waarom altijd inverteerbaar?** De eigenwaarden van $A^TWA + \lambda I$ zijn $\sigma_i^2 + \lambda \geq \lambda > 0$ voor alle $i$. De matrix is dus altijd strikt positief definiet voor $\lambda > 0$.

---

## V4 — Overfitting en regularisatie

### (a) Overfitting: definitie en detectie

**Overfitting** treedt op wanneer een model te nauw aansluit bij de **ruizige trainingsdata** in plaats van bij de onderliggende functie $f$:
- Het model **leert de ruis** $\varepsilon_i$ in plaats van de structuur van $f$.
- Het presteert goed op de trainingspunten: lage $E_{\text{train}}$.
- Het presteert **slecht op nieuwe punten**: hoge $E_{\text{test}}$.

> ⚠️ **Examentip:** Onderscheid overfitting van onderfitting (oversmoothing). Bij onderfitting zijn zowel $E_{\text{train}}$ als $E_{\text{test}}$ hoog.

**Vaststellen via generalisatiefout:**
- Splits data in trainingsset en testset.
- Evalueer de testfout $E_{\text{test}}$ als functie van de modelcomplexiteit (graad, aantal knooppunten, etc.).
- Overfitting treedt op zodra $E_{\text{test}}$ begint te stijgen terwijl $E_{\text{train}}$ blijft dalen.

**Alternatief:** Leave-one-out kruisvalidatie — stel $N$ modellen op, elk met één punt weggelaten.

---

### (b) Tikhonov vs. LASSO

Beide methodes regulariseren via een strafterm $\lambda F(c)$ in de kostfunctie, maar de keuze van $F$ geeft heel andere oplossingen.

**Tikhonov (L2-regularisatie / ridge regression):**
$$F(c) = \|c\|_2^2 = \sum_k c_k^2 \quad \Rightarrow \quad (A^TWA + \lambda I)\hat{c} = A^TWf.$$

**LASSO (L1-regularisatie):**
$$F(c) = \|c\|_1 = \sum_k |c_k| \quad \Rightarrow \quad \text{niet-lineair stelsel.}$$

**Verschil:**

| Eigenschap | Tikhonov (L2) | LASSO (L1) |
|---|---|---|
| Coëfficiënten | Klein maar **nooit exact nul** | Vaak **exact nul** (schaarse oplossing) |
| Optimalisatieprobleem | Lineair stelsel | Niet-lineair (iteratief oplossen) |
| Differentieerbaar? | Ja | Niet in $c_k = 0$ |
| Selectie van relevante features | Nee | Ja (automatische variabelenselectie) |

**Geometrisch argument (waarom LASSO schaars is):**

Het regularisatieprobleem $\min E(c) + \lambda F(c)$ is equivalent met $\min E(c)$ s.t. $F(c) \leq K$.

```
   L2-bol (Tikhonov)          L1-bol (LASSO)

         ___                       /\
        /   \                     /  \
       | c*  |   c* op rand      / c* \   c* in hoek
        \___/    (geen nulpunt)  \    /   (c₂ = 0)
                                  \  /
                                   \/
```

- **L2-bol** is glad en rond → raakpunt elliptische niveaukrommen ligt willekeurig op de rand → geen coëfficiënt precies nul.
- **L1-bol** heeft scherpe hoeken op de assen → niveaukrommen raken de bol het vaakst in een hoekpunt → één (of meer) coëfficiënten exact nul.

> ⚠️ **Examentip:** Ken dit geometrisch argument. Het staat bijna zeker op het examen.

---

### (c) Regularisatie bij diepe neurale netwerken

**Geregulariseerde kostfunctie voor een neuraal netwerk met $L$ lagen:**

$$\text{Kost}_{\text{reg}}(c) = \frac{1}{N}\sum_{i=1}^N \left\|\mathbf{y}(\mathbf{x}^{(i)}) - \mathbf{a}^{[L]}(\mathbf{x}^{(i)})\right\|_2^2 + \frac{\lambda}{N}\sum_{l=2}^L \|W^{[l]}\|_2^2.$$

**Motivatie:** Grote gewichten $W^{[l]}$ maken neuronen gevoelig voor kleine perturbaties in de invoer — dat is een teken van overfitting. Door grote gewichten te penaliseren, wordt het netwerk robuuster en generaliseert het beter.

**Praktisch:** De gewichten worden "kleiner gehouden", maar de biases $\mathbf{b}^{[l]}$ worden **niet** geregulariseerd (zij bepalen de verschuiving, niet de gevoeligheid).

**Andere technieken:**
- **Dropout:** willekeurig neuronen tijdelijk uitschakelen tijdens het trainen
- **Early stopping:** stop zodra $E_{\text{test}}$ begint te stijgen
- **Data augmentatie:** kunstmatig meer trainingsdata aanmaken (spiegelen, roteren)

---

## V5 — Diep neuraal netwerk

### (a) Neuraal netwerk: 3 lagen, 2 inputs, 3 verborgen, 2 outputs

```
Invoerlaag [1]       Verborgen laag [2]       Uitvoerlaag [3]

    x₁ ○————————→ ○ n₂₁ ——→ ○ n₃₁ → y₁
       |        ↗ ○ n₂₂ ——→ ○ n₃₂ → y₂
       |       ↗  ○ n₂₃
    x₂ ○——————

(volledig verbonden tussen elke laag)
```

**Dimensies:**
- Invoerlaag: $\mathbf{x} \in \mathbb{R}^2$ (2 inputs)
- Verborgen laag [2]: $W^{[2]} \in \mathbb{R}^{3 \times 2}$, $\mathbf{b}^{[2]} \in \mathbb{R}^3$ → **9 parameters**
- Uitvoerlaag [3]: $W^{[3]} \in \mathbb{R}^{2 \times 3}$, $\mathbf{b}^{[3]} \in \mathbb{R}^2$ → **8 parameters**
- **Totaal: 17 parameters**

---

### (b) Voorwaartse propagatie

**Stap 1 — Verborgen laag:**
$$\mathbf{z}^{[2]} = W^{[2]} \mathbf{x} + \mathbf{b}^{[2]} \in \mathbb{R}^3$$
$$\mathbf{a}^{[2]} = \sigma(\mathbf{z}^{[2]}) \in \mathbb{R}^3 \quad \text{(componentgewijs)}$$

**Stap 2 — Uitvoerlaag:**
$$\mathbf{z}^{[3]} = W^{[3]} \mathbf{a}^{[2]} + \mathbf{b}^{[3]} \in \mathbb{R}^2$$
$$\mathbf{a}^{[3]} = \sigma(\mathbf{z}^{[3]}) = F(\mathbf{x}) \in \mathbb{R}^2$$

**Dimensies expliciet:**
$$\underbrace{W^{[2]}}_{3 \times 2} \cdot \underbrace{\mathbf{x}}_{2 \times 1} + \underbrace{\mathbf{b}^{[2]}}_{3 \times 1} = \underbrace{\mathbf{z}^{[2]}}_{3 \times 1} \xrightarrow{\sigma} \underbrace{\mathbf{a}^{[2]}}_{3 \times 1}$$

$$\underbrace{W^{[3]}}_{2 \times 3} \cdot \underbrace{\mathbf{a}^{[2]}}_{3 \times 1} + \underbrace{\mathbf{b}^{[3]}}_{2 \times 1} = \underbrace{\mathbf{z}^{[3]}}_{2 \times 1} \xrightarrow{\sigma} \underbrace{\mathbf{a}^{[3]}}_{2 \times 1}$$

**Volledige samengestelde formule:**
$$F(\mathbf{x}) = \sigma\!\Bigl(W^{[3]} \sigma\!\bigl(W^{[2]} \mathbf{x} + \mathbf{b}^{[2]}\bigr) + \mathbf{b}^{[3]}\Bigr).$$

> ⚠️ **Examentip:** Op het examen moet je alle dimensies expliciet benoemen. Schrijf ze bij elke stap.

---

### (c) Trainingsalgoritme bij veel datapunten

**Probleem:** Bij grote $N$ is het berekenen van de volledige gradiënt $\nabla \text{Kost}(c) = \frac{1}{N}\sum_{i=1}^N \nabla C_{\mathbf{x}^{(i)}}(c)$ te duur.

**Aanbevolen methode: Mini-batch stochastische gradiëntafdaling (mini-batch SGD)**

⚙️ **Algoritme:**

```
Gegeven: trainingsdata, learning rate α, batchgrootte m
voor elke epoch:
    permuteer de N trainingspunten willekeurig
    voor k = 1, 2, ..., N/m:
        kies mini-batch {x^{k₁}, ..., x^{kₘ}} (m willekeurige punten)
        bereken gradiënt via backpropagation:
            g = (1/m) Σᵢ ∇C_{x^{kᵢ}}(c^(k))
        update: c^(k+1) = c^(k) - α · g
```

**Hyperparameters die te kiezen zijn:**
- **Learning rate $\alpha > 0$:** beïnvloedt de stapgrootte; te groot → divergentie, te klein → trage convergentie
- **Batchgrootte $m$:** compromis tussen nauwkeurigheid van de gradiëntschatting ($m$ groot) en snelheid per stap ($m$ klein)
- **Aantal epochs:** hoeveel keer de volledige dataset doorlopen wordt
- **Regularisatieparameter $\lambda$:** voor L2-regularisatie

**Voordelen van mini-batch SGD vs. volledige gradiëntafdaling:**
- Goedkoper per stap: één matrix-vectorproduct per mini-batch ipv $N$
- Kan beter ontsnappen aan lokale minima (de ruis in de gradiëntschatting helpt)
- Geschikt voor GPU's (parallelle berekening per mini-batch)

> ⚠️ **Examentip:** De gradiënt wordt berekend via **backpropagation** (de kettingregel achterwaarts doorheen de lagen). Ken de terminologie: epoch, mini-batch, learning rate, backpropagation.

---

## V6 — Methode van de toegevoegde gradiënten (Conjugate Gradients)

**Probleemstelling:** Los op $Ax = b$ met $A \in \mathbb{R}^{m \times m}$ **symmetrisch positief definiet** (SPD).

### (a) Algoritme in pseudocode + uitleg

⚙️ **Algoritme 12 — Conjugate Gradients (standaard):**

```
Invoer: A (SPD), b, kmax, ε
p⁽⁰⁾ = r⁽⁰⁾ = b,   x⁽⁰⁾ = 0
voor k = 1, 2, ..., kmax:
    α⁽ᵏ⁾ = (r⁽ᵏ⁻¹⁾ᵀ r⁽ᵏ⁻¹⁾) / (p⁽ᵏ⁻¹⁾ᵀ A p⁽ᵏ⁻¹⁾)
    x⁽ᵏ⁾ = x⁽ᵏ⁻¹⁾ + α⁽ᵏ⁾ · p⁽ᵏ⁻¹⁾
    r⁽ᵏ⁾ = r⁽ᵏ⁻¹⁾ - α⁽ᵏ⁾ · A p⁽ᵏ⁻¹⁾
    als ‖r⁽ᵏ⁾‖₂ < ε: geef x⁽ᵏ⁾ terug
    β⁽ᵏ⁾ = (r⁽ᵏ⁾ᵀ r⁽ᵏ⁾) / (r⁽ᵏ⁻¹⁾ᵀ r⁽ᵏ⁻¹⁾)
    p⁽ᵏ⁾ = r⁽ᵏ⁾ + β⁽ᵏ⁾ · p⁽ᵏ⁻¹⁾
```

**Uitleg van elke stap:**

- **$r^{(k)} = b - Ax^{(k)}$:** het **residu** (= de fout op de vergelijking). Gelijk aan de negatieve gradiënt van $f(x) = \frac{1}{2}x^TAx - x^Tb$.
- **$\alpha^{(k)}$:** de **optimale stapgrootte** langs de zoekrichting $p^{(k-1)}$, gevonden via exacte lijnzoeking.
- **$p^{(k)} = r^{(k)} + \beta^{(k)} p^{(k-1)}$:** de nieuwe **zoekrichting** is een lineaire combinatie van het huidige residu en de vorige zoekrichting.
- **$\beta^{(k)}$:** de **conjugeringsfactor** die de zoekrichtingen A-toegevoegd maakt.
- **Stopcriterium:** $\|r^{(k)}\|_2 < \varepsilon$.

**Verband met optimalisatie:** $Ax = b$ is gelijkwaardig met het minimaliseren van $f(x) = \frac{1}{2}x^TAx - x^Tb$ (want $\nabla f = Ax - b = 0 \Leftrightarrow Ax = b$).

---

### (b) De methode bouwt een Krylov-ruimte op

**Stelling (Eigenschap 11.2):**

Na $k-1$ stappen geldt:
$$\mathcal{K}_k(A, b) = \langle b, Ab, \ldots, A^{k-1}b \rangle = \langle p^{(0)}, \ldots, p^{(k-1)} \rangle = \langle r^{(0)}, \ldots, r^{(k-1)} \rangle.$$

**Bewijs (door inductie):**

**Basisgeval $k=1$:** $p^{(0)} = r^{(0)} = b$, zodat $\mathcal{K}_1(A,b) = \langle b \rangle = \langle p^{(0)} \rangle = \langle r^{(0)} \rangle$. ✓

**Inductiestap:** Veronderstel dat de stelling geldt voor stap $k-1$, d.w.z.
$$\langle p^{(0)}, \ldots, p^{(k-2)} \rangle = \langle r^{(0)}, \ldots, r^{(k-2)} \rangle = \mathcal{K}_{k-1}(A,b).$$

Uit de update $r^{(k-1)} = r^{(k-2)} - \alpha^{(k-1)} A p^{(k-2)}$ volgt:
$$r^{(k-1)} \in r^{(k-2)} + A\, \langle p^{(k-2)} \rangle \subset \mathcal{K}_{k-1}(A,b) + A\, \mathcal{K}_{k-1}(A,b) = \mathcal{K}_k(A,b).$$

Uit de update $p^{(k-1)} = r^{(k-1)} + \beta^{(k-1)} p^{(k-2)}$ volgt $p^{(k-1)} \in \mathcal{K}_k(A,b)$.

Omgekeerd: $A^{k-1}b \in \mathcal{K}_k(A,b)$ is een lineaire combinatie van $r^{(0)}, \ldots, r^{(k-1)}$ (via de inductiehypothese uitgebreid). $\square$

> ⚠️ **Examentip:** De sleutel is dat $r^{(k)} = r^{(k-1)} - \alpha^{(k)} A p^{(k-1)}$: vermenigvuldiging met $A$ verhoogt de Krylov-graad met 1.

---

### (c) CG vindt de beste benadering in de Krylov-ruimte

**Stelling 11.3.3 — A-norm minimaliteit:**

Definieer de **fout** $e^{(k)} = x^{(k)} - x^*$ en de **A-norm** $\|x\|_A = \sqrt{x^T A x}$. Dan is $x^{(k)}$ het element van $\mathcal{K}_k(A,b)$ dat de A-norm van de fout minimaliseert:
$$\|e^{(k)}\|_A^2 = \min_{x \in \mathcal{K}_k(A,b)} \|x - x^*\|_A^2.$$

**Bewijs:**

Elke $x \in \mathcal{K}_k(A,b)$ schrijven we als $x = x^{(k)} + \delta x$ met $\delta x \in \mathcal{K}_k(A,b)$:
$$\|x - x^*\|_A^2 = \|\delta x + e^{(k)}\|_A^2 = \|e^{(k)}\|_A^2 + 2\underbrace{e^{(k)T} A\, \delta x}_{?} + \|\delta x\|_A^2.$$

Nu is $Ae^{(k)} = A(x^{(k)} - x^*) = Ax^{(k)} - b = -r^{(k)}$.

Uit Eigenschap 11.3 (orthogonaliteit van residu's) geldt: $r^{(k)} \perp \mathcal{K}_k(A,b)$, d.w.z. $r^{(k)T} \delta x = 0$ voor alle $\delta x \in \mathcal{K}_k(A,b)$.

Dus:
$$e^{(k)T} A \delta x = -r^{(k)T} \delta x = 0.$$

De gekruiste term is altijd nul, en $\|\delta x\|_A^2 \geq 0$. Het minimum wordt bereikt voor $\delta x = 0$, d.w.z. $x = x^{(k)}$. $\square$

---

### (d) Verband met veeltermbenadering + versnelling

**Stelling 11.3.4:**

$$\|e^{(k)}\|_A^2 = \min_{p \in \mathcal{N}_k} \|p(A) e^{(0)}\|_A^2,$$

waarbij $\mathcal{N}_k = \{p \in \mathbb{R}[t],\ \deg p \leq k,\ p(0) = 1\}$.

**Bewijs (schets):** Elke $x \in \mathcal{K}_k(A,b)$ geeft $x - x^* = \hat{p}(A) e^{(0)}$ voor een $\hat{p} \in \mathcal{N}_k$ (want $b = Ax^* = -Ae^{(0)}$ en Krylov-vectoren zijn veeltermen in $A$ toegepast op $b$). Het gestelde volgt uit de A-norm-minimaliteitsresultaat.

**Interpretatie:** Via de eigenwaardenontbinding $A = X\Lambda X^{-1}$:
$$\|e^{(k)}\|_A^2 \approx \min_{p \in \mathcal{N}_k} \sum_{i=1}^m |p(\lambda_i)|^2 \|x_i^* \text{-component}\|_A^2.$$

De convergentie hangt af van hoe goed het spectrum van $A$ **afgedekt** kan worden door een veelterm $p \in \mathcal{N}_k$ die klein is op alle eigenwaarden.

**Implicaties voor convergentie:**
- Als $A$ slechts **$k$ verschillende eigenwaarden** heeft, convergeert CG exact in $k$ stappen.
- Als eigenwaarden **geclusterd** zijn, is snelle convergentie te verwachten.

**Hoe exploiteren? — Preconditioning:**

Kies een inverteerbare matrix $D$ zodat $D^{-1}AD^{-T}$ beter geclusterde eigenwaarden heeft. Los in plaats van $Ax = b$ het equivalente stelsel op:
$$(D^{-1}AD^{-T})\underbrace{D^Tx}_{y} = D^{-1}b.$$

Pas CG toe op dit **geconditioneerde stelsel**. Als $D^{-1}AD^{-T}$ eigenwaarden dicht bij 1 heeft, convergeert CG in weinig iteraties.

> ⚠️ **Examentip:** Preconditioning is de praktische toepassing van het veeltermbenaderingsverband. Ken de boodschap: cluster de eigenwaarden → snellere convergentie.

---

## V7 — Niet-lineair optimalisatieprobleem

**Probleemstelling:** Los $\min_{x \in \mathbb{R}^n} f(x)$ op met $f$ tweemaal continu differentieerbaar.

### (a) Sterkste daling in de richting van de negatieve gradiënt

**Claim:** Van alle richtingen $s$ met $\|s\|_2 = 1$ daalt $f$ het sterkst in de richting $-\nabla f(x)$.

**Bewijs:**

De **richtingsafgeleide** van $f$ in punt $x$ in richting $s$ (met $\|s\|_2 = 1$) is:
$$\frac{d}{d\alpha} f(x + \alpha s)\bigg|_{\alpha=0} = s^T \nabla f(x).$$

We zoeken de richting $s$ waarvoor dit het **meest negatief** is, d.w.z.:
$$\min_{\|s\|_2 = 1} s^T \nabla f(x).$$

Via de **Cauchy-Schwarz-ongelijkheid**:
$$s^T \nabla f(x) \geq -\|s\|_2 \cdot \|\nabla f(x)\|_2 = -\|\nabla f(x)\|_2,$$

met gelijkheid wanneer:
$$s = -\frac{\nabla f(x)}{\|\nabla f(x)\|_2}.$$

**Conclusie:** De richtingsafgeleide is het meest negatief (= sterkste daling) in de richting $s = -\nabla f(x) / \|\nabla f(x)\|_2$, d.w.z. de **negatieve gradiëntrichting**. $\square$

> ⚠️ **Examentip:** Dit is Eigenschap 11.1. Ken het bewijs: richtingsafgeleide + Cauchy-Schwarz.

💡 **Intuïtie:** De gradiënt $\nabla f(x)$ staat loodrecht op de niveauverzameling door $x$. Door in de richting $-\nabla f(x)$ te gaan, beweeg je loodrecht door de niveauverzamelingen, d.w.z. "het steilst omlaag".

---

### (b) Lijnzoekmethode vs. Armijo-backtracking

Beide methodes dienen voor het bepalen van de **stapgrootte** $\alpha^{(k)}$ in de iteratie:
$$x^{(k)} = x^{(k-1)} + \alpha^{(k)} p^{(k)}, \qquad p^{(k)} = -\nabla f(x^{(k-1)}).$$

#### Exacte lijnzoekmethode

**Definitie:**
$$\alpha^{(k)} = \arg\min_{\alpha > 0} f(x^{(k-1)} + \alpha p^{(k)}).$$

De stapgrootte minimaliseert de doelfunctie **exact** langs de gekozen richting $p^{(k)}$.

**Eigenschappen:**
- **Optimaal langs de richting:** men maakt de best mogelijke stap in de richting $p^{(k)}$
- **Duur:** vereist een bijkomend (scalair) optimalisatieprobleem per iteratie
- **Zigzag-convergentie bij steilste afdaling:** opeenvolgende richtingen staan loodrecht op elkaar → trage convergentie

**Bewijs van loodrechtheid:** Na exacte lijnzoeking is $\frac{d}{d\alpha}f(x^{(k)}) = \nabla f(x^{(k)})^T p^{(k)} = 0$. Maar $p^{(k+1)} = -\nabla f(x^{(k)})$, dus $p^{(k+1)} \perp p^{(k)}$ — elke nieuwe richting staat loodrecht op de vorige.

#### Backtracking met de Armijo-voorwaarde

**Idee:** Geen exacte minimalisatie langs de richting, maar slechts een **voldoende daling** eisen.

**Armijo-voorwaarde:**
$$f\!\left(x^{(k-1)} + \alpha p^{(k)}\right) \leq f\!\left(x^{(k-1)}\right) + c \cdot \alpha \cdot \nabla f(x^{(k-1)})^T p^{(k)}, \qquad c \in (0, 1).$$

⚙️ **Algoritme (backtracking):**
```
α = 1  (startwaarde)
terwijl Armijo-voorwaarde niet voldaan:
    α = q · α   (typisch q = 0.5)
```

**Intuïtie:**
```
f
|  f(x) + α·∇f·p  (raaklijn: helling ∇f·p < 0)
|   \
|    \  f(x) + c·α·∇f·p  (minder steil, c < 1)
|     \___
|         \
|      f(x + αp)   ← moet hier onder liggen
|___________________________ α
```

De Armijo-voorwaarde eist dat de nieuwe functiewaarde **onder de gecorrigeerde raaklijn** ligt, wat een voldoende daling garandeert.

**Drie geometrische eigenschappen die Samaey expliciet test:**

1. **Gegarandeerde stop:** De eis is *minder streng* dan de raaklijn van de doelfunctie (die helling $\nabla f \cdot p$ heeft). Omdat $f$ continu is en de raaklijn neerwaarts gaat, is er altijd een $\alpha > 0$ klein genoeg zodat de functiewaarde onder de gecorrigeerde raaklijn (helling $c \cdot \nabla f \cdot p$) valt. Backtracking stopt dus altijd in eindig aantal stappen.

2. **Strenger voor grotere stappen:** Het rechterlid $f(x) + c\alpha \nabla f \cdot p$ groeit lineair in $\alpha$ (daalt lineair, want $\nabla f \cdot p < 0$). Voor grote $\alpha$ is de eis dat $f(x + \alpha p)$ ver genoeg onder $f(x)$ ligt — strenger dan voor kleine $\alpha$. Dit vermijdt te grote stappen die over het minimum springen.

3. **Groot aanvaardingsgebied ver van het optimum:** Ver van het optimum is $|\nabla f \cdot p|$ groot (steile helling), zodat de gecorrigeerde raaklijn snel daalt. Er is een groot bereik aan $\alpha$-waarden die de Armijo-eis voldoen → veel minder halveerstappe nodig → rekenwerk bespaart. Dicht bij het optimum is de daling klein en worden de stappen kleiner, wat verhindert dat het algoritme stopt.

**Vergelijking:**

| | Exacte lijnzoeking | Armijo-backtracking |
|---|---|---|
| Rekenkost | Hoog (1 extra minimaliseringsprobleem) | Laag (slechts functie-evaluaties) |
| Garantie | Optimale stap langs richting | Voldoende daling |
| Convergentiesnelheid | Theoretisch optimaal | Praktisch snel genoeg |
| Zigzag | Ja (bij steilste afdaling) | Minder uitgesproken |
| Gebruik | Theoretische analyse | Praktische implementaties |

> ⚠️ **Examentip:** Ken de Armijo-voorwaarde letterlijk en weet de intuïtie uitleggen (gecorrigeerde raaklijn). Weet ook het verschil met de exacte lijnzoeking.

💡 **Waarom $c < 1$?** De raaklijn $f(x) + \alpha \nabla f \cdot p$ (met helling $\nabla f \cdot p < 0$) heeft een neerwaartse helling. De Armijo-conditie eist dat we minstens een fractie $c$ van die helling omlaag gaan — vermijdt te kleine stappen (stagnatie) en garandeert monotoniciteit $f(x^{(k)}) < f(x^{(k-1)})$.

---

## Samenvatting Sleutelformules

| Concept | Formule |
|---|---|
| Dimensie splineruimte | $n + k$ |
| B-splinerecursie | $N_{i,k+1} = \frac{x-t_i}{t_{i+k}-t_i}N_{i,k} + \frac{t_{i+k+1}-x}{t_{i+k+1}-t_{i+1}}N_{i+1,k}$ |
| Partitie van eenheid | $\sum_i N_{i,k+1}(x) = 1$ |
| 2D-orthogonale coëfficiënten | $a_{k,l} = \frac{\langle f, p_k q_l \rangle}{\|p_k\|^2 \|q_l\|^2}$ |
| Tikhonov normaalstelsel | $(A^TWA + \lambda I)\hat{c} = A^TWf$ |
| Voorwaartse propagatie | $\mathbf{a}^{[l]} = \sigma(W^{[l]}\mathbf{a}^{[l-1]} + \mathbf{b}^{[l]})$ |
| CG-stapgrootte | $\alpha^{(k)} = \frac{r^{(k-1)T}r^{(k-1)}}{p^{(k-1)T}Ap^{(k-1)}}$ |
| CG-conjugering | $\beta^{(k)} = \frac{r^{(k)T}r^{(k)}}{r^{(k-1)T}r^{(k-1)}}$ |
| Richtingsafgeleide | $\nabla_s f(x) = s^T \nabla f(x)$ |
| Armijo-voorwaarde | $f(x + \alpha p) \leq f(x) + c\alpha\nabla f^T p$, $c \in (0,1)$ |
