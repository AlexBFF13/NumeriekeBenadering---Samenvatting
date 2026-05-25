# Deel 6 – Discrete benadering op basis van meetdata

## Overzicht en rode draad

Tot nu toe benaderden we **gekende** functies door veeltermen of splines. In de praktijk heeft men echter geen gesloten formule voor $f$: men heeft enkel een stel **meetpunten** $(x_i, f_i)$, waarbij de $f_i$ meetfouten kunnen bevatten. We willen dan een continue functie vinden die "goed genoeg" door die punten loopt — niet exact interpoleren (dat zou de meetfouten mee-interpoleren), maar een gladde, representatieve benadering.

Dit hoofdstuk brengt de theorie van Hoofdstukken 3–5 samen in drie concrete situaties:
1. **Discrete veeltermbenadering** (§6.1): kleinste-kwadraten met veeltermen, methode van Forsythe.
2. **Discrete splinebenadering** (§6.2): kleinste-kwadraten met B-splines, normaalstelsel en QR-factorisatie.
3. **Tweedimensionale benadering** (§6.3): benadering op een rechthoekig rooster, tensorproduct-splines, herhaalde eendimensionale benadering.

---

## 6.1 Discrete veeltermbenadering

### De discrete unitaire ruimte

Gegeven een stel abscissen $a \leq x_1 < x_2 < \cdots < x_N \leq b$ met bijhorende functiewaarden $f_1, f_2, \ldots, f_N$ (gemeten, dus met fouten). De ruimte van discrete functies $\ell^2(N)$ is een unitaire ruimte met het **discrete gewogen scalair product**:

$$
(f, g) = \sum_{i=1}^N w_i\, f_i\, g_i, \qquad w_i > 0.
$$

**Isomorfisme:** $\ell^2(N)$ is isomorf met $\mathbb{R}^N$ (gewogen Euclidisch scalair product) én met $P_{N-1}[a,b]$ (veeltermen van graad $\leq N-1$). Met elk punt $(x_1,f_1), \ldots, (x_N, f_N)$ correspondeert namelijk een unieke interpolerende veelterm van graad $\leq N-1$.

### Het benaderingsprobleem

Zoek een veelterm $y_n(x)$ van graad $n \ll N$ die de uitdrukking

$$
\sum_{i=1}^N w_i\, (f_i - y_n(x_i))^2
$$

minimaliseert. Dit is de orthogonale projectie van de discrete functie $(f_1, \ldots, f_N)$ op de deelruimte van veeltermen van graad $\leq n$, in de norm van $\ell^2(N)$.

---

### Het normaalstelsel: drie keuzen voor de basis

**Keuze 1 — Monomiale basis** $\{1, x, x^2, \ldots, x^n\}$:

$$
y_n(x) = c_0 + c_1 x + \cdots + c_n x^n.
$$

Het normaalstelsel $Gc = b$ heeft $G_{jk} = \sum_{i=1}^N w_i x_i^{j+k}$.

**Probleem:** Voor punten op een beperkt interval is dit opnieuw een bijna-Hilbert-matrix, met hoog conditiegetal. Bovendien: bij het verhogen van de graad van $n$ naar $n+1$ moet het hele stelsel worden herbouwd.

**Keuze 2 — Continue orthogonale veeltermen** (bijv. Legendre $P_k$):

$$
y_n(x) = b_0 P_0(x) + b_1 P_1(x) + \cdots + b_n P_n(x).
$$

De grammatrix heeft $G_{kl} = \sum_{i=1}^N w_i P_k(x_i) P_l(x_i)$.

Dit is **geen** diagonaalmatrix: de Legendre-veeltermen zijn orthogonaal voor het **continue** scalair product $(f,g) = \int_{-1}^1 fg\, dx$, maar niet voor het **discrete** scalair product. Voor regelmatig verspreide punten is de matrix wel sterk diagonaaldominant → tamelijk goed geconditioneerd.

**Keuze 3 — Discrete orthogonale veeltermen** (methode van Forsythe):

Zoek veeltermen $p_0(x), p_1(x), \ldots$ die orthogonaal zijn t.o.v. **sommatie**:

$$
\sum_{i=1}^N w_i\, p_k(x_i)\, p_l(x_i) = 0 \quad \text{voor } l \neq k.
$$

Dan is de grammatrix **exact diagonaal** en worden de coëfficiënten direct gegeven door:

$$
\boxed{y_n(x) = \sum_{k=0}^n a_k p_k(x), \quad a_k = \frac{\sum_{i=1}^N w_i\, p_k(x_i)\, f_i}{\sum_{i=1}^N w_i\, p_k^2(x_i)}.}
$$

Dit is de **methode van Forsythe** (1957).

---

### Het overgedetermineerd stelsel

Een alternatieve formulering: zoek $c = (c_0, \ldots, c_n)^T$ die $\|D(Ac - f)\|_2$ minimaliseert, met $A$ de Vandermonde-matrix $A_{ij} = x_i^j$ en $D = \text{diag}(\sqrt{w_1}, \ldots, \sqrt{w_N})$.

Dit is equivalent met het normaalstelsel:

$$
A^T W A\, c = A^T W\, f, \quad W = \text{diag}(w_1, \ldots, w_N).
$$

We hebben de methodes voor het oplossen hiervan al gezien in §2.4 (QR-factorisatie, normaalvergelijkingen).

---

### 6.1.4 De methode van Forsythe in detail

#### Bepalen van de discrete orthogonale veeltermen

Doordat de discrete veeltermen $p_k$ aan dezelfde eigenschappen voldoen als de continue orthogonale veeltermen uit §4.2, geldt ook hier een **drietermsrecursiebetrekking**:

$$
p_{-1}(x) \equiv 0, \quad p_0(x) \equiv \lambda_0,
$$

$$
p_k(x) = \lambda_k\left[(x - \alpha_k)\, p_{k-1}(x) - \beta_k\, p_{k-2}(x)\right] \quad k = 1, 2, \ldots
$$

met **discrete** recursiecoëfficiënten:

$$
\alpha_k = \frac{\sum_{i=1}^N w_i\, x_i\, p_{k-1}^2(x_i)}{\sum_{i=1}^N w_i\, p_{k-1}^2(x_i)}, \qquad \beta_k = \frac{\sum_{i=1}^N w_i\, x_i\, p_{k-1}(x_i)\, p_{k-2}(x_i)}{\sum_{i=1}^N w_i\, p_{k-2}^2(x_i)}.
$$

Per stap $k$ zijn slechts **drie sommen** nodig (de noemer van $\alpha_k$ is gelijk aan de teller van $a_k$).

#### Bewaren van de veeltermwaarden, niet de coëfficiënten

**Idee van Forsythe:** Sla niet de algebraïsche coëfficiënten $b_{kj}$ van $p_k(x) = \sum_j b_{kj} x^j$ op, maar enkel de **functiewaarden** $p_k(x_i)$ voor $i = 1, \ldots, N$.

**Waarom:** De coëfficiënten $b_{kj}$ kunnen enorm groot worden en van wisselend teken, wat numerieke instabiliteit geeft. De functiewaarden $p_k(x_i)$ zijn stabiel en bevatten alle benodigde informatie.

**Berekening:** Via de recursiebetrekking:
$$
p_k(x_i) = \lambda_k\left[(x_i - \alpha_k)\, p_{k-1}(x_i) - \beta_k\, p_{k-2}(x_i)\right], \quad i = 1, \ldots, N.
$$

Zodra $p_k(x_i)$ berekend is, kan $p_{k-2}(x_i)$ worden overschreven — slechts twee kolommen van $N$ waarden nodig.

#### Numeriek stabiele berekening van de coëfficiënten $a_k$

De naïeve formule $a_k = \frac{\sum w_i p_k(x_i) f_i}{\sum w_i p_k^2(x_i)}$ kan numeriek onstabiel zijn voor grote $k$ (grote $p_k$-waarden). Beter is het **residu**-formule:

$$
a_k = \frac{\sum_{i=1}^N w_i\, p_k(x_i)\, r_{k-1,i}}{\sum_{i=1}^N w_i\, p_k^2(x_i)},
$$

waarbij $r_{k-1,i} = f_i - \sum_{j=0}^{k-1} a_j p_j(x_i)$ het residu is van de benadering van graad $k-1$. De residuen $r_{k-1,i}$ nemen af naarmate $k$ toeneemt, wat de foutvoortplanting sterk beperkt.

#### Keuze van de gewichten

De gewichten $w_i$ kunnen statistisch gemotiveerd worden. Als men de functiewaarde $f_i$ meerdere keren gemeten heeft met spreiding $s_i$, dan is de statistisch optimale keuze:

$$
w_i = \frac{1}{s_i^2}.
$$

Andere keuzen:
- $w_i = 1$: alle punten even zwaar (geen a priori kennis van fouten).
- $w_i = (1 - x_i^2)^{-1/2}$: meer gewicht aan de randen (Chebyshev-gewicht → gelijkmatiger fout).
- $w_i = w_i / f_i^2$: minimaliseer **relatieve** fout $\sum w_i \left(\frac{f_i - y_n(x_i)}{f_i}\right)^2$ (mits $f_i \neq 0$).

---

## 6.2 Discrete benadering met splinefuncties

### Probleemstelling

Gegeven meetpunten $(x_r, f_r)$, $r = 1, \ldots, N$. Zoek een spline $s(x)$ van graad $k$ met gegeven knooppunten $t_0, t_1, \ldots, t_n$ zodat

$$
\sum_{r=1}^N w_r\, (f_r - s(x_r))^2 \quad \text{minimaal.}
$$

**Te kiezen parameters:**
- **Graad** van de spline: in de praktijk bijna altijd $k = 3$ (kubisch). Voor hogere afgeleiden: graad 5.
- **Aantal en ligging** van de knooppunten: dit is de moeilijkste keuze.
- **Coëfficiënten** in de B-splinevoorstelling $s(x) = \sum_{i=-k}^{n-1} c_i N_{i,k+1}(x)$.

---

### Het normaalstelsel via B-splines

In B-splinevoorstelling wordt het minimalisatieprobleem het normaalstelsel $Gc = b$ met:

$$
G_{i,j} = \sum_{r=1}^N w_r\, N_{i,k+1}(x_r)\, N_{j,k+1}(x_r), \qquad b_i = \sum_{r=1}^N w_r\, N_{i,k+1}(x_r)\, f_r.
$$

**Dimensie:** $(n+k) \times (n+k)$.

#### Schoenberg-Whitney-voorwaarden

> ⚠️ **Belangrijk voor examen:** Regulariteitscriterium voor de normaalmatrix.

De normaalmatrix $G$ is **symmetrisch positief definiet** (en dus inverteerbaar) als en slechts als er een deelrij meetpunten $u_{-k} < u_{-k+1} < \cdots < u_{n-1}$ bestaat zodat:

$$
t_j < u_j < t_{j+k+1}, \qquad j = -k, -k+1, \ldots, n-1.
$$

**Intuïtie:** In het dragergebied $[t_j, t_{j+k+1}]$ van elke B-spline $N_{j,k+1}$ moet minstens één meetpunt liggen. Anders is de B-spline niet "zichtbaar" in de data en is de corresponderende coëfficiënt onbepaald.

#### Bandstructuur van de normaalmatrix

Door de **lokale-drager-eigenschap** van B-splines:

$$
G_{i,j} = 0 \quad \text{als } |i - j| > k.
$$

De normaalmatrix heeft bandbreedte $2k+1$ — een **bandmatrix**. Dit is een grote numerieke voordeel: het stelsel kan efficiënt worden opgelost met de **Cholesky-factorisatie voor symmetrische bandmatrices**.

---

### Het overgedetermineerd stelsel en QR-factorisatie

Het normaalstelsel kan ook geschreven worden als $A^T W A\, c = A^T W\, f$ waarbij $A_{r,i} = N_{i,k+1}(x_r)$ een $N \times (n+k)$ matrix is.

**Bandstructuur van $A$:** Rij $r$ van matrix $DA$ heeft ten hoogste $k+1$ niet-nul elementen (dankzij de lokale drager van B-splines), in opeenvolgende kolommen.

**Numeriek stabiele aanpak:** QR-factorisatie van $DA$ via **Givens-rotaties**, zonder het expliciete uitrekenen van de normaalmatrix. Na factorisatie $DA = QR$:

$$
\tilde{R}\, c = Q^T D f,
$$

waarbij $\tilde{R}$ een bovendriehoeksmatrix is met bovenbandbreedte $k+1$.

**Voordeel t.o.v. het normaalstelsel:** Het conditioneren van de QR-aanpak is beter (conditiegetal van $A$ in plaats van $A^T A$), en de bandstructuur wordt volledig benut.

---

### Praktische opmerkingen over knooppuntenkeuze

- **Te weinig knooppunten:** de spline sluit onvoldoende aan bij de data.
- **Te veel knooppunten:** de spline volgt de meetfouten op — **overfitting**.
- **Samenvallende knooppunten:** nuttig als de te benaderen functie zelf discontinuïteiten in afgeleiden heeft. Een $l$-voudig knooppunt verlaagt de continuïteitsgraad naar $C^{k-l}$.
- **Automatische knooppuntenbepaling:** men kan ook de knooppuntposities $t_j$ als vrije parameters behandelen. Dan wordt het benaderingsprobleem **niet-lineair** — veel moeilijker, met veel lokale extrema. Goede startwaarden zijn essentieel.
- **Gegevensreductie:** het aantal knooppunten is doorgaans veel kleiner dan het aantal meetpunten, wat leidt tot een compacte representatie.

---

## 6.3 Discrete benadering in twee ruimtelijke dimensies

### Het benaderingsprobleem

Gegeven functiewaarden $f_{i,j}$ op een **regelmatig rechthoekig rooster** $(x_i, y_j)$, $i = 1, \ldots, M$, $j = 1, \ldots, N$. Zoek een tweedimensionale benadering $z(x,y)$ die

$$
\sum_{i=1}^M \sum_{j=1}^N w_{i,j}\, (f_{i,j} - z(x_i, y_j))^2
$$

minimaliseert.

**Separabele gewichtsfunctie:** We veronderstellen $w_{i,j} = w_i^{(1)} w_j^{(2)}$. Dit is cruciaal voor de efficiënte aanpak die volgt.

---

### 6.3.1 Tweedimensionale veeltermbenadering

Schrijf $z(x,y) = \sum_{k=0}^m \sum_{l=0}^n c_{k,l} x^k y^l$ (monomiale basis). Het normaalstelsel heeft $(m+1)(n+1)$ onbekenden en een slecht geconditioneerde coëfficiëntenmatrix.

**Oplossing:** Gebruik **orthogonale veeltermen in twee veranderlijken**.

---

### 6.3.2 Orthogonale veeltermen in twee veranderlijken

De gezochte tweedimensionale veeltermen $\phi_{k,l}(x,y)$ moeten voldoen aan:

$$
\sum_{i=1}^M \sum_{j=1}^N w_{i,j}\, \phi_{k,l}(x_i, y_j)\, \phi_{u,v}(x_i, y_j) = 0 \quad \text{tenzij } k=u \text{ en } l=v.
$$

**Sleutelobservatie:** Door de separabiliteit van het gewicht werkt het **tensorproduct**:

$$
\phi_{k,l}(x, y) = p_k(x)\, q_l(y),
$$

waarbij $\{p_k\}$ orthogonaal zijn t.o.v. sommatie in $x$ en $\{q_l\}$ t.o.v. sommatie in $y$:

$$
\sum_{i=1}^M w_i^{(1)}\, p_k(x_i)\, p_u(x_i) = 0 \text{ voor } k \neq u, \qquad \sum_{j=1}^N w_j^{(2)}\, q_l(y_j)\, q_v(y_j) = 0 \text{ voor } l \neq v.
$$

**Verificatie van de 2D-orthogonaliteit:**

$$
\sum_{i,j} w_{i,j}\, p_k(x_i) q_l(y_j)\, p_u(x_i) q_v(y_j) = \underbrace{\left(\sum_i w_i^{(1)} p_k(x_i) p_u(x_i)\right)}_{=0 \text{ als } k \neq u} \cdot \underbrace{\left(\sum_j w_j^{(2)} q_l(y_j) q_v(y_j)\right)}_{=0 \text{ als } l \neq v}.
$$

De coëfficiënten zijn direct:

$$
\boxed{a_{k,l} = \frac{\sum_i \sum_j w_i^{(1)} w_j^{(2)}\, f_{i,j}\, p_k(x_i)\, q_l(y_j)}{\|p_k\|^2 \|q_l\|^2}.}
$$

De benadering zelf:

$$
z(x, y) = \sum_{k=0}^m \sum_{l=0}^n a_{k,l}\, p_k(x)\, q_l(y).
$$

---

### 6.3.3 Herhaalde eendimensionale benadering

> ⚠️ **Belangrijk voor examen:** Dit is de efficiënte implementatiestrategie.

De formule voor $a_{k,l}$ kan worden gefactoriseerd:

$$
a_{k,l} = \frac{\sum_{j=1}^N w_j^{(2)}\, b_{k,j}\, q_l(y_j)}{\|q_l\|^2}, \qquad \text{met} \quad b_{k,j} = \frac{\sum_{i=1}^M w_i^{(1)}\, f_{i,j}\, p_k(x_i)}{\|p_k\|^2}.
$$

**Algoritme in twee stappen:**

**Stap 1 — Benadering langs de $x$-richting:**

Voor elke $j = 1, \ldots, N$: bereken de kleinste-kwadratenbenadering van de $M$ waarden $f_{1,j}, f_{2,j}, \ldots, f_{M,j}$ (een doorsnede met het vlak $y = y_j$) door een veelterm van graad $m$ in $x$. Dit levert coëfficiënten $b_{k,j}$, voor $k = 0, \ldots, m$.

$$
g_j(x) = \sum_{k=0}^m b_{k,j}\, p_k(x).
$$

**Stap 2 — Benadering langs de $y$-richting:**

Voor elke $k = 0, \ldots, m$: beschouw de $N$ waarden $b_{k,1}, b_{k,2}, \ldots, b_{k,N}$ als een discrete functie in $y$. Benader die door een veelterm van graad $n$ in $y$:

$$
h_k(y) = \sum_{l=0}^n a_{k,l}\, q_l(y).
$$

**Uitvoering met de methode van Forsythe:**
- Stap 1: $N$ keer Forsythe in de $x$-richting ($M$ punten, graad $m$).
- Stap 2: $(m+1)$ keer Forsythe in de $y$-richting ($N$ punten, graad $n$).

**Voordeel:** Geen nieuw algoritme nodig. Bestaande eendimensionale implementatie tweemaal toepassen.

**Evaluatie** van $z(x,y)$ voor een punt $(x,y)$: ook via herhaalde eendimensionale evaluatie — eerst Smith-schema langs $x$, dan langs $y$.

**Uitbreiding naar hogere dimensies:** de aanpak veralgemeent direct naar 3D en meer.

---

### 6.3.4 Tensorproduct-splinefuncties

#### Definitie

Beschouw een rechthoekig gebied $R = [a,b] \times [c,d]$ met knooppuntenrijen $\{\tau_i\}_{i=0}^g$ in de $x$-richting en $\{\mu_j\}_{j=0}^h$ in de $y$-richting.

### Definitie 6.3.1

Een **tensorproduct-splinefunctie** van graad $k$ in $x$ en graad $l$ in $y$ is een functie $s(x,y)$ op $R$ zodat:

1. Op elke deelrechthoek $[\tau_i, \tau_{i+1}] \times [\mu_j, \mu_{j+1}]$ is $s(x,y)$ een veelterm van graad $k$ in $x$ en graad $l$ in $y$.
2. Alle partiële afgeleiden $\frac{\partial^{p+q} s}{\partial x^p \partial y^q}$ voor $0 \leq p < k$ en $0 \leq q < l$ zijn continu op $R$.

**Dimensie van de splineruimte:** $(g+k) \times (h+l)$.

**Speciale naam:** Graad $k = l = 3$ geeft **bikubische splines**.

#### B-splinebasis

Na uitbreiding van de knooppuntenrijen (elk $2k$ resp. $2l$ extra knooppunten), vormen de **kruisproducten** van eendimensionale B-splines een basis:

$$
\boxed{s(x, y) = \sum_{i=-k}^{g-1} \sum_{j=-l}^{h-1} c_{i,j}\, N_{i,k+1}^\tau(x)\, N_{j,l+1}^\mu(y).}
$$

De superscripten $\tau$ en $\mu$ duiden aan welke knooppuntenrij gebruikt wordt.

**Totaal aantal coëfficiënten:** $(g+k) \times (h+l)$ — gelijk aan de dimensie van de splineruimte.

#### Herhaalde eendimensionale evaluatie

Om $s(x,y)$ te evalueren in $(x,y) \in [\tau_r, \tau_{r+1}) \times [\mu_s, \mu_{s+1})$:

**Stap 1 — Koëfficiënten als functies van $y$:**

$$
a_i = \sum_{j=s-l}^{s} c_{i,j}\, N_{j,l+1}^\mu(y), \qquad i = r-k, \ldots, r.
$$

Dit zijn $k+1$ eendimensionale spline-evaluaties in $y$ (de Boor-algoritme).

**Stap 2 — Evalueer de eendimensionale spline in $x$:**

$$
s(x,y) = \sum_{i=r-k}^r a_i\, N_{i,k+1}^\tau(x).
$$

Dit is één de Boor-evaluatie in $x$.

**Totaal:** $k+1 + 1 = k+2$ toepassingen van het de Boor-algoritme. De tweedimensionale evaluatie reduceert volledig tot herhaalde eendimensionale evaluatie.

---

## Samenvatting

### Discrete veeltermbenadering

| Aanpak | Matrix | Conditionering | Graadaanpassing |
|---|---|---|---|
| Monomiale basis | vol, Hilbert-achtig | zeer slecht | volledig herbouw |
| Continue orthogonale basis | bijna diagonaal | goed | gedeeltelijk herbouw |
| **Forsythe (discrete orth. basis)** | **diagonaal** | **optimaal** | **enkel nieuwe term** |

### Discrete splinebenadering

| Aanpak | Systeem | Structuur | Oplossingsmethode |
|---|---|---|---|
| Normaalvergelijkingen | $(n+k)^2$ | bandmatrix breedte $2k+1$ | Cholesky band |
| QR van overgedetermineerd | $N \times (n+k)$ | band $k+1$ per rij | Givens-rotaties |

**Schoenberg-Whitney:** de normaalmatrix is SPD als en slechts als in de drager van elke B-spline minstens één meetpunt ligt.

### Tweedimensionale benadering

**Sleuteltruc — Tensorproduct + separabel gewicht → herhaalde 1D-benadering:**

$$
z(x,y) = \sum_{k,l} a_{k,l}\, p_k(x)\, q_l(y), \quad a_{k,l} \text{ via 2 × Forsythe}
$$

```
Stap 1: N × (M punten, graad m in x)  → b[k,j]
Stap 2: (m+1) × (N punten, graad n in y) → a[k,l]
```

Dezelfde structuur geldt voor tensorproduct-splines: herhaalde eendimensionale de Boor-evaluatie.
