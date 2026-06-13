# Ontbrekende Examenvragen — Uitgewerkte Antwoorden

> Gebaseerd op: *Examenvragen tot 2025 (Maxim Bijnens)* en *proofs (1).pdf*

---

## Inhoudstafel

### DEEL 1 — Ontbrekende Examenvragen

- **[H2: Benadering van vectoren](#hoofdstuk-2-benadering-van-vectoren1)**
- **[H3: Benaderen in een unitaire ruimte](#hoofdstuk-3-benaderen-in-een-unitaire-ruimte)**
- **[H5: Benaderen met splinefuncties](#hoofdstuk-5-benaderen-met-splinefuncties)**
- **[H8: Grafen en eigenwaarden](#hoofdstuk-8-grafen-en-eigenwaarden)**
- **[H9: Eigenwaardenalgoritmes](#hoofdstuk-9-eigenwaardenalgoritmes)**
- **[H12: Ijle representaties en benaderingen](#hoofdstuk-12-ijle-representaties-en-benaderingen)**

### DEEL 2 — Ontbrekende Bewijzen

- **[H2: Benadering van vectoren](#hoofdstuk-2-benadering-van-vectoren)**
- **[H3: Benaderen in een unitaire ruimte — Bewijzen](#hoofdstuk-3-benaderen-in-een-unitaire-ruimte--bewijzen)**
- **[H4: Veeltermbenaderingen](#hoofdstuk-4-veeltermbenaderingen)**
- **[H5: Splinefuncties — Bewijzen](#hoofdstuk-5-benaderen-met-splinefuncties--bewijzen)**
- **[H8: Grafen en eigenwaarden — Bewijzen](#hoofdstuk-8-grafen-en-eigenwaarden--bewijzen)**
- **[H9: Eigenwaardenalgoritmes — Bewijzen](#hoofdstuk-9-eigenwaardenalgoritmes--bewijzen)**
- **[H11: Optimalisatie-algoritmes](#hoofdstuk-11-optimalisatie-algoritmes)**
- **[H12: Ijle representaties — Bewijzen](#hoofdstuk-12-ijle-representaties--bewijzen)**

---

## DEEL 1 — Ontbrekende Examenvragen

---

## Hoofdstuk 2: Benadering van vectoren1
---

## V-A — QR-factorisatie: vergelijking van drie methodes

**Vraag (H2):** Bespreek de QR-factorisatie en hoe je die kan opstellen.
- (a) Bespreek de geneste structuur van de QR-factorisatie en waarom deze structuur kan gebruikt worden in de Arnoldi-iteratie.
- (d) Bespreek kort hoe men met Givensrotaties een QR-factorisatie kan opstellen.
- (f) Vergelijk de drie methodes op basis van stabiliteit, orthogonaliteit en rekenkost. Bij welk type matrices verkies je elke methode?

---

### (a) Geneste structuur van de QR-factorisatie

De QR-factorisatie $A = QR$ (of $A = \hat{Q}\hat{R}$ in de onvolledige versie) heeft een **geneste structuur**:

$$\langle \mathbf{q}_1, \ldots, \mathbf{q}_k \rangle = \langle \mathbf{a}_1, \ldots, \mathbf{a}_k \rangle, \quad k = 1, \ldots, n.$$

Dit betekent dat na $k$ stappen van het Gram-Schmidt algoritme de ruimte opgespannen door de eerste $k$ orthonormale vectoren gelijk is aan de ruimte opgespannen door de eerste $k$ kolomvectoren van $A$.

**Verband met Arnoldi:** De Arnoldi-iteratie bouwt stap voor stap een orthonormale basis van de Krylov-ruimte $\mathcal{K}_k(A,\mathbf{b}) = \langle \mathbf{b}, A\mathbf{b}, \ldots, A^{k-1}\mathbf{b} \rangle$ op. Precies zoals bij QR groeit de basis stap voor stap: na iteratie $k$ heeft men een orthonormale basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_k\}$ voor $\mathcal{K}_k$. De geneste structuur is essentieel: de basis voor $\mathcal{K}_{k+1}$ bevat de volledige basis voor $\mathcal{K}_k$. Dit laat toe om de Arnoldi-iteratie on-the-fly uit te breiden zonder vroegere berekeningen te herhalen.

> ⚠️ **Examentip:** De Krylov-ruimte heeft de geneste structuur $\mathcal{K}_1 \subset \mathcal{K}_2 \subset \cdots$, precies zoals de geprojecteerde deelruimten bij QR. Arnoldi is dus een QR-factorisatie van de matrix $[\mathbf{b},\, A\mathbf{b},\, A^2\mathbf{b},\, \ldots]$.

---

### (d) Givensrotaties

Een **Givensrotatie** $G(i,j,\theta)$ is een unitaire matrix die in een 2×2 deelmatrix een rotatie uitvoert:
$$G(i,j,\theta) = \begin{pmatrix} \ddots & & & \\ & c & \cdots & s \\ & \vdots & \ddots & \vdots \\ & -\bar{s} & \cdots & c \\ & & & & \ddots \end{pmatrix}, \quad c = \cos\theta,\; s = \sin\theta.$$

**Toepassing op QR:** Voor elk element $A_{ij}$ onder de diagonaal kiest men $\theta$ zodat $G^* A$ dat element op nul brengt:
$$G(k, k+1, \theta)^* \begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} r \\ 0 \end{pmatrix}, \quad r = \sqrt{|a|^2 + |b|^2}.$$

Door alle subelementen één voor één te elimineren (van links naar rechts, van boven naar onder) bekomt men $R = G_N^* \cdots G_1^* A$ en $Q = G_1 \cdots G_N$.

**Voordeel:** Elke rotatie raakt slechts 2 rijen aan → ideaal voor **spaarse matrices** of wanneer slechts enkele elementen worden toegevoegd/gewijzigd (bijv. toevoegen van rijen in een bestaande QR-factorisatie).

---

### (f) Vergelijking van de drie methodes

| Eigenschap | Gram-Schmidt | Givensrotatie | Householder |
|---|---|---|---|
| **Rekenkost** (enkel R) | $\sim 2mn^2$ | $3mn^2 - n^3$ | $2mn^2 - \tfrac{2}{3}n^3$ |
| **Rekenkost** (R én Q) | $\sim 2mn^2$ | $6m^2n - n^3$ | $4m^2n - \tfrac{2}{3}n^3$ |
| **Orthogonaliteit** | Gevaarlijke aftrekkingen → fouten | Unitair → machineprecisie | Unitair → machineprecisie |
| **Stabiliteit** | Slechter (verbetering met MGS + herortho) | Zeer goed | Zeer goed |
| **Type matrices** | Kleine/middelgrote, goed gecond.; onvolledige QR voldoende; Q expliciet nodig | Spaarse matrices; rijen/kolommen toevoegen in bestaande QR | Grote, dense, slecht gecond.; precisie belangrijk; Q niet expliciet nodig |

> ⚠️ **Examentip:** Givensrotaties zijn inefficiënt voor grote dense matrices (elk element apart elimineren). Householder is de standaardkeuze voor stabiliteit bij grote problemen. Gram-Schmidt kost $O(mn^2)$ voor de volledige $\hat{Q}$ — voordelig wanneer $m \gg n$ en enkel $\hat{Q}$ nodig is.

---

## Hoofdstuk 3: Benaderen in een unitaire ruimte

---

## V-B — Afstandsfunctie en strikte norm

**Vraag (H3):** Beschouw de functie $f(x,y) = \dfrac{|x-y|}{1+|x-y|}$; $x, y \in \mathbb{R}$.
- (a) Is $f$ een afstand?
- (b) Kunnen we de geïnduceerde norm definiëren voor $f$?
- (c) Geef een voorbeeld van een strikte norm. Waarom is deze norm strikt?
- (d) Bewijs dat de beste benadering in een strikt genormeerde ruimte bestaat en uniek is.

---

### (a) Is $f$ een afstand?

We verifiëren de vier voorwaarden voor een afstandsfunctie $\rho(x,y) = f(x,y)$:

**1. Niet-negativiteit:** $|x-y| \geq 0$, dus $f(x,y) \geq 0$. ✓

**2. Nulvector:** $f(x,y) = 0 \Leftrightarrow |x-y| = 0 \Leftrightarrow x = y$. ✓

**3. Symmetrie:** $f(x,y) = \dfrac{|x-y|}{1+|x-y|} = \dfrac{|y-x|}{1+|y-x|} = f(y,x)$. ✓

**4. Driehoeksongelijkheid:** Te bewijzen: $f(x,z) \leq f(x,y) + f(y,z)$.

Stel $a = |x-y|$ en $b = |y-z|$. De gewone driehoeksongelijkheid geeft $|x-z| \leq a + b$. De functie $t \mapsto \frac{t}{1+t}$ is strikt stijgend, dus:
$$f(x,z) = \frac{|x-z|}{1+|x-z|} \leq \frac{a+b}{1+a+b} = \frac{a}{1+a+b} + \frac{b}{1+a+b} \leq \frac{a}{1+a} + \frac{b}{1+b} = f(x,y) + f(y,z). \checkmark$$

**Conclusie:** $f$ is een afstand. $\square$

---

### (b) Geïnduceerde norm definiëren?

**Antwoord: Nee.**

Een norm moet **homogeen** zijn: $\|ax\| = |a| \cdot \|x\|$. Controleer $f(\alpha x, 0)$ voor $y = 0$:
$$f(\alpha x, 0) = \frac{|\alpha x|}{1 + |\alpha x|} = \frac{|\alpha||x|}{1 + |\alpha||x|}.$$

Dit is **niet** gelijk aan $|\alpha| \cdot f(x, 0) = \frac{|\alpha||x|}{1 + |x|}$ (tenzij $|\alpha| = 1$).

$f$ is dus niet homogeen en definieert geen norm. Bijgevolg kan er ook geen geïnduceerde norm bestaan. ($f$ is wel supermetrisch, maar beide eigenschappen — homogeniteit én supermetriciteit — zijn vereist.)

---

### (c) Voorbeeld van een strikte norm

**Definitie:** Een norm $\|\cdot\|$ is **strikt** als voor alle $x_1 \neq x_2$ met $\|x_1\| = \|x_2\| = 1$ geldt: $\|x_1 + x_2\| < 2$.

**Voorbeeld: de Euclidische norm** $\|x\|_2 = \sqrt{\sum_i x_i^2}$ op $\mathbb{R}^n$.

**Waarom strikt?** Uit de parallellogramgelijkheid:
$$\|x_1 + x_2\|_2^2 = -\|x_1 - x_2\|_2^2 + 2(\|x_1\|_2^2 + \|x_2\|_2^2) = 4 - \|x_1 - x_2\|_2^2.$$

Als $x_1 \neq x_2$ en $\|x_1\|_2 = \|x_2\|_2 = 1$, dan $\|x_1 - x_2\|_2 > 0$, dus $\|x_1 + x_2\|_2 < 2$. $\square$

**Intuïtie:** De eenheidsbol van de Euclidische norm is een gladde bol (strikt convex). De maximumnorm $\|\cdot\|_\infty$ is **niet** strikt: $f(x)=1$, $g(x)=x$ op $[0,1]$ geeft $\|f+g\|_\infty = 2$.

---

### (d) Bewijs: beste benadering bestaat en is uniek in strikt genormeerde ruimte

**Stelling:** Zij $\mathcal{D}$ een eindigdimensionale deelruimte van een strikt genormeerde ruimte $V$ en $v \in V$. Dan bestaat de beste benadering van $v$ in $\mathcal{D}$ en is deze uniek.

#### Existentie

Definieer $d = \inf\{\|v - w\| : w \in \mathcal{D}\}$. Volgens de definitie van infimum bestaat er een rij $\{w_k\}$ in $\mathcal{D}$ zodat $\|v - w_k\| \to d$.

De rij is uniform begrensd: $\|w_k\| \leq \|w_k - v\| + \|v\| \leq \|w_1 - v\| + \|v\|$.

Stel $n = \dim(\mathcal{D})$ en beschouw een basis $\{a_1, \ldots, a_n\}$. Dan $w_k = \sum c_{ki} a_i$ en de rij $\{(c_{k1},\ldots,c_{kn})\}$ is begrensd. Volgens **Weierstrass-Bolzano** bestaat een convergente deelrij met limiet $(\hat{c}_1, \ldots, \hat{c}_n)$. Definieer $\hat{w} = \sum \hat{c}_i a_i \in \mathcal{D}$. Dan:
$$\|v - \hat{w}\| \leq \underbrace{\|v - w_k\|}_{\to d} + \underbrace{\|w_k - \hat{w}\|}_{\to 0} \to d.$$
Dus $\|v - \hat{w}\| = d$: de beste benadering bestaat.

#### Uniciteit

Bewijs uit het ongerijmde. Veronderstel twee beste benaderingen $v_1 \neq v_2$ met $\|v - v_1\| = \|v - v_2\| = d$.

Definieer $e_i = \frac{1}{d}(v - v_i)$ op de eenheidsbol, $e_1 \neq e_2$. Omdat de norm **strikt** is:
$$\left\|v - \frac{v_1 + v_2}{2}\right\| = d \left\|\frac{1}{2}e_1 + \frac{1}{2}e_2\right\| < d.$$

Maar $\frac{v_1+v_2}{2} \in \mathcal{D}$ is een **betere** benadering dan $v_1$ — tegenspraak. $\square$

> ⚠️ **Examentip:** Existentie via Weierstrass-Bolzano; uniciteit via het middelpunt in een strikte norm.

---

## Hoofdstuk 5: Benaderen met splinefuncties

---

## V-C — B-splines (volledige vraag)

**Vraag (H5):** Bespreek B-splines.

- (b) Hoe worden B-splines gecombineerd tot de algemene splinefunctie $s(x)$?
- (c) Leg uit en bewijs hoe we de splinefunctie efficiënt kunnen evalueren.
- (d) Wat verandert er bij samenvallende knooppunten? Wanneer gebruiken?

---

### (b) Combineren tot algemene splinefunctie

Elke splinefunctie $s(x)$ van graad $k$ op knooppunten $t_0 < \cdots < t_n$ kan **uniek** geschreven worden als:
$$\boxed{s(x) = \sum_{i=-k}^{n-1} c_i \, N_{i,k+1}(x),}$$
waarbij $\{c_i\}$ de **controlepunten** (B-splinecoëfficiënten) zijn.

**Dimensietelling:** Er zijn $n + k$ basisfuncties (klopt met dim. splineruimte $= n + k$).

**Eigenschappen:**
- **Lokale controle:** De waarde van $s$ in $[t_j, t_{j+1})$ hangt enkel af van $c_{j-k}, \ldots, c_j$ (slechts $k+1$ coëfficiënten).
- **Convex-hul-eigenschap:** $s(x)$ ligt in de convexe hul van de naburige controlepunten.

---

### (c) Efficiënte evaluatie: het algoritme van de Boor

**Probleem:** Gegeven $x \in [t_j, t_{j+1})$, bereken $s(x)$ efficiënt.

**Algoritme van de Boor:** Via herhaalde lineaire interpolatie.

**Bewijs/afleiding:**

Start van $s(x) = \sum_{i=-k}^{n-1} c_i N_{i,k+1}(x)$. Pas de recursiebetrekking toe:
$$s(x) = \sum_{i=-(k-1)}^{n-1} c_i^{[1]} N_{i,k}(x),$$
waarbij $c_i^{[1]} = \alpha_{i,1} c_i + (1 - \alpha_{i,1}) c_{i-1}$ met $\alpha_{i,1} = \dfrac{x - t_i}{t_{i+k} - t_i}$.

Herhaal dit $k$ keer:
$$s(x) = \sum_{i=0}^{n-1} c_i^{[k]} N_{i,1}(x).$$

Voor $x \in [t_j, t_{j+1})$ is enkel $N_{j,1}(x) = 1$, alle andere zijn nul. Dus:
$$\boxed{s(x) = c_j^{[k]},}$$
met de algemene recursie:
$$c_i^{[r]} = \alpha_{i,r} \, c_i^{[r-1]} + (1 - \alpha_{i,r}) \, c_{i-1}^{[r-1]}, \quad \alpha_{i,r} = \frac{x - t_i}{t_{i+k+1-r} - t_i}. \quad \square$$

**Schema voor $k=3$, $x \in [t_j, t_{j+1})$:**
```
c_{j-3}^[0]
             → c_{j-2}^[1]
c_{j-2}^[0]              → c_{j-1}^[2]
             → c_{j-1}^[1]              → c_j^[3] = s(x)
c_{j-1}^[0]              → c_j^[2]
             → c_j^[1]
c_j^[0]
```

**Rekenkost:** $O(k^2)$ bewerkingen, slechts $k+1$ coëfficiënten nodig.

> ⚠️ **Examentip:** De Boor is de B-spline analogon van Horner. Telkens interpoleren tussen twee naburige coëfficiënten met gewicht $\alpha_{i,r}$.

---

### (d) Samenvallende knooppunten

**Definitie:** Als een knooppunt $t_i$ met multipliciteit $\mu > 1$ voorkomt ($t_i = t_{i+1} = \cdots = t_{i+\mu-1}$), spreekt men van **samenvallende knooppunten**.

**Effect op continuïteit:** Een knooppunt met multipliciteit $\mu$ verlaagt de continuïteit van $C^{k-1}$ naar $C^{k-\mu}$. Bij multipliciteit $k+1$: discontinu.

**Wanneer gebruiken?**
- **Interpolatie-eis:** Dwingt de spline door een bepaald punt te gaan (multipliciteit $k$).
- **Randcondities:** Samenvallende randknooppunten zorgen dat de spline het eerste en laatste controlepunt interpoleert.
- **Scherpe knik:** Multipliciteit $k$ laat een scherpe hoek toe (bv. in CAD-ontwerp).

---

## Hoofdstuk 8: Grafen en eigenwaarden

---

## V-D — Graafpartitionering

**Vraag (H8):** Bespreek graafpartitionering.
- (a) Hoe komen we tot de formule van dit probleem?
- (b) Hoe 'relaxeren' we dit probleem?
- (c) Toon aan dat de oplossing de eigenvector is van de op een na kleinste eigenwaarde.
- (d) Welk algoritme/methode gebruiken we?

---

### (a) Formule van het partitioneringsprobleem

**Probleemstelling:** Verdeel de $N$ knopen van een graaf in twee gelijke deelgrafen $K_1$ en $K_2$ door zo weinig mogelijk bogen door te knippen.

**Graaf-Laplaciaan:** $L = D - A$, waarbij $D$ de graadmatrix (diagonaal) en $A$ de verbindingsmatrix.

**Partitievector:** $\mathbf{x} \in \{-1, +1\}^N$ met $x_i = +1$ als $i \in K_1$ en $x_i = -1$ als $i \in K_2$.

**Gelijke deelgrafen:** $\mathbf{1}^T\mathbf{x} = 0$.

**Aantal doorgeknippen bogen:**
$$\frac{1}{4} \mathbf{x}^T L \mathbf{x} = \frac{1}{4} \sum_{i,j} L_{ij} x_i x_j = \frac{1}{2} \sum_{(i,j) \in B} (-x_i x_j + 1) = \#(B_\mathbf{x}).$$

**Optimaliseringsprobleem:**
$$\min_{\mathbf{x} \in \{-1,+1\}^N} \frac{1}{4} \mathbf{x}^T L \mathbf{x} \quad \text{o.v.v.} \quad \mathbf{1}^T\mathbf{x} = 0.$$

---

### (b) Relaxatie

De binaire beperking maakt het probleem NP-hard. **Relaxatie:**
$$\mathbf{x} \in \{-1,+1\}^N \;\longrightarrow\; \|\mathbf{x}\|_2^2 = N.$$

Het gerelaxeerde probleem:
$$\min_{\mathbf{x} \in \mathbb{R}^N} \frac{1}{4} \mathbf{x}^T L \mathbf{x} \quad \text{o.v.v.} \quad \|\mathbf{x}\|_2^2 = N,\; \mathbf{1}^T\mathbf{x} = 0.$$

---

### (c) Oplossing is eigenvector van $\lambda_{N-1}$ (Fiedler-vector)

**Bewijs:**

Omdat $L$ symmetrisch is: $L = V\,\text{diag}(\lambda_1,\ldots,\lambda_N)\,V^T$, met $V$ orthogonaal.

**Eigenwaarde 0:** $L\mathbf{1} = \mathbf{0}$, dus $\lambda_N = 0$ met $\mathbf{v}_N = \frac{1}{\sqrt{N}}\mathbf{1}$.

**Substitutie** $\mathbf{x} = V\mathbf{y}$: de eis $\mathbf{1}^T\mathbf{x} = 0$ wordt $\mathbf{v}_N^T V\mathbf{y} = y_N = 0$.

Het probleem transformeert:
$$\min \frac{1}{4} \sum_{i=1}^N \lambda_i y_i^2 \quad \text{o.v.v.} \quad \|\mathbf{y}\|_2^2 = N,\; y_N = 0.$$

Om $\sum \lambda_i y_i^2$ te minimaliseren onder $\sum y_i^2 = N$ en $y_N = 0$: geef al het gewicht aan de kleinste vrije eigenwaarde $\lambda_{N-1}$:
$$\mathbf{y} = [0 \cdots 0 \; \sqrt{N} \; 0]^T \implies \mathbf{x} = V\mathbf{y} = \sqrt{N}\,\mathbf{v}_{N-1}. \quad \square$$

**Praktische partitie:**
$$K_1 = \{i : (\mathbf{v}_{N-1})_i \geq 0\}, \quad K_2 = \{i : (\mathbf{v}_{N-1})_i < 0\}.$$

---

### (d) Algoritme: Lanczos

Omdat $L$ **symmetrisch en spaars** is en we slechts de tweede-kleinste eigenwaarde nodig hebben, gebruiken we het **Lanczos-algoritme** (variant van Arnoldi voor symmetrische matrices):

- Enkel matrix-vector producten $L\mathbf{v}$ vereist ($O(nnz)$ per stap).
- Na $k \ll N$ iteraties benadert de Fiedler-eigenwaarde van $T_k = Q_k^T L Q_k$.

---

## V-E — PageRank

**Vraag (H8):** Bespreek het PageRank probleem.
- (a) Componenten van de graaf + rank toekennen
- (b) Benaderend eigenwaardeprobleem
- (c) Omvormen tot volwaardig eigenwaardeprobleem
- (d) Algoritme

---

### (a) Componenten en rank

**Graaf:** Knopen = webpagina's, gerichte bogen = hyperlinks.

**Rank-definitie:** Een pagina is centraal als ze veel inkomende links heeft van zelf centrale pagina's. De **link-matrix** $A$:
$$A_{ij} = \frac{1}{N_i} \quad \text{als pagina } P_i \text{ een link heeft naar } P_j \text{ (met } N_i \text{ uitlinks van } P_i\text{)}, \quad 0 \text{ anders.}$$

$A$ is **rij-stochastisch**. De rank-vector $\boldsymbol{\pi}$ voldoet aan $A^T \boldsymbol{\pi} = \boldsymbol{\pi}$ (eigenvector met eigenwaarde 1).

---

### (b) Benaderend eigenwaardeprobleem

De rank-vector $\boldsymbol{\pi}$ is een **eigenvector van $A^T$ met eigenwaarde 1**.

**Problemen:**

1. **Dangling nodes:** Pagina's zonder uitlinks → nulrij → matrix niet stochastisch.
2. **Niet-samenhangende graaf:** Eigenwaarde 1 niet uniek.

---

### (c) Omvorming tot volwaardig eigenwaardeprobleem

**Stap 1 — Dangling nodes:** Vervang nulrijen door $\frac{1}{N}\mathbf{1}^T$:
$$\tilde{A}_{ij} = A_{ij} + \frac{d_i}{N}, \quad d_i = 1 \text{ als pagina } i \text{ geen uitlinks heeft, anders } 0.$$
$\tilde{A}$ is rij-stochastisch.

**Stap 2 — Teleportatie (Google-matrix):**
$$\hat{A} = (1-\alpha)\tilde{A} + \frac{\alpha}{N}\mathbf{1}\mathbf{1}^T, \quad \alpha \approx 0.1.$$

$\hat{A}$ is nu een **irreduceerbare rij-stochastische matrix**: eigenwaarde 1 is enkelvoudig en dominant (Perron-Frobenius). De aanpassingen hebben weinig effect op de ranking ($\alpha$ is klein).

---

### (d) Algoritme: methode van de machten

$$\boldsymbol{\pi}^{(k+1)} = \hat{A}^T \boldsymbol{\pi}^{(k)}, \quad \boldsymbol{\pi}^{(0)} = \frac{1}{N}\mathbf{1}.$$

**Convergentiefactor:** $|\lambda_2(\hat{A}^T)| \leq 1 - \alpha$ — typisch $50$–$100$ iteraties.

---

## Hoofdstuk 9: Eigenwaardenalgoritmes

---

## V-F — Methode van de machten en deelruimte-iteratie

**Vraag (H9):** Bespreek hoe de deelruimte-iteratie de methode van de machten veralgemeent.
- (a) Basisiteratie + convergentie + convergentiefactor
- (b) Dominante eigenwaarde schatten
- (c) Kleinste modulus / dichtst bij $\sigma$
- (d) Deelruimte-iteratie
- (e) Eigenwaarden en -vectoren verkrijgen

---

### (a) Methode van de machten

**Doel:** Bereken de dominante eigenwaarde $\lambda_1$ (grootste $|\lambda_i|$) van $A \in \mathbb{C}^{m \times m}$.

Stel de eigenwaarden geordend: $|\lambda_1| \geq |\lambda_2| \geq \cdots \geq |\lambda_m|$, met bijbehorende eigenvectoren $\mathbf{v}_1, \ldots, \mathbf{v}_m$.

**Iteratie:** start met willekeurige $\mathbf{x}^{(0)}$, herhaal:
$$\mathbf{p}^{(k)} = A\mathbf{x}^{(k-1)}, \qquad \mathbf{x}^{(k)} = \frac{\mathbf{p}^{(k)}}{\|\mathbf{p}^{(k)}\|_2}.$$

**Convergentie:** Schrijf $\mathbf{x}^{(0)} = \sum_i c_i \mathbf{v}_i$ (eigenbasis, $c_1 \neq 0$). Dan:
$$A^k \mathbf{x}^{(0)} = \lambda_1^k \left(c_1 \mathbf{v}_1 + \sum_{i \geq 2} c_i \left(\frac{\lambda_i}{\lambda_1}\right)^k \mathbf{v}_i\right).$$

Als $|\lambda_1| > |\lambda_2|$: elke term $(\lambda_i/\lambda_1)^k \to 0$, dus $\mathbf{x}^{(k)} \to \mathbf{v}_1$.

**Convergentiefactor:** $\rho = \left|\dfrac{\lambda_2}{\lambda_1}\right|$ — lineaire convergentie.

---

### (b) Dominante eigenwaarde schatten

Via het **Rayleigh-quotiënt**:
$$r(\mathbf{x}^{(k)}) = \frac{(\mathbf{x}^{(k)})^* A\, \mathbf{x}^{(k)}}{(\mathbf{x}^{(k)})^* \mathbf{x}^{(k)}} \xrightarrow{k\to\infty} \lambda_1.$$

Convergentie is **kwadratisch** voor Hermitische matrices (want $\mathbf{x}^{(k)}$ convergeert lineair en het Rayleigh-quotiënt is een tweedeordesbenadering).

---

### (c) Kleinste modulus / dichtst bij $\sigma$

**Kleinste eigenwaarde — inverse iteratie:** Pas de methode van de machten toe op $A^{-1}$: los in elke stap $A\mathbf{p}^{(k)} = \mathbf{x}^{(k-1)}$ op (via vooraf berekende LU). De dominante eigenwaarde van $A^{-1}$ is $1/\lambda_{\min}$.

**Eigenwaarde dichtst bij $\sigma$ — verschuiving:** Pas inverse iteratie toe op $(A - \sigma I)^{-1}$:
$$\rho = \left|\frac{\lambda_j - \sigma}{\lambda_{\text{tweede dichtste}} - \sigma}\right|.$$

Hoe dichter $\sigma$ bij $\lambda_j$, hoe sneller de convergentie.

---

### (d) Deelruimte-iteratie

**Doel:** Bereken de $p$ dominante eigenvectoren van $A \in \mathbb{C}^{m \times m}$ tegelijk ($p \ll m$; veralgemening van de methode van de machten naar een deelruimte).

**Iteratie:** start met $Q_0 \in \mathbb{C}^{m \times p}$ (orthonormale kolommen):
$$P_k = A Q_{k-1}, \qquad Q_k R_k = P_k \quad \text{(QR-factorisatie)}, \quad k = 1, 2, \ldots$$

De QR-stap vervangt de normalisatie uit de methode van de machten: ze herstelt zowel lengte als orthogonaliteit van alle $p$ vectoren tegelijk.

**Convergentie:** De kolomruimte van $Q_k$ convergeert naar $\langle \mathbf{v}_1, \ldots, \mathbf{v}_p \rangle$ (de deelruimte van de $p$ dominante eigenvectoren).

**Convergentiefactor:** $\left|\dfrac{\lambda_{p+1}}{\lambda_p}\right|$ per iteratie — lineaire convergentie.

---

### (e) Eigenwaarden en -vectoren verkrijgen

Na convergentie ($Q_k \approx [\mathbf{v}_1 \cdots \mathbf{v}_p]$): **Rayleigh-Ritz projectie**:
$$Q_k^* A Q_k \in \mathbb{C}^{p \times p}.$$

Bereken de eigenontbinding van dit kleine stelsel (via het volledige QR-algoritme). De benaderingen zijn:

- **Eigenwaarden van $A$:** eigenwaarden van $Q_k^* A Q_k$.
- **Eigenvectoren van $A$:** $Q_k \mathbf{y}_i$, met $\mathbf{y}_i$ de eigenvectoren van $Q_k^* A Q_k$.

---

## V-G — Impliciet verschoven QR-algoritme

**Vraag (H9):** Bespreek het impliciet verschoven QR-algoritme.
- (a) Deelruimte-iteratie op $p(A)$ + interpretatie $\rho_i$
- (b) Waarom inefficiënt?
- (c) Werking voor Hessenbergmatrix
- (d) Bewijs equivalentie met deelruimte-iteratie op $p(A)$
- (e) Keuze van $\rho_i$

---

### (a) Deelruimte-iteratie op $p(A)$

**Idee:** Pas de deelruimte-iteratie toe op:
$$p(A) = (A - \rho_1 I)(A - \rho_2 I) \cdots (A - \rho_\ell I).$$

**Interpretatie van $\rho_i$:** De $\rho_i$ zijn **verschuivingen** (shifts). Als $\rho_i \approx \lambda_j$, dan $|p(\lambda_j)| \approx 0$ — die eigenvector wordt gesupprimeerd. De convergentiefactor wordt:
$$\left|\frac{p(\lambda_{n+1})}{p(\lambda_n)}\right| \ll \left|\frac{\lambda_{n+1}}{\lambda_n}\right|.$$

Goede keuze van $\rho_i$ geeft **kubische** of betere convergentie.

---

### (b) Inefficiëntie

Het expliciet berekenen van $p(A) = (A-\rho_1 I)\cdots(A-\rho_\ell I)$ kost $O(\ell m^2)$ per iteratiestap. Bovendien gaat de spaarse structuur van $A$ verloren. Voor grote matrices volledig onhaalbaar.

---

### (c) Werking voor Hessenbergmatrix

**Preprocessing:** Breng $A^{(0)}$ naar **Hessenbergvorm** via Householder ($O(m^3)$, éénmalig).

**Eén iteratiestap** (met $\ell$ verschuivingen):

$$p(A^{(k-1)}) = Q^{(k)} R^{(k)} \implies A^{(k)} = (Q^{(k)})^* A^{(k-1)} Q^{(k)}.$$

**Eigenschappen van $A^{(k)}$:**
1. **Gelijkvormig** met $A^{(k-1)}$ → zelfde eigenwaarden.
2. **Opnieuw in Hessenbergvorm** (volgt uit structuur van $Q^{(k)}$).
3. **Eerste kolom** van $Q^{(k)}$ is $\propto p(A^{(k-1)}) e_1$.

**Impliciet:** $p(A^{(k-1)})$ wordt nooit expliciet berekend; $Q^{(k)}$ wordt opgebouwd via Givensrotaties.

---

### (d) Bewijs equivalentie met deelruimte-iteratie op $p(A)$

**Schets (Bewijs 25, proofs.pdf):**

Omdat $A^{(k)}$ een Hessenbergmatrix is, kan $A^{(k-1)} Q^{(k)} = Q^{(k)} A^{(k)}$ geïnterpreteerd worden als de recursiebetrekking van een **Arnoldi-iteratie** op $A^{(k-1)}$. Daardoor:
$$\langle q_1^{(k)}, \ldots, q_j^{(k)} \rangle = \mathcal{K}_j(A^{(k-1)}, q_1^{(k)}).$$

Omdat $q_1^{(k)} \propto p^{(k)}(A^{(k-1)}) e_1$:
$$\mathcal{K}_j(A^{(k-1)}, q_1^{(k)}) = p^{(k)}(A^{(k-1)}) \mathcal{K}_j(A^{(k-1)}, e_1) = p^{(k)}(A^{(k-1)}) \langle e_1, \ldots, e_j \rangle,$$
waarbij gebruik gemaakt werd van de Hessenbergstructuur. De kolomruimten van $Q^{(k)}$ komen exact overeen met de deelruimten van een deelruimte-iteratie op $p^{(k)}(A^{(k-1)})$. $\square$

---

### (e) Keuze van de verschuivingen $\rho_i$

**Symmetrische matrices — Wilkinson-verschuiving:**

Beschouw het $2 \times 2$ rechterbenedenblok van de actuele Hessenbergmatrix. Kies $\rho$ als de **eigenwaarde die het dichtst bij $a_m$ ligt**. Geeft kubische convergentie.

**Willekeurige matrices — Francis dubbele verschuiving:**

Gebruik twee complexe geconjugeerde verschuivingen $\rho_1 = \bar{\rho}_2$ (eigenwaarden van het $2 \times 2$ rechterbenedenblok). Laat reële rekenkunde toe voor complexe eigenwaarden. Geeft ook kubische convergentie.

---

## Hoofdstuk 12: Ijle representaties en benaderingen

---

## V-H — Spectrale norm, Frobeniusnorm en PCA

**Vraag (H12):**
- Hoe bereken je $\|A\|_2$? Hoe bereken je $\|A\|_F$?
- Hoe wordt PCA toegepast voor datacompressie?

---

### Berekening van $\|A\|_2$

**Definitie:** $\|A\|_2 = \sigma_1(A)$ (de grootste singuliere waarde).

**Aanpak:**
1. **Bi-diagonalisatie:** Breng $A$ via unitaire transformaties naar bi-diagonale vorm $A_0$ (linker- en rechtermatrices onafhankelijk te kiezen, stabiel).
2. **Structuurmatrix:** Bereken eigenwaarden van:
   $$H = \begin{pmatrix} 0 & A_0^* \\ A_0 & 0 \end{pmatrix}.$$
   De eigenwaarden van $H$ zijn $\pm\sigma_i(A)$. Dit is stabieler dan eigenwaardenontbinding van $A^*A$ (die $\kappa^2$ geeft).
3. **Methode van de machten** op $H$ met $n=2$ (deelruimte-iteratie, want $H$ heeft zowel $+\sigma_1$ als $-\sigma_1$).

---

### Berekening van $\|A\|_F$

$$\|A\|_F = \sqrt{\sum_{i,j} |A_{ij}|^2} = \sqrt{\sum_i \sigma_i^2(A)}.$$

**Alle singuliere waarden nodig** → gebruik het **impliciet verschoven QR-algoritme** op de bi-diagonale matrix $A_0$. Rekenkost $O(mn^2)$ totaal.

---

### PCA voor datacompressie

**Setting:** $X \in \mathbb{R}^{m \times n}$ — $m$ observaties, elk een vector in $\mathbb{R}^n$ ($n$ variabelen). Typisch $m \gg n$ of omgekeerd, maar altijd $\mathrm{rang}(X) = r \leq \min(m,n)$.

**Stap 0 — Centreren:** Trek het gemiddelde per kolom af zodat elke kolom gemiddelde 0 heeft. Dit is noodzakelijk: PCA zoekt richtingen van maximale *variantie*, niet van maximale absolute waarden.

**Stap 1 — SVD:**
$$X = U\Sigma V^T, \quad \Sigma = \mathrm{diag}(\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r > 0).$$

De rechtersinguliere vectoren $\mathbf{v}_1, \ldots, \mathbf{v}_r \in \mathbb{R}^n$ zijn de **principale componenten** — richtingen in de variabelenruimte $\mathbb{R}^n$ geordend naar afnemende variantie. De variantie langs $\mathbf{v}_j$ is evenredig met $\sigma_j^2$.

**Stap 2 — Projectie op $k$ componenten:**

Definieer $V_k = [\mathbf{v}_1 \mid \cdots \mid \mathbf{v}_k] \in \mathbb{R}^{n \times k}$. Dan:
$$\tilde{X} = X V_k \in \mathbb{R}^{m \times k}.$$

Elke rij van $\tilde{X}$ is de lage-dimensionale coördinaatvector van de overeenkomstige observatie: $\tilde{X}_{ij}$ is de coördinaat van observatie $i$ langs principale component $\mathbf{v}_j$.

**Stap 3 — Reconstructie:**
$$\hat{X} = \tilde{X}V_k^T = XV_kV_k^T = U_k\Sigma_k V_k^T =: X_k^{\mathrm{svo}} \in \mathbb{R}^{m \times n}.$$

Dit is de **getrunceerde SVD van rang $k$** — dezelfde matrix als in de beste rang-$k$ benadering (Eckart–Young). PCA is dus de *optimale* lage-rang compressie in de Frobeniusnorm.

**Reconstructiefout:**
$$\|X - \hat{X}\|_F^2 = \sigma_{k+1}^2 + \cdots + \sigma_r^2.$$

**Keuze van $k$:** kies $k$ zodat de **verklaarde variantie** voldoende is:
$$\frac{\sum_{j=1}^k \sigma_j^2}{\sum_{j=1}^r \sigma_j^2} \geq \text{drempel (bv. 95\%)}.$$

**Motivatie:** Fysische systemen zijn inherent laag-dimensionaal; de data liggen nabij een laag-dimensionale deelruimte. De singuliere waarden $\sigma_j$ vallen snel af zodra dat het geval is.

---

## DEEL 2 — Ontbrekende Bewijzen

---

## Prioriteitstabel bewijzen

> Gebaseerd op **wiki_examenvragen** (echte examens 2023–2025) + Modelvragen_Samaey.md en Modelvragen_WM.md.

| Prio | Bewijs | Naam | H | Bewijs gevraagd in |
|:---:|:---:|---|:---:|---|
| 1 | 30 | Eckart–Young–Mirsky (spectraalnorm) | 12 | **Juni 2024 + Juni 2025** — elk jaar opnieuw, hint gegeven in opgave |
| 2 | 4 | Normaalstelsel $\Leftrightarrow$ overgedetermineerd stelsel | 2 | **Juni 2023 + 2024 (beide sessies)** — "toon aan dat oplossingen overeenkomen" |
| 3 | 13 | Drietermsrecursiebetrekking | 4 | **Juni 2025** — bewijs m.b.v. orthogonaliteit op lagere-graadsveeltermen |
| 4 | 24 | Graafpartitionering: Fiedler-vector | 8 | **7 Juni 2024** — "toon aan dat oplossing de eigenvector is bij $\lambda_{N-1}$" |
| 5 | 27 | CG-orthogonaliteitseigenschappen | 11 | Modelvraag Samaey V6 ("bijna exact dezelfde als modelvraag") |
| 6 | 28 | CG minimaliseert $A$-norm over Krylov | 11 | Modelvraag Samaey V6(c) — volledig bewijs |
| 7 | 11 | Orthogonale projector: $P = P^*$ | 3 | Modelvraag WM V1(c) ("exact dezelfde als modelvraag") |
| 8 | 8 | Norm unitair $\Leftrightarrow$ parallellogramgelijkheid | 3 | **Juni 2023** (is $I(f)$ een scalaire norm?) + Modelvraag WM V4(c) |
| 9 | 2 | Karakterisatie beste benadering | 2 | Impliciet in elk kleinste-kwadraten bewijs — Pythagoras, 3 regels |
| 10 | 29 | CG als veeltermbenadering | 11 | Modelvraag Samaey V6(d) — bijectie $\mathcal{K}_k \leftrightarrow \mathcal{N}_k$ |
| 11 | 14 | $\phi_k$ heeft $k$ enkelvoudige reële nulpunten | 4 | **Juni 2025** — "bewijs is niet nodig" maar structuur wordt gevraagd; WM V6(a) wél bewijs |
| 12 | 3 | Beste benadering bestaat en is uniek via $P_\mathcal{D}$ | 2 | Modelvraag WM V1 — uniciteit via Pythagoras |
| 13 | 9 | Eenheidsbol unitaire ruimte is strikt convex | 3 | **Juni 2023** (is $I(f)$ strikt?) + Modelvraag WM V4(b) |
| 14 | 26 | Lokaal minimum via Hessiaan | 11 | Modelvraag Samaey V7 — maar Armijo zelf is elk jaar gevraagd |
| 15 | 16 | Equioscillatiestelling | 4 | Modelvraag WM V6(c) context — zelfde structuur als Bewijs 14 |
| 16 | 10 | Grammatrix HPD $\Leftrightarrow$ lin. onafhankelijk | 3 | Modelvraag WM V5 — $v^*Gv = \|\cdot\|^2$ truc |
| 17 | 25 | Impliciet QR $\Leftrightarrow$ deelruimte-iteratie | 9 | Staat in V-G(d) — nog niet gezien in wiki, maar conceptueel zwaar |
| 18 | 7 | Unitaire ruimte $\Rightarrow$ genormeerde ruimte | 3 | Modelvraag WM V4 context — Cauchy-Schwarz |
| 19 | 6 | Beste benadering in strikt genormeerde ruimte | 3 | Staat al volledig in V-B(d) |
| 20 | 22 | De Boor-algoritme | 5 | Staat al volledig in V-C(c) |
| 21 | 5 | Norm $\Rightarrow$ metriek | 3 | Bijna triviaal, enkel driehoeksongelijkheid telt |
| 22 | 12 | Orthogonale veeltermen vormen basis | 4 | Triviaal: dim $= n+1$, lineaire onafhankelijkheid direct |
| 23 | 20 | Positiviteit: $M_{i,k+1} > 0$ | 5 | Inductie via recursie — zelden apart gevraagd |
| 24 | 30 | Eckart–Young–Mirsky (Frobeniusnorm) | 12 | Hulpstelling $\sigma_i(B-A) \geq \sigma_{k+i}(A)$ — technischer; soms gevraagd naast spectraalnorm |
| 25 | 17 | $M_{i,k+1}$ en $N_{i,k+1}$ zijn splinefuncties | 5 | Technisch, definitieel — nooit apart gezien |
| 26 | 19 | Lokale drager: $M_{i,k+1} = 0$ buiten $[t_i, t_{i+k+1}]$ | 5 | Technisch, bijna definitie |
| 27 | 23 | Afgeleide van splinefunctie | 5 | Technisch indexwerk, nooit apart gevraagd |
| 28 | 21 | Grensafgeleiden zijn nul | 5 | Volgt onmiddellijk uit continuïteit — laagste prioriteit |

> **Legende:** Prio 1–4 = op examen geweest in 2023–2025; 5–8 = staat in modelvragen als "exact dezelfde vraag"; 9–16 = sterk verwant; 17–28 = achtergrond of technisch detail.

---

## Hoofdstuk 2: Benadering van vectoren

---

## Bewijs 2 — Karakterisatie van de beste benadering

**Stelling:** $\hat{\mathbf{y}}$ is een beste benadering in deelruimte $\mathcal{D}$ voor $\mathbf{b}$ als en alleen als $\mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}$.

**Bewijs:** Neem willekeurige $\mathbf{y} \in \mathcal{D}$. Vermits $\mathbf{y} - \hat{\mathbf{y}} \in \mathcal{D}$ en $\mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}$, staan $(\mathbf{b}-\hat{\mathbf{y}})$ en $(\mathbf{y}-\hat{\mathbf{y}})$ loodrecht. Via **Pythagoras**:
$$\|\mathbf{b} - \mathbf{y}\|_2^2 = \|\mathbf{b} - \hat{\mathbf{y}}\|_2^2 + \|\hat{\mathbf{y}} - \mathbf{y}\|_2^2 \geq \|\mathbf{b} - \hat{\mathbf{y}}\|_2^2.$$
Dus $\mathbf{y}$ benadert $\mathbf{b}$ niet beter dan $\hat{\mathbf{y}}$. Omwille van het volgende resultaat kunnen we spreken van de beste benadering.

---

## Bewijs 3 — Beste benadering bestaat en is uniek via $P_\mathcal{D}$

**Stelling:** De beste benadering in $\mathcal{D}$ voor $\mathbf{b}$ is $\hat{\mathbf{y}} = P_\mathcal{D} \mathbf{b}$.

**Bewijs:** Ontbind $\mathbf{b} = P_\mathcal{D} \mathbf{b} + (I - P_\mathcal{D}) \mathbf{b}$. Dan:
$$\mathbf{b} - \hat{\mathbf{y}} = (I - P_\mathcal{D}) \mathbf{b} \in \mathcal{N}(P_\mathcal{D}) = \mathcal{D}^\perp \implies \mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}.$$
Uit Bewijs 2: $\hat{\mathbf{y}}$ is een beste benadering.

**Uniciteit:** Als ook $\mathbf{y} \in \mathcal{D}$ een beste benadering is:
$$\|\mathbf{b} - \mathbf{y}\|_2^2 = \|\mathbf{b} - \hat{\mathbf{y}}\|_2^2 + \|\mathbf{y} - \hat{\mathbf{y}}\|_2^2 \geq \|\mathbf{b} - \hat{\mathbf{y}}\|_2^2,$$
met gelijkheid enkel als $\|\mathbf{y} - \hat{\mathbf{y}}\|_2 = 0$, dus $\mathbf{y} = \hat{\mathbf{y}}$. $\square$

---

## Bewijs 4 — Normaalstelsel

**Stelling:** $\hat{\mathbf{x}}$ is oplossing van $\min\|A\mathbf{x}-\mathbf{b}\|_2$ $\Leftrightarrow$ $A^*A\hat{\mathbf{x}} = A^*\mathbf{b}$.

**Bewijs ($\Rightarrow$):** Stel $A\hat{\mathbf{x}} = \hat{\mathbf{y}} = P_\mathcal{D} \mathbf{b}$. Dan $(\mathbf{b} - A\hat{\mathbf{x}}) \perp \mathcal{D} = \mathcal{R}(A)$, dus voor alle $i$:
$$\mathbf{a}_i^*(\mathbf{b} - A\hat{\mathbf{x}}) = 0 \implies A^*(\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0} \implies A^*A\hat{\mathbf{x}} = A^*\mathbf{b}.$$

**Bewijs ($\Leftarrow$):** Uit $A^*A\hat{\mathbf{x}} = A^*\mathbf{b}$ volgt $A^*(A\hat{\mathbf{x}} - \mathbf{b}) = \mathbf{0}$, dus $(A\hat{\mathbf{x}} - \mathbf{b}) \perp \mathcal{R}(A)$. Via Bewijzen 2 en 3: $A\hat{\mathbf{x}}$ is de beste benadering. $\square$

---

## Hoofdstuk 3: Benaderen in een unitaire ruimte — Bewijzen

---

## Bewijs 5 — Genormeerde ruimte impliceert metrische ruimte

**Stelling:** $\rho(x, y) = \|x - y\|$ is een afstand.

**Bewijs:** Niet-negativiteit, nulvector en symmetrie volgen onmiddellijk uit de normeigenschappen. **Driehoeksongelijkheid:** stel $\alpha = x - z$ en $\beta = z - y$:
$$\|x - y\| = \|\alpha + \beta\| \leq \|\alpha\| + \|\beta\| = \|x - z\| + \|z - y\|. \quad \square$$

---

## Bewijs 6 — Beste benadering in strikt genormeerde ruimte

Zie **V-B, deelvraag (d)** — volledig bewijs staat daar.

---

## Bewijs 7 — Unitaire ruimte impliceert genormeerde ruimte

**Stelling:** $\|x\| = \sqrt{(x, x)}$ is een norm.

**Bewijs:** Normeigenschappen 1–3 volgen onmiddellijk. **Driehoeksongelijkheid** (voor $x+y \neq 0$) via **Cauchy-Schwarz**:
$$0 \leq (x+y, x+y) = (x, x+y) + (y, x+y)$$
$$\leq \sqrt{(x,x)}\sqrt{(x+y,x+y)} + \sqrt{(y,y)}\sqrt{(x+y,x+y)}.$$
Delen door $\|x+y\| > 0$: $\|x+y\| \leq \|x\| + \|y\|$. $\square$

---

## Bewijs 8 — Genormeerde ruimte is unitair $\Leftrightarrow$ parallellogramgelijkheid

**Stelling:** Een genormeerde vectorruimte is unitair (met scalair product dat de norm induceert) als en slechts als de norm voldoet aan de **parallellogramgelijkheid**:
$$\|x+y\|^2 + \|x-y\|^2 = 2\!\left(\|x\|^2 + \|y\|^2\right).$$

**Bewijs ($\Rightarrow$, noodzakelijkheid):**

Als de norm geïnduceerd wordt door een scalair product, dan:
$$\|x+y\|^2 + \|x-y\|^2 = (x+y,x+y) + (x-y,x-y) = 2(x,x) + 2(y,y). \quad \square$$

**Bewijs ($\Leftarrow$, voldoenheid, reëel geval):**

Men definieert het kandidaat-scalair product via de **polarisatie-identiteit**:
$$(x,y) = \tfrac{1}{4}\!\left\{\|x+y\|^2 - \|x-y\|^2\right\}.$$

Men toont aan dat dit inderdaad een scalair product is en dat de bijbehorende geïnduceerde norm gelijk is aan de oorspronkelijke norm. Het bewijs is vrij technisch en wordt in de cursus achterwege gelaten.

> ⚠️ **Examentip:** De noodzakelijkheidsrichting is examenstof en eenvoudig: expansie van de twee normen via het scalair product. De voldoendheidsrichting hoef je niet volledig te kennen. Het corollarium is belangrijk: $\|\cdot\|_p$ voor $p \neq 2$ is **geen** scalair product (parallellogramgelijkheid geldt er niet).

---

## Bewijs 9 — Eenheidsbol in unitaire ruimte is strikt convex

**Stelling:** Als $x_1 \neq x_2$ en $\|x_1\| = \|x_2\| = 1$, dan $\|x_1 + x_2\| < 2$.

**Bewijs:** Via de parallellogramgelijkheid:
$$\|x_1 + x_2\|^2 = 4 - \|x_1 - x_2\|^2 < 4 \quad (\text{want } x_1 \neq x_2). \quad \square$$

---

## Bewijs 10 — Grammatrix is HPD $\Leftrightarrow$ lineaire onafhankelijkheid

**Stelling:** $G_{jk} = (a_j, a_k)$ is Hermitisch positief definiet $\Leftrightarrow$ $\{a_1, \ldots, a_n\}$ zijn lineair onafhankelijk.

**Bewijs:** Voor willekeurige $v$:
$$v^* G v = \left\|\sum_k v_k a_k\right\|^2.$$

**($\Rightarrow$)** Lineair onafhankelijk: als $v \neq 0$ dan $\sum v_k a_k \neq 0$, dus $v^*Gv > 0$.

**($\Leftarrow$)** Lineair afhankelijk: $\exists\, \hat{v} \neq 0$ met $\sum \hat{v}_k a_k = 0$, dus $\hat{v}^*G\hat{v} = 0$ → $G$ singulier. $\square$

---

## Bewijs 11 — Orthogonale projector in unitaire ruimte: $(Pw, v) = (w, Pv)$

**Stelling:** Een projector $P$ is orthogonaal als en alleen als
$$\forall\, v, w \in V: \quad (Pw, v) = (w, Pv).$$

**Bewijs ($\Rightarrow$, noodzakelijkheid):**

Als $R(P) \perp N(P)$, dan geldt voor alle $v, w \in V$:
$$\bigl(Pw,\ (I-P)v\bigr) = 0 \quad \text{en} \quad \bigl((I-P)w,\ Pv\bigr) = 0.$$

Uitschrijven van de eerste gelijkheid: $(Pw, v) - (Pw, Pv) = 0$, dus $(Pw, v) = (Pw, Pv)$.

Uitschrijven van de tweede gelijkheid: $(w, Pv) - (Pw, Pv) = 0$, dus $(w, Pv) = (Pw, Pv)$.

Samenvoegen: $(Pw, v) = (Pw, Pv) = (w, Pv)$. $\square$

**Bewijs ($\Leftarrow$, voldoenheid):**

Neem willekeurige $x = Pu \in R(P)$ en $y \in N(P)$. Dan:
$$(x, y) = (Pu, y) = (u, Py) = (u, 0) = 0.$$

De ruimten $R(P)$ en $N(P)$ zijn dus orthogonaal, zodat $P$ een orthogonale projector is. $\square$

> ⚠️ **Examentip:** Dit is de abstracte versie van Bewijs 1 (in $\mathbb{C}^m$: $P$ orthogonaal $\Leftrightarrow$ $P = P^*$). Hier, in een unitaire ruimte, vervangt de zelfgeadjungeerdheidseis $P = P^*$ door de equivalent: $(Pw, v) = (w, Pv)$.

---

## Hoofdstuk 4: Veeltermbenaderingen

---

## Bewijs 12 — Orthogonale veeltermen vormen basis van $P_n[a,b]$

**Stelling:** $\{\phi_0, \phi_1, \ldots, \phi_n\}$ is een basis van $P_n[a,b]$.

**Bewijs:** $P_n[a,b]$ heeft dimensie $n+1$. Het stel telt $n+1$ elementen. Stel $\sum c_k \phi_k = 0$. Neem scalair product met $\phi_j$:
$$0 = c_j (\phi_j, \phi_j) \implies c_j = 0 \quad \forall j.$$
Lineair onafhankelijk → basis. $\square$

---

## Bewijs 13 — Drietermsrecursiebetrekking voor orthogonale veeltermen

**Stelling:** $\phi_k(x) = \lambda_k[(x - \alpha_k)\phi_{k-1}(x) - \beta_k\phi_{k-2}(x)]$.

**Bewijs:** Ontbind $x\phi_{k-1}(x) = \sum_{l=0}^k b_l \phi_l(x)$ met $b_l = \dfrac{(x\phi_{k-1}, \phi_l)}{(\phi_l, \phi_l)}$.

Omdat $x\phi_l$ graad $l+1$ heeft en $\phi_{k-1} \perp$ alle veeltermen van graad $< k-1$: $b_l = 0$ voor $l < k-2$.

Slechts drie termen blijven over:
$$x\phi_{k-1} = b_{k-2}\phi_{k-2} + b_{k-1}\phi_{k-1} + b_k\phi_k$$
$$\implies \phi_k = \frac{1}{b_k}\left[(x - b_{k-1})\phi_{k-1} - b_{k-2}\phi_{k-2}\right]. \quad \square$$

---

## Bewijs 14 — Orthogonale veelterm $\phi_k$ heeft $k$ enkelvoudige reële nulpunten in $(a,b)$

**Bewijs uit het ongerijmde:** Veronderstel dat $\phi_k$ slechts $m < k$ tekenwisselingen heeft in $(a,b)$, in punten $x_1, \ldots, x_m$. Dan wisselt
$$\psi(x) = \phi_k(x)\,(x-x_1)\cdots(x-x_m)$$
nergens van teken (elke tekenwisseling van $\phi_k$ wordt opgeheven door de corresponderende factor). Dus:
$$\int_a^b w(x)\,\phi_k(x)\,(x-x_1)\cdots(x-x_m)\,dx \neq 0.$$
Maar $(x-x_1)\cdots(x-x_m)$ heeft graad $m < k$, dus is deze integraal nul door orthogonaliteit van $\phi_k$ — tegenspraak. $\square$

---

## Bewijs 16 — Equioscillatiestelling

**Stelling:** De fout van de $n$-de graadsbenadering (minimaxbenadering) wordt nul in minstens $n+1$ punten van $(a,b)$.

**Bewijs per contradictie:** Veronderstel dat $r_n(x)$ slechts $m \leq n$ tekenwisselingen ondergaat op de punten $x_1, \ldots, x_m$. Dan wisselt $r_n(x) \cdot V(x)$ met $V(x) = (x-x_1)\cdots(x-x_m)$ nooit van teken, zodat:
$$\int_a^b w(x) r_n(x) V(x)\, dx \neq 0.$$

Maar $y_n$ is de beste benadering → $r_n \perp P_n[a,b]$. Omdat $V \in P_n[a,b]$ (graad $m \leq n$): de integraal moet nul zijn — **tegenspraak**. $\square$

---

## Hoofdstuk 5: Benaderen met splinefuncties — Bewijzen

---

## Bewijs 17 — $M_{i,k+1}$ en $N_{i,k+1}$ zijn splinefuncties

Via gedeelde differentie:
$$M_{i,k+1}(x) = \sum_{s=i}^{i+k+1} \lambda_s (t_s - x)_+^k.$$
Dit is een lineaire combinatie van afgeknotte-machtsfuncties van graad $k$ → splinefunctie. Hetzelfde geldt voor $N_{i,k+1} = (t_{i+k+1}-t_i) M_{i,k+1}$. $\square$

---

## Bewijs 19 — Lokale drager: $M_{i,k+1}(x) = 0$ buiten $[t_i, t_{i+k+1}]$

Voor $x \geq t_{i+k+1}$: alle $(t_s - x)_+^k = 0$.

Voor $x \leq t_i$: het $+$-teken kan weggelaten worden, zodat $M_{i,k+1}(x) = \Delta_t^{k+1}(t-x)^k$ — een differentie van orde $k+1$ van een veelterm van graad $k$ → identiek nul. $\square$

---

## Bewijs 20 — Positiviteit: $M_{i,k+1}(x) > 0$ voor $t_i < x < t_{i+k+1}$

**Inductie** op $k$ via de recursiebetrekking:
$$M_{i,k+1}(x) = \underbrace{\frac{x-t_i}{t_{i+k+1}-t_i}}_{>0} M_{i,k}(x) + \underbrace{\frac{t_{i+k+1}-x}{t_{i+k+1}-t_i}}_{>0} M_{i+1,k}(x).$$
Per inductiehypothese zijn beide termen positief op hun drager, die samen $(t_i, t_{i+k+1})$ bedekken. $\square$

---

## Bewijs 21 — Grensafgeleiden zijn nul

$$M_{i,k+1}^{(j)}(t_i) = M_{i,k+1}^{(j)}(t_{i+k+1}) = 0 \quad \text{voor } j = 0, \ldots, k-1.$$

Volgt uit de continuïteit van $M_{i,k+1}$ (splinefunctie, $C^{k-1}$) en het feit dat $M_{i,k+1} = 0$ voor $x \leq t_i$ en $x \geq t_{i+k+1}$. $\square$

---

## Bewijs 22 — De Boor-algoritme

Zie **V-C, deelvraag (c)** — het volledige bewijs staat daar.

---

## Bewijs 23 — Afgeleide van splinefunctie

**Stelling:** $s^{(r)}(x) = \sum_{i=-(k-r)}^{n-1} c_i^{(r)} N_{i,k+1-r}(x)$ met $c_i^{(r)} = (k+1-r)\dfrac{c_i^{(r-1)} - c_{i-1}^{(r-1)}}{t_{i+k+1-r} - t_i}$.

**Bewijs voor $r=1$:** Gebruik $N'_{i,k+1}(x) = k\!\left(\dfrac{N_{i,k}}{t_{i+k}-t_i} - \dfrac{N_{i+1,k}}{t_{i+k+1}-t_{i+1}}\right)$:
$$s'(x) = \sum_{i=-k}^{n-1} c_i N'_{i,k+1}(x) = \sum_{i=-(k-1)}^{n-1} \underbrace{k\frac{c_i - c_{i-1}}{t_{i+k}-t_i}}_{c_i^{(1)}} N_{i,k}(x). \quad \square$$

---

## Hoofdstuk 8: Grafen en eigenwaarden — Bewijzen

---

## Bewijs 24 — Graafpartitionering: oplossing is eigenvector van $\lambda_{N-1}$

Zie **V-D, deelvraag (c)** — het volledige bewijs staat daar.

---

## Hoofdstuk 9: Eigenwaardenalgoritmes — Bewijzen

---

## Bewijs 25 — Impliciet verschoven QR $\Leftrightarrow$ deelruimte-iteratie op $p(A)$

Zie **V-G, deelvraag (d)** — het bewijs staat daar.

---

## Hoofdstuk 11: Optimalisatie-algoritmes

---

## Bewijs 26 — Lokaal minimum via Hessiaan

**Stelling:** Als $\nabla f(x^*) = 0$ en $\nabla^2 f(x^*) > 0$, dan is $x^*$ een geïsoleerd lokaal minimum.

**Bewijs:** Uit $\nabla f(x^*) = 0$: richtingsafgeleide nul in alle richtingen. Uit $\nabla^2 f(x^*) > 0$: $s^T \nabla^2 f(x^*) s > 0$ voor alle $s \neq 0$. Via Taylor:
$$f(x^* + \epsilon s) = f(x^*) + \underbrace{\epsilon \nabla f(x^*)^T s}_{=0} + \frac{\epsilon^2}{2} s^T \nabla^2 f(x^* + \xi\epsilon s) s > f(x^*)$$
voor voldoende kleine $\epsilon$ (continuïteit van $\nabla^2 f$). $\square$

---

## Bewijs 27 — CG-orthogonaliteitseigenschappen

**Stelling:** Als $r^{(k-1)} \neq 0$:
$$r^{(k)T} r^{(j)} = 0 \quad \text{en} \quad p^{(k)T} A p^{(j)} = 0, \quad j = 1, \ldots, k-1.$$

**Bewijs door inductie op $k$:**

**Eigenschap 1** — Uit de residu-update $r^{(k)} = r^{(k-1)} - \alpha^{(k)} Ap^{(k-1)}$:
$$r^{(k)T} r^{(j)} = r^{(k-1)T} r^{(j)} - \alpha^{(k)} r^{(j)T} Ap^{(k-1)}.$$
Voor $j < k-1$: beide termen nul via inductiehypothese.
Voor $j = k-1$: nul als $\alpha^{(k)} = \dfrac{r^{(k-1)T}r^{(k-1)}}{r^{(k-1)T}Ap^{(k-1)}}$ — consistent met het algoritme.

**Eigenschap 2** — Uit $p^{(k)} = r^{(k)} + \beta^{(k)} p^{(k-1)}$:
$$p^{(k)T} Ap^{(j)} = r^{(k)T} Ap^{(j)} + \beta^{(k)} p^{(k-1)T} Ap^{(j)}.$$
Voor $j < k-1$: nul via Eigenschap 1 en inductiehypothese.
Voor $j = k-1$: nul als $\beta^{(k)} = -\dfrac{p^{(k-1)T} Ar^{(k)}}{p^{(k-1)T} Ap^{(k-1)}}$ — consistent met het algoritme. $\square$

> ⚠️ **Examentip:** De residu's zijn onderling orthogonaal; de zoekrichtingen zijn $A$-toegevoegd. Dit is wat "conjugate gradients" betekent.

---

## Bewijs 28 — CG minimaliseert de fout in de $A$-norm over de Krylov-ruimte

> **Staat al in Modelvragen_Samaey.md (§ CG, deelvraag c) — volledig bewijs aanwezig.** Hieronder een samenvatting.

**Stelling:** Als $r^{(k-1)} \neq 0$, dan geldt:
$$\|e^{(k)}\|_A^2 = \min_{x \in \mathcal{K}_k(A,b)} \|x - x^*\|_A^2.$$

**Bewijs:**

Schrijf elke $x \in \mathcal{K}_k(A,b)$ als $x = x^{(k)} + \delta x$ met $\delta x \in \mathcal{K}_k(A,b)$:
$$\|x - x^*\|_A^2 = \|\delta x + e^{(k)}\|_A^2 = \|e^{(k)}\|_A^2 + 2\, e^{(k)T} A\, \delta x + \|\delta x\|_A^2.$$

Merk op dat $Ae^{(k)} = Ax^{(k)} - b = -r^{(k)}$. Uit Bewijs 27 geldt $r^{(k)} \perp \mathcal{K}_k(A,b)$, dus:
$$e^{(k)T} A\, \delta x = -r^{(k)T} \delta x = 0.$$

De kruisterm is altijd nul en $\|\delta x\|_A^2 \geq 0$. Het minimum wordt bereikt voor $\delta x = 0$, d.w.z. $x = x^{(k)}$. $\square$

**Gevolg:** $\|e^{(k)}\|_A \leq \|e^{(k-1)}\|_A$ (monotone daling) en $e^{(n)} = 0$ voor een $n \leq m$.

---

## Bewijs 29 — CG als veeltermbenadering

**Stelling:** Als $r^{(k-1)} \neq 0$, dan geldt:
$$\|e^{(k)}\|_A^2 = \min_{p \in \mathcal{N}_k} \|p(A)\, e^{(0)}\|_A^2,$$
met $\mathcal{N}_k = \{p \in \mathbb{R}[t] : \deg p \leq k,\ p(0) = 1\}$.

**Bewijs:**

Elke $x \in \mathcal{K}_k(A,b)$ kan geschreven worden als:
$$x = c_1 b + c_2 Ab + \cdots + c_k A^{k-1}b$$
voor bepaalde $c_1, \ldots, c_k \in \mathbb{R}$. Omdat $b = Ax^* = -Ae^{(0)}$ (want $e^{(0)} = x^{(0)} - x^*$ en $x^{(0)} = 0$ impliceert $b = -Ae^{(0)}$), geldt:
$$x - x^* = \bigl(I - c_1 A - c_2 A^2 - \cdots - c_k A^k\bigr) e^{(0)} =: \hat{p}(A)\, e^{(0)},$$
waarbij $\hat{p}(t) = 1 - c_1 t - c_2 t^2 - \cdots - c_k t^k \in \mathcal{N}_k$ (graad $\leq k$ en $\hat{p}(0) = 1$). Omgekeerd geeft elke $\hat{p} \in \mathcal{N}_k$ een $x = x^* + \hat{p}(A)e^{(0)} \in \mathcal{K}_k(A,b)$.

Er is dus een bijectie tussen $\mathcal{K}_k(A,b)$ en $\mathcal{N}_k$, en het gestelde volgt direct uit Bewijs 28:
$$\|e^{(k)}\|_A^2 = \min_{x \in \mathcal{K}_k(A,b)} \|x - x^*\|_A^2 = \min_{\hat{p} \in \mathcal{N}_k} \|\hat{p}(A)\, e^{(0)}\|_A^2. \quad \square$$

> ⚠️ **Examentip:** Dit verband verklaart de convergentie van CG via het spectrum van $A$. Zijn alle eigenwaarden geclusterd of zijn er slechts $k$ verschillende eigenwaarden, dan convergeert CG snel (in $\leq k$ stappen). Het Chebyshev-veelterm-argument geeft de scherpste bovengrens op de convergentiesnelheid.

---

## Hoofdstuk 12: Ijle representaties — Bewijzen

---

## Bewijs 30 — Beste rang-$k$ benadering (Eckart–Young–Mirsky)

> **Spectraalnorm-bewijs staat ook in Modelvragen_WM.md (§ rang-k benadering).** Hieronder het volledige bewijs voor beide normen.

**Stelling:** De getrunceerde SVD $A_k^{svo} = \sum_{i=1}^k \sigma_i u_i v_i^*$ voldoet aan:
$$\|A_k^{svo} - A\|_2 = \min_{B \in \mathcal{M}_k} \|B - A\|_2 \qquad \text{en} \qquad \|A_k^{svo} - A\|_F = \min_{B \in \mathcal{M}_k} \|B - A\|_F,$$
waarbij $\mathcal{M}_k$ de verzameling van alle rang-$\leq k$ matrices is.

---

### Bewijs voor de spectraalnorm

Neem een willekeurige $B \in \mathcal{M}_k$. Schrijf $B = XY^*$ met $X \in \mathbb{C}^{m \times k}$, $Y \in \mathbb{C}^{n \times k}$.

**Sleutelargument:** De ruimte $\langle v_1, \ldots, v_{k+1} \rangle$ heeft dimensie $k+1$ en $\mathcal{N}(Y^*)$ heeft dimensie $n - k$. Omdat $(k+1) + (n-k) = n+1 > n$, bestaat er een $w$ met $\|w\|_2 = 1$ zodat:
$$w \in \langle v_1, \ldots, v_{k+1} \rangle \cap \mathcal{N}(Y^*).$$

Schrijf $w = \sum_{i=1}^{k+1} c_i v_i$ met $\sum_{i=1}^{k+1} c_i^2 = 1$. Omdat $w \in \mathcal{N}(Y^*)$ geldt $Bw = XY^*w = 0$, dus:
$$\|B - A\|_2 \geq \|(B-A)w\|_2 = \|Aw\|_2 = \left\|\sum_{i=1}^{k+1} c_i \sigma_i u_i\right\|_2 = \sqrt{\sum_{i=1}^{k+1} c_i^2 \sigma_i^2} \geq \sigma_{k+1} \underbrace{\sqrt{\sum_{i=1}^{k+1} c_i^2}}_{=1} = \sigma_{k+1}.$$

Vergelijking met $\|A_k^{svo} - A\|_2 = \sigma_{k+1}$ toont dat $B$ niet beter kan zijn. $\square$

---

### Bewijs voor de Frobeniusnorm

We tonen eerst de hulpstelling $\sigma_i(B - A) \geq \sigma_{k+i}(A)$ voor $i = 1, \ldots, p-k$ (met $p = \min(m,n)$).

Noteer met superscript $\mathrm{svo}$ de beste rang-$\ell$ benadering in de spectraalnorm. Dan:
$$\sigma_i(B - A) = \|(A - B) - (A - B)^{\mathrm{svo}}_{i-1}\|_2 = \|A - \underbrace{(B + (A-B)^{\mathrm{svo}}_{i-1})}_{\text{rang} \leq k + (i-1)}\|_2 \geq \|A - A^{\mathrm{svo}}_{k+i-1}\|_2 = \sigma_{k+i}(A).$$

Hieruit volgt:
$$\|B - A\|_F^2 = \sum_{i=1}^p \sigma_i(B-A)^2 \geq \sum_{i=1}^{p-k} \sigma_i(B-A)^2 \geq \sum_{i=1}^{p-k} \sigma_{k+i}(A)^2 = \sum_{i=k+1}^p \sigma_i(A)^2 = \|A_k^{svo} - A\|_F^2.$$

Bijgevolg is $B$ geen betere Frobenius-benadering dan $A_k^{svo}$. $\square$

> ⚠️ **Examentip:** De spectraalnormbewijzen vragen dimensietelling om een geschikte $w$ te vinden. Het Frobenius-bewijs gebruikt de hulpstelling $\sigma_i(B-A) \geq \sigma_{k+i}(A)$ recursief.

---