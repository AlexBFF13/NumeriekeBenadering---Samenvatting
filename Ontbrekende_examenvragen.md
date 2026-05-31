# Ontbrekende Examenvragen — Uitgewerkte Antwoorden

> Gebaseerd op: *Examenvragen tot 2025 (Maxim Bijnens)* en *proofs (1).pdf*
> Enkel volledig ontbrekende examenvragen; deelvragen bij bestaande vragen komen later.

---

## DEEL 1 — Ontbrekende Examenvragen

---

## Hoofdstuk 2 & 3: Benadering van vectoren / Benaderen in een unitaire ruimte

---

## V-A — QR-factorisatie: vergelijking van drie methodes

**Vraag (H2/H3):** Bespreek de QR-factorisatie en hoe je die kan opstellen.
- (a) Bespreek de geneste structuur van de QR-factorisatie en waarom deze structuur kan gebruikt worden in de Arnoldi-iteratie.
- (d) Bespreek kort hoe men met Givensrotaties een QR-factorisatie kan opstellen.
- (f) Vergelijk de drie methodes op basis van stabiliteit, orthogonaliteit en rekenkost. Bij welk type matrices verkies je elke methode?

---

### (a) Geneste structuur van de QR-factorisatie

De QR-factorisatie $A = QR$ (of $A = \hat{Q}\hat{R}$ in de onvolledige versie) heeft een **geneste structuur**:

$$\langle q_1, \ldots, q_k \rangle = \langle a_1, \ldots, a_k \rangle, \quad k = 1, \ldots, n.$$

Dit betekent dat na $k$ stappen van het Gram-Schmidt algoritme de ruimte opgespannen door de eerste $k$ orthonormale vectoren gelijk is aan de ruimte opgespannen door de eerste $k$ kolomvectoren van $A$.

**Verband met Arnoldi:** De Arnoldi-iteratie bouwt stap per stap een orthonormale basis van de Krylov-ruimte $\mathcal{K}_k(A,b) = \langle b, Ab, \ldots, A^{k-1}b \rangle$ op. Precies zoals bij QR groeit de basis stap voor stap: na iteratie $k$ heeft men een orthonormale basis $\{q_1, \ldots, q_k\}$ voor $\mathcal{K}_k$. De geneste structuur is essentieel: de basis voor $\mathcal{K}_{k+1}$ bevat de volledige basis voor $\mathcal{K}_k$. Dit laat toe om de Arnoldi-iteratie on-the-fly uit te breiden zonder vroegere berekeningen te herhalen.

> ⚠️ **Examentip:** De Krylov-ruimte heeft de geneste structuur $\mathcal{K}_1 \subset \mathcal{K}_2 \subset \cdots$, precies zoals de geprojecteerde deelruimten bij QR. Arnoldi is dus een QR-factorisatie van de matrix $[b, Ab, A^2b, \ldots]$.

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

**Stelling:** Zij $\mathcal{D}$ een eindigdimensionale deelruimte van een strikt genormeerde ruimte $V$ en $\vec{v} \in V$. Dan bestaat de beste benadering van $\vec{v}$ in $\mathcal{D}$ en is deze uniek.

#### Existentie

Definieer $d = \inf\{\|\vec{v} - \vec{w}\| : \vec{w} \in \mathcal{D}\}$. Volgens de definitie van infimum bestaat er een rij $\{\vec{w}_k\}$ in $\mathcal{D}$ zodat $\|\vec{v} - \vec{w}_k\| \to d$.

De rij is uniform begrensd: $\|\vec{w}_k\| \leq \|\vec{w}_k - \vec{v}\| + \|\vec{v}\| \leq \|\vec{w}_1 - \vec{v}\| + \|\vec{v}\|$.

Stel $n = \dim(\mathcal{D})$ en beschouw een basis $\{\vec{a}_1, \ldots, \vec{a}_n\}$. Dan $\vec{w}_k = \sum c_{ki} \vec{a}_i$ en de rij $\{(c_{k1},\ldots,c_{kn})\}$ is begrensd. Volgens **Weierstrass-Bolzano** bestaat een convergente deelrij met limiet $(\hat{c}_1, \ldots, \hat{c}_n)$. Definieer $\vec{\zeta} = \sum \hat{c}_i \vec{a}_i \in \mathcal{D}$. Dan:
$$\|\vec{v} - \vec{\zeta}\| \leq \underbrace{\|\vec{v} - \vec{w}_k\|}_{\to d} + \underbrace{\|\vec{w}_k - \vec{\zeta}\|}_{\to 0} \to d.$$
Dus $\|\vec{v} - \vec{\zeta}\| = d$: de beste benadering bestaat.

#### Uniciteit

Bewijs uit het ongerijmde. Veronderstel twee beste benaderingen $\vec{v}_1 \neq \vec{v}_2$ met $\|\vec{v} - \vec{v}_1\| = \|\vec{v} - \vec{v}_2\| = d$.

Definieer $\vec{e}_i = \frac{1}{d}(\vec{v} - \vec{v}_i)$ op de eenheidsbol, $\vec{e}_1 \neq \vec{e}_2$. Omdat de norm **strikt** is:
$$\left\|\vec{v} - \frac{\vec{v}_1 + \vec{v}_2}{2}\right\| = d \left\|\frac{1}{2}\vec{e}_1 + \frac{1}{2}\vec{e}_2\right\| < d.$$

Maar $\frac{\vec{v}_1+\vec{v}_2}{2} \in \mathcal{D}$ is een **betere** benadering dan $\vec{v}_1$ — tegenspraak. $\square$

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

**Partitievector:** $x \in \{-1, +1\}^N$ met $x_i = +1$ als $i \in K_1$ en $x_i = -1$ als $i \in K_2$.

**Gelijke deelgrafen:** $\mathbf{1}^T x = 0$.

**Aantal doorgeknippen bogen:**
$$\frac{1}{4} x^T L x = \frac{1}{4} \sum_{i,j} L_{ij} x_i x_j = \frac{1}{2} \sum_{(i,j) \in B} (-x_i x_j + x_i^2) = \#(B_x).$$

**Optimaliseringsprobleem:**
$$\min_{x \in \{-1,+1\}^N} \frac{1}{4} x^T L x \quad \text{o.v.v.} \quad \mathbf{1}^T x = 0.$$

---

### (b) Relaxatie

De binaire beperking maakt het probleem NP-hard. **Relaxatie:**
$$x \in \{-1,+1\}^N \;\longrightarrow\; \|x\|_2^2 = N.$$

Het gerelaxeerde probleem:
$$\min_{x \in \mathbb{R}^N} \frac{1}{4} x^T L x \quad \text{o.v.v.} \quad \|x\|_2^2 = N,\; \mathbf{1}^T x = 0.$$

---

### (c) Oplossing is eigenvector van $\lambda_{N-1}$ (Fiedler-vector)

**Bewijs:**

Omdat $L$ symmetrisch is: $L = V\Lambda V^T$, met $V$ orthogonaal.

**Eigenwaarde 0:** $L\mathbf{1} = \mathbf{0}$, dus $\lambda_N = 0$ met $v_N = \frac{1}{\sqrt{N}}\mathbf{1}$.

**Substitutie** $x = Vy$: de eis $\mathbf{1}^T x = 0$ wordt $v_N^T Vy = [0 \cdots 0 \; 1]y = y_N = 0$.

Het probleem transformeert:
$$\min \frac{1}{4} \sum_{i=1}^N \lambda_i y_i^2 \quad \text{o.v.v.} \quad \|y\|_2^2 = N,\; y_N = 0.$$

Om $\sum \lambda_i y_i^2$ te minimaliseren onder $\sum y_i^2 = N$ en $y_N = 0$: geef al het gewicht aan de kleinste vrije eigenwaarde $\lambda_{N-1}$:
$$y = [0 \cdots 0 \; \sqrt{N} \; 0]^T \implies x = Vy = \sqrt{N}\, v_{N-1}. \quad \square$$

**Praktische partitie:** $x_i = \text{sign}\bigl((v_{N-1})_i\bigr)$.

---

### (d) Algoritme: Lanczos

Omdat $L$ **symmetrisch en spaars** is en we slechts de tweede-kleinste eigenwaarde nodig hebben, gebruiken we het **Lanczos-algoritme** (variant van Arnoldi voor symmetrische matrices):
- Enkel matrix-vector producten $Lv$ vereist ($O(nnz)$ per stap).
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

**Rank-definitie:** Een pagina is centraal als ze veel inkomende links heeft van zelf centrale pagina's. De **link-matrix** $\hat{A}$:
$$\hat{A}_{ij} = \frac{1}{k_j} \quad \text{als } j \to i \text{ (met } k_j \text{ totale uitlinks van } j\text{)}, \quad 0 \text{ anders.}$$

De rank voldoet aan $\hat{A} r = r$ (eigenvector met eigenwaarde 1).

---

### (b) Benaderend eigenwaardeprobleem

De rank-vector $r$ is een **eigenvector van $\hat{A}$ met eigenwaarde 1**.

**Problemen:**
1. **Dangling nodes:** Pagina's zonder uitlinks → nulkolom → matrix niet stochastisch.
2. **Niet-samenhangende graaf:** Eigenwaarde 1 niet uniek.

---

### (c) Omvorming tot volwaardig eigenwaardeprobleem

**Stap 1 — Dangling nodes:** Vervang nulkolommen door $\frac{1}{N}\mathbf{1}$:
$$\tilde{A}_{ij} = \hat{A}_{ij} + \frac{d_j}{N}, \quad d_j = 1 \text{ als } j \text{ geen uitlinks heeft, anders } 0.$$

**Stap 2 — Google-matrix (teleportatie):**
$$G = (1-\alpha)\tilde{A} + \frac{\alpha}{N}\mathbf{1}\mathbf{1}^T, \quad \alpha \approx 0.15.$$

$G$ is nu een **primitieve stochastische matrix**: eigenwaarde 1 is enkelvoudig en dominant (Perron-Frobenius). De aanpassingen hebben weinig effect op de ranking ($\alpha$ is klein).

---

### (d) Algoritme: methode van de machten

```
r^(0) = (1/N) * 1       ← uniforme beginvector
herhaal:
    r^(k+1) = G * r^(k) ← matrix-vector product (O(nnz) voor ijle matrices)
    normaliseer r^(k+1)
totdat convergentie
```

**Convergentiefactor:** $|\lambda_2(G)| \leq 1 - \alpha$ — typisch $50$–$100$ iteraties.

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

**Doelstelling:** Bereken de dominante eigenwaarde $\lambda_1$ (grootste $|\lambda_i|$) van $A$.

**Basisiteratie:**
```
kies b^(0) met ||b^(0)||₂ = 1
voor k = 1, 2, ...:
    z^(k) = A b^(k-1)
    b^(k) = z^(k) / ||z^(k)||₂
```

**Convergentiebewijs:**

Schrijf $b^{(0)} = \sum_i c_i x_i$ met $\{x_i\}$ de eigenbasis. Dan:
$$A^k b^{(0)} = \lambda_1^k \left(c_1 x_1 + \sum_{i \geq 2} c_i \left(\frac{\lambda_i}{\lambda_1}\right)^k x_i\right).$$

Als $|\lambda_1| > |\lambda_2|$ en $c_1 \neq 0$: $\left(\frac{\lambda_i}{\lambda_1}\right)^k \to 0$ voor $i \geq 2$, dus $b^{(k)} \to \frac{x_1}{\|x_1\|}$.

**Convergentiefactor:** $\left|\dfrac{\lambda_2}{\lambda_1}\right|$ — lineaire convergentie.

---

### (b) Dominante eigenwaarde schatten

Via het **Rayleigh-quotiënt**:
$$\mu^{(k)} = (b^{(k)})^* A b^{(k)} \xrightarrow{k\to\infty} \lambda_1.$$

De convergentie is **kwadratisch** voor Hermitische matrices.

---

### (c) Kleinste modulus / dichtst bij $\sigma$

**Kleinste eigenwaarde (inverse iteratie):** Pas de methode van de machten toe op $A^{-1}$. De dominante eigenwaarde van $A^{-1}$ is $\lambda_m^{-1}$ (kleinste modulus van $A$). In elke stap: los $Az^{(k)} = b^{(k-1)}$ op via LU-factorisatie.

**Eigenwaarde dichtst bij $\sigma$ (verschuiving):** Pas inverse iteratie toe op $(A - \sigma I)^{-1}$. Hoe dichter $\sigma$ bij $\lambda_j$, hoe sneller de convergentie:
$$\text{convergentiefactor} = \left|\frac{\lambda_j - \sigma}{\lambda_{\text{tweede dichtste}} - \sigma}\right|.$$

---

### (d) Deelruimte-iteratie

**Doelstelling:** Bereken de $n$ dominante eigenvectoren van $A$.

**Algoritme:**
```
kies Z^(0) ∈ ℂ^(m×n) met orthonormale kolommen
voor k = 1, 2, ...:
    Z̃^(k) = A Z^(k-1)            ← matrix-matrixproduct
    Z̃^(k) = Z^(k) R^(k)          ← QR-factorisatie (herorthogonalisatie)
```

**Convergentie:** De kolomruimte van $Z^{(k)}$ convergeert naar $\langle x_1, \ldots, x_n \rangle$.

**Convergentiefactor:** $\left|\dfrac{\lambda_{n+1}}{\lambda_n}\right|^k$.

---

### (e) Eigenwaarden en -vectoren verkrijgen

Na convergentie: **Rayleigh-Ritz projectie**:
$$H^{(k)} = (Z^{(k)})^* A Z^{(k)} \in \mathbb{C}^{n \times n}.$$

Bereken eigenwaarden van $H^{(k)}$ via het volledige QR-algoritme (klein stelsel). De eigenvectoren van $A$ zijn $Z^{(k)} y_i$ met $y_i$ de eigenvectoren van $H^{(k)}$.

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

**Doel:** Reduceer de dimensie van een dataset door de belangrijkste richtingen te behouden.

**Procedure:**
1. **Datamatrix** $X \in \mathbb{R}^{n \times N}$ (gecentreerd: $\sum x^{(i)} = 0$).
2. **SVD:** $X = U\Sigma V^T$.
3. **Principale componenten** $r_i = v_i$ (rechtersinguliere vectoren — richtingen van maximale variantie).
4. **Projectie op $k$ componenten:** $\tilde{X} = V_k^T X \in \mathbb{R}^{k \times N}$ (lage-dimensionale representatie).
5. **Reconstructie:** $\hat{X} = V_k\tilde{X} = U_k\Sigma_k V_k^T = X_k^{svo}$ (getrunceerde SVD).

**Reconstructiefout:** $\|X - \hat{X}\|_F^2 = \sigma_{k+1}^2 + \cdots + \sigma_r^2$.

**Motivatie:** Fysische systemen zijn inherent laag-dimensionaal; de data liggen nabij een laag-dimensionale manifold.

---

## DEEL 2 — Ontbrekende Bewijzen

---

## Hoofdstuk 2: Benadering van vectoren

---

## Bewijs 2 — Karakterisatie van de beste benadering

**Stelling:** $\hat{y}$ is een beste benadering in deelruimte $\mathcal{D}$ voor $b$ als en alleen als $b - \hat{y} \perp \mathcal{D}$.

**Bewijs ($\Rightarrow$):** Neem willekeurige $y \in \mathcal{D}$. Vermits $y - \hat{y} \in \mathcal{D}$ en $b - \hat{y} \perp \mathcal{D}$, staan $(b-\hat{y})$ en $(y-\hat{y})$ loodrecht. Via **Pythagoras**:
$$\|b - y\|_2^2 = \|b - \hat{y}\|_2^2 + \|\hat{y} - y\|_2^2 \geq \|b - \hat{y}\|_2^2.$$
Dus $\hat{y}$ benadert $b$ minstens even goed als elke $y \in \mathcal{D}$.

**Bewijs ($\Leftarrow$):** Als $\hat{y}$ een beste benadering is maar $b - \hat{y} \not\perp \mathcal{D}$, bestaat $z \in \mathcal{D}$ met $(b-\hat{y}, z) \neq 0$. Dan biedt $\hat{y} + tz$ voor gepaste $t$ een kleinere afstand — tegenspraak. $\square$

---

## Bewijs 3 — Beste benadering bestaat en is uniek via $P_\mathcal{D}$

**Stelling:** De beste benadering in $\mathcal{D}$ voor $b$ is $\hat{y} = P_\mathcal{D} b$.

**Bewijs:** Ontbind $b = P_\mathcal{D} b + (I - P_\mathcal{D}) b$. Dan:
$$b - \hat{y} = (I - P_\mathcal{D}) b \in \mathcal{N}(P_\mathcal{D}) = \mathcal{D}^\perp \implies b - \hat{y} \perp \mathcal{D}.$$
Uit Bewijs 2: $\hat{y}$ is een beste benadering.

**Uniciteit:** Als ook $y \in \mathcal{D}$ een beste benadering is:
$$\|b - y\|_2^2 = \|b - \hat{y}\|_2^2 + \|y - \hat{y}\|_2^2 \geq \|b - \hat{y}\|_2^2,$$
met gelijkheid enkel als $\|y - \hat{y}\|_2 = 0$, dus $y = \hat{y}$. $\square$

---

## Bewijs 4 — Normaalstelsel

**Stelling:** $\hat{x}$ is oplossing van $\min\|Ax-b\|_2$ $\Leftrightarrow$ $A^*A\hat{x} = A^*b$.

**Bewijs ($\Rightarrow$):** Stel $A\hat{x} = \hat{y} = P_\mathcal{D} b$. Dan $(b - A\hat{x}) \perp \mathcal{D} = \mathcal{R}(A)$, dus voor alle $i$:
$$a_i^*(b - A\hat{x}) = 0 \implies A^*(b - A\hat{x}) = 0 \implies A^*A\hat{x} = A^*b.$$

**Bewijs ($\Leftarrow$):** Uit $A^*A\hat{x} = A^*b$ volgt $A^*(A\hat{x} - b) = 0$, dus $(A\hat{x} - b) \perp \mathcal{R}(A)$. Via Bewijzen 2 en 3: $A\hat{x}$ is de beste benadering. $\square$

---

## Hoofdstuk 3: Benaderen in een unitaire ruimte — Bewijzen

---

## Bewijs 5 — Genormeerde ruimte impliceert metrische ruimte

**Stelling:** $\rho(\vec{x}, \vec{y}) = \|\vec{x} - \vec{y}\|$ is een afstand.

**Bewijs:** Niet-negativiteit, nulvector en symmetrie volgen onmiddellijk uit de normeigenschappen. **Driehoeksongelijkheid:** stel $\vec{\alpha} = \vec{x} - \vec{z}$ en $\vec{\beta} = \vec{z} - \vec{y}$:
$$\|\vec{x} - \vec{y}\| = \|\vec{\alpha} + \vec{\beta}\| \leq \|\vec{\alpha}\| + \|\vec{\beta}\| = \|\vec{x} - \vec{z}\| + \|\vec{z} - \vec{y}\|. \quad \square$$

---

## Bewijs 6 — Beste benadering in strikt genormeerde ruimte

Zie **V-B, deelvraag (d)** — volledig bewijs staat daar.

---

## Bewijs 7 — Unitaire ruimte impliceert genormeerde ruimte

**Stelling:** $\|\vec{x}\| = \sqrt{(\vec{x}, \vec{x})}$ is een norm.

**Bewijs:** Normeigenschappen 1–3 volgen onmiddellijk. **Driehoeksongelijkheid** (voor $\vec{x}+\vec{y} \neq \vec{0}$) via **Cauchy-Schwarz**:
$$0 \leq (\vec{x}+\vec{y}, \vec{x}+\vec{y}) = (\vec{x}, \vec{x}+\vec{y}) + (\vec{y}, \vec{x}+\vec{y})$$
$$\leq \sqrt{(\vec{x},\vec{x})}\sqrt{(\vec{x}+\vec{y},\vec{x}+\vec{y})} + \sqrt{(\vec{y},\vec{y})}\sqrt{(\vec{x}+\vec{y},\vec{x}+\vec{y})}.$$
Delen door $\|\vec{x}+\vec{y}\| > 0$: $\|\vec{x}+\vec{y}\| \leq \|\vec{x}\| + \|\vec{y}\|$. $\square$

---

## Bewijs 9 — Eenheidsbol in unitaire ruimte is strikt convex

**Stelling:** Als $\vec{x}_1 \neq \vec{x}_2$ en $\|\vec{x}_1\| = \|\vec{x}_2\| = 1$, dan $\|\vec{x}_1 + \vec{x}_2\| < 2$.

**Bewijs:** Via de parallellogramgelijkheid:
$$\|\vec{x}_1 + \vec{x}_2\|^2 = 4 - \|\vec{x}_1 - \vec{x}_2\|^2 < 4 \quad (\text{want } \vec{x}_1 \neq \vec{x}_2). \quad \square$$

---

## Bewijs 10 — Grammatrix is HPD $\Leftrightarrow$ lineaire onafhankelijkheid

**Stelling:** $G_{jk} = (\vec{a}_j, \vec{a}_k)$ is Hermitisch positief definiet $\Leftrightarrow$ $\{\vec{a}_1, \ldots, \vec{a}_n\}$ zijn lineair onafhankelijk.

**Bewijs:** Voor willekeurige $v$:
$$v^* G v = \left\|\sum_k v_k \vec{a}_k\right\|^2.$$

**($\Rightarrow$)** Lineair onafhankelijk: als $v \neq 0$ dan $\sum v_k \vec{a}_k \neq 0$, dus $v^*Gv > 0$.

**($\Leftarrow$)** Lineair afhankelijk: $\exists\, \hat{v} \neq 0$ met $\sum \hat{v}_k \vec{a}_k = 0$, dus $\hat{v}^*G\hat{v} = 0$ → $G$ singulier. $\square$

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