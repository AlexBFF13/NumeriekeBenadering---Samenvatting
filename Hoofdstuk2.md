# Deel 2 – Benadering van vectoren

## Overzicht

Dit hoofdstuk behandelt één centrale vraag: *hoe vind je de beste benadering van een vector $\mathbf{b}$ in een deelruimte $\mathcal{D}$?* Dat klinkt abstract, maar het is precies wat er achter curve fitting, signaaldecompositie en kleinste-kwadratenregressie zit.

De rode draad is het gebruik van **orthogonale basissen**. Ze maken berekeningen eenvoudiger, coëfficiënten onafhankelijker van de keuze van deelruimte, en algoritmen numeriek stabieler. De vier secties bouwen logisch op elkaar:

1. **Waarom orthogonale basissen?** (§2.1)
2. **Wat zijn projectoren?** (§2.2) — de sleuteloperatie voor benadering
3. **Hoe bouw je een orthogonale basis?** (§2.3) — Gram-Schmidt, Givens, Householder
4. **Hoe bereken je de beste benadering?** (§2.4) — normaalstelsel vs. QR-aanpak

De resultaten zijn geldig voor $\mathbb{C}^m$ (complex) maar zijn rechtstreeks toepasbaar op $\mathbb{R}^m$.

---

## 2.1 Orthogonale vs. scheve basissen

### Basisnotatie

Een **deelruimte** $\mathcal{D} \subset \mathbb{C}^m$ van dimensie $n$ heeft een basis $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$, en elke vector $\mathbf{v} \in \mathcal{D}$ kan geschreven worden als:

$$\mathbf{v} = A\mathbf{c} = \sum_{i=1}^n c_i \mathbf{a}_i, \quad A = [\mathbf{a}_1 \cdots \mathbf{a}_n] \in \mathbb{C}^{m \times n}.$$

De **grammatrix** van een basis is:

$$G = A^*A = \begin{pmatrix} \mathbf{a}_1^*\mathbf{a}_1 & \cdots & \mathbf{a}_1^*\mathbf{a}_n \\ \vdots & \ddots & \vdots \\ \mathbf{a}_n^*\mathbf{a}_1 & \cdots & \mathbf{a}_n^*\mathbf{a}_n \end{pmatrix}. \tag{2.3}$$

Voor twee vectoren $\mathbf{v} = A\mathbf{c}$ en $\mathbf{w} = A\mathbf{d}$ in $\mathcal{D}$ geldt:

$$(\mathbf{v}, \mathbf{w}) = \mathbf{v}^*\mathbf{w} = \mathbf{c}^* G \mathbf{d}.$$

### Drie voordelen van een orthogonale basis

> ⚠️ **Belangrijk voor examen:** Ken deze drie voordelen en kun je ze uitleggen/illustreren.

**Voordeel 1 — Inwendig product en norm via coëfficiënten:**

Als de basis **orthonormaal** is ($G = I_n$), dan geldt:

$$(\mathbf{v}, \mathbf{w}) = \mathbf{c}^*\mathbf{d}, \qquad \|\mathbf{v}\|_2^2 = \|\mathbf{c}\|_2^2. \tag{2.4}$$

Het inwendig product en de norm kunnen volledig bepaald worden op het niveau van de coëfficiëntenvectoren — zonder de basisvectoren zelf te kennen.

**Voordeel 2 — Eenvoudig berekenen van coëfficiënten:**

De coëfficiënten worden bepaald door het normaalstelsel $G\mathbf{c} = A^*\mathbf{v}$ (zie 2.5). Bij een **orthogonale basis** wordt $G$ diagonaal en valt het stelsel uiteen in onafhankelijke vergelijkingen:

$$c_j = \frac{\mathbf{a}_j^*\mathbf{v}}{\|\mathbf{a}_j\|_2^2}, \quad j = 1, \ldots, n. \tag{2.5'}$$

Bij een **orthonormale** basis is het nog eenvoudiger: $c_j = \mathbf{q}_j^* \mathbf{v}$.

Dit heeft een belangrijk gevolg: als je de deelruimte uitbreidt (een extra basisvector toevoegt), blijven alle eerder berekende coëfficiënten ongewijzigd. Bij een scheve basis moet je alles herberekenen.

**Voordeel 3 — Numerieke stabiliteit:**

Een scheve basis met vectoren die bijna in dezelfde richting wijzen geeft een slecht geconditioneerd probleem.

### Voorbeeld: conditiegetal explosioneert bij scheve basis

Beschouw $\mathcal{D} = \{(v_1,v_2,v_3) \in \mathbb{R}^3 : v_2 = v_3\}$ met basis:

$$\mathbf{a}_1 = \begin{pmatrix}1\\0\\0\end{pmatrix}, \quad \mathbf{a}_2 = \begin{pmatrix}1\\\epsilon\\\epsilon\end{pmatrix}, \quad \epsilon > 0 \text{ klein}.$$

Voor $\mathbf{v} = [1\ 1\ 1]^T$ zijn de coëfficiënten: $c_1 = 1 - 1/\epsilon$, $c_2 = 1/\epsilon$. Voor $\epsilon = 10^{-4}$ is $c_1 = -9999$ en $c_2 = 10000$: twee enorm grote getallen met tegengesteld teken. Een relatieve fout van $1\%$ op $c_1$ geeft een relatieve fout van bijna $100\%$ op $\mathbf{v}$ — **gevaarlijke aftrekking**.

Het conditiegetal $\kappa(A) \approx 1/(\sqrt{2}\,\epsilon)$ en $\kappa(G) = \kappa(A)^2$ exploderen voor kleine $\epsilon$:

| $\epsilon$ | $\kappa(A)$ | $\kappa(G)$ |
|---|---|---|
| $10^{-1}$ | $1.4 \times 10^1$ | $2.0 \times 10^2$ |
| $10^{-2}$ | $1.4 \times 10^2$ | $2.0 \times 10^4$ |
| $10^{-4}$ | $1.4 \times 10^4$ | $2.0 \times 10^8$ |

De orthonormale basis $\{\mathbf{q}_1, \mathbf{q}_2\}$ van dezelfde ruimte heeft $\kappa([\mathbf{q}_1\ \mathbf{q}_2]) = 1$ en $\kappa(G) = 1$ — perfect geconditioneerd.

---

## 2.2 Begrip projector

### Definitie

> ⚠️ **Belangrijk voor examen:** Weet de definitie, eigenschappen en de karakterisering van orthogonale projectoren.

### Definitie
Een **projector** is een matrix $P \in \mathbb{C}^{m \times m}$ die **idempotent** is:
$$P^2 = P.$$

**Intuïtie:** Als je een vector eenmaal projecteert en daarna nog een keer, kom je op hetzelfde punt uit — projecteren tweemaal is hetzelfde als eenmaal.

### Eigenschappen van projectoren

Voor een projector $P$ geldt:

- Elke $\mathbf{v} \in \mathbb{C}^m$ wordt uniek ontbonden als: $\mathbf{v} = \underbrace{P\mathbf{v}}_{\in \mathcal{R}(P)} + \underbrace{(I-P)\mathbf{v}}_{\in \mathcal{N}(P)}$
- $\mathcal{R}(P) \cap \mathcal{N}(P) = \{0\}$
- $\dim(\mathcal{R}(P)) + \dim(\mathcal{N}(P)) = m$
- Als $\mathbf{v} \in \mathcal{R}(P)$, dan is $P\mathbf{v} = \mathbf{v}$ (projectie op eigen subspace geeft identiteit)
- $\tilde{P} = I - P$ is ook een projector (de **complementaire projector**)

**Geometrische betekenis:** $P$ projecteert *op* $\mathcal{R}(P)$, *langs* de richting van $\mathcal{N}(P)$.

### Constructieformule voor een scheve projector

Gegeven deelruimten $S_1$ (dimensie $n$, basis in kolommen van $A$) en $S_2$ (het orthogonaal complement van de projectierichting, basis in kolommen van $B$), is de projector op $S_1$ langs $S_2$:

$$P_{S_1, S_2} = A(B^*A)^{-1}B^*. \tag{2.12}$$

### Orthogonale projector

### Definitie
Een projector $P$ is **orthogonaal** als $\mathcal{R}(P) \perp \mathcal{N}(P)$.

### Stelling
**Een projector $P$ is orthogonaal als en alleen als $P = P^*$.**

### Bewijs

**($\Rightarrow$) De voorwaarde is nodig:**
Stel dat $P$ orthogonaal is. Kies een orthonormale basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$ van $\mathcal{R}(P)$ en $\{\mathbf{q}_{n+1}, \ldots, \mathbf{q}_m\}$ van $\mathcal{N}(P)$.

**Stap 1:** Omdat $\mathcal{R}(P) \perp \mathcal{N}(P)$, is $Q = [\mathbf{q}_1 \cdots \mathbf{q}_m]$ een unitaire matrix.

**Stap 2:** Bereken $Q^* P Q$. Elke $\mathbf{q}_i$ met $i \leq n$ zit in $\mathcal{R}(P)$, dus $P\mathbf{q}_i = \mathbf{q}_i$. Elke $\mathbf{q}_j$ met $j > n$ zit in $\mathcal{N}(P)$, dus $P\mathbf{q}_j = 0$.

$$Q^*PQ = \text{diag}(1, \ldots, 1, 0, \ldots, 0).$$

**Stap 3:** Deze diagonaalmatrix is reëel en symmetrisch, dus $Q^*PQ = (Q^*PQ)^* = Q^*P^*Q$, waaruit volgt $P = P^*$.

**($\Leftarrow$) De voorwaarde is voldoende:**
Neem willekeurig $\mathbf{x} = P\mathbf{u} \in \mathcal{R}(P)$ en $\mathbf{y} \in \mathcal{N}(P)$.

$$\mathbf{x}^*\mathbf{y} = (P\mathbf{u})^*\mathbf{y} = \mathbf{u}^* P^* \mathbf{y} \overset{P^*=P}{=} \mathbf{u}^* P\mathbf{y} = \mathbf{u}^* \cdot 0 = 0.$$

Dus $\mathcal{R}(P) \perp \mathcal{N}(P)$. $\square$

**Dit zegt ons:** Een projector is orthogonaal precies wanneer de matrix zelftoegevoegd is. Dit is de wiskundige vertaling van "de projectie staat loodrecht op de geprojecteerde ruimte."

### Constructieformule voor een orthogonale projector

Gegeven een deelruimte $\mathcal{D}$ met kolommen van $A$ als basis:

$$P_\mathcal{D} = A(A^*A)^{-1}A^*. \tag{2.13}$$

**Speciale gevallen:**

- Als $\mathcal{D} = \langle \mathbf{a} \rangle$ (één vector): $P_{\langle \mathbf{a} \rangle} = \dfrac{\mathbf{a}\mathbf{a}^*}{\|\mathbf{a}\|_2^2}.$

- Als de basis **orthonormaal** is ($Q = [\mathbf{q}_1 \cdots \mathbf{q}_n]$ met $Q^*Q = I_n$):

$$P_\mathcal{D} = QQ^*, \qquad P_\mathcal{D}\mathbf{v} = \sum_{i=1}^n (\mathbf{q}_i^*\mathbf{v})\,\mathbf{q}_i. \tag{2.15-2.16}$$

**Intuïtie van (2.16):** De projectie is gewoon de som van de componenten van $\mathbf{v}$ in elke basisrichting — dat is precies wat je intuïtief verwacht van een loodrechte projectie.

---

## 2.3 Orthogonalisatieprocedures

### 2.3.1 Link met de QR-factorisatie

We vertrekken van een basis $\{\mathbf{a}_1, \ldots, \mathbf{a}_n\}$ en willen een **geneste** orthonormale basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$ die voldoet aan:

$$\langle \mathbf{a}_1, \ldots, \mathbf{a}_k \rangle = \langle \mathbf{q}_1, \ldots, \mathbf{q}_k \rangle, \quad k = 1, \ldots, n. \tag{2.17}$$

Dit vertaalt zich in een **bovendriehoekige** relatie:

$$\mathbf{a}_j = r_{1j}\mathbf{q}_1 + r_{2j}\mathbf{q}_2 + \cdots + r_{jj}\mathbf{q}_j,$$

wat in matrixvorm is:

$$A = \hat{Q}\hat{R}, \tag{2.18}$$

met $\hat{Q} \in \mathbb{C}^{m \times n}$ (orthonormale kolommen: $\hat{Q}^*\hat{Q} = I_n$) en $\hat{R} \in \mathbb{C}^{n \times n}$ bovendriehoeks. Dit is de **onvolledige QR-factorisatie**.

Door $\hat{Q}$ uit te breiden met $m - n$ extra orthonormale kolommen en $\hat{R}$ met $m - n$ nulrijen, krijgen we de **volledige QR-factorisatie**:

$$A = QR, \quad Q \in \mathbb{C}^{m \times m} \text{ unitair}, \quad R \in \mathbb{C}^{m \times n} \text{ pseudo-driehoeks}. \tag{2.21}$$

> ⚠️ **Belangrijk voor examen:** Ken het verschil tussen volledige en onvolledige QR-factorisatie. Weet dat de onvolledige factorisatie uniek is (als de diagonaalelementen van $\hat{R}$ reëel positief zijn).

### 2.3.2 Gram-Schmidt orthogonalisatie

**Idee:** Verwijder uit elke nieuwe vector $\mathbf{a}_j$ de componenten in de richting van de eerder gevonden orthonormale vectoren $\mathbf{q}_1, \ldots, \mathbf{q}_{j-1}$.

**Afleiding stap $j$:**

$$r_{ij} = \mathbf{q}_i^* \mathbf{a}_j \quad (i < j), \qquad \mathbf{v}_j = \mathbf{a}_j - \sum_{i=1}^{j-1} r_{ij}\mathbf{q}_i = \mathbf{a}_j - P_{\langle \mathbf{q}_1, \ldots, \mathbf{q}_{j-1}\rangle}\mathbf{a}_j,$$

$$r_{jj} = \|\mathbf{v}_j\|_2, \qquad \mathbf{q}_j = \mathbf{v}_j / r_{jj}.$$

**Klassiek Gram-Schmidt algoritme (KGS):**
```
for j = 1 to n:
    v_j = a_j
    for i = 1 to j-1:
        r_ij = q_i* a_j      ← inwendig product met originele a_j
        v_j  = v_j - r_ij * q_i
    r_jj = ||v_j||_2
    q_j  = v_j / r_jj
```

**Gewijzigd Gram-Schmidt algoritme (GGS):**
```
for j = 1 to n:
    v_j = a_j
    for i = 1 to j-1:
        r_ij = q_i* v_j      ← inwendig product met bijgewerkte v_j
        v_j  = v_j - r_ij * q_i
    r_jj = ||v_j||_2
    q_j  = v_j / r_jj
```

**Verschil:** KGS verwijdert alle eerdere componenten in één keer uit $\mathbf{a}_j$; GGS doet dit stapsgewijs. Wiskundig equivalent, maar GGS is numeriek gunstiger omdat gevaarlijke aftrekkingen minder sterk doorwerken.

**Geometrische interpretatie:** In stap $j$ projecteren we $\mathbf{a}_j$ op het orthogonaal complement van $\langle \mathbf{q}_1, \ldots, \mathbf{q}_{j-1}\rangle$:

$$\mathbf{v}_j = (I - P_{\langle \mathbf{q}_1,\ldots,\mathbf{q}_{j-1}\rangle})\mathbf{a}_j. \tag{2.22}$$

**Numerieke stabiliteit:** Gram-Schmidt is niet achterwaarts stabiel. Hoe slechter de originele basis geconditioneerd is (hoe schuiner), hoe sneller de orthogonaliteit van de berekende $\hat{Q}$ verloren gaat. Oplossing: **herorthogonalisatie** — de procedure wordt een tweede maal uitgevoerd op de bekomen $\hat{Q}$.

**Rekenkost:** $\mathcal{O}(2mn^2)$ bewerkingen (zowel KGS als GGS).

**Toepassingsdomein:** Bijzonder geschikt voor grote schalige problemen ($m \gg n$) en voor situaties waarbij de kolommen van $A$ één voor één beschikbaar komen (bijv. Krylov-methodes in Hoofdstuk 9).

### 2.3.3 QR-factorisatie met Givens-rotaties

**Idee:** Maak nullen onder de diagonaal van $A$ door vooraan te vermenigvuldigen met elementaire orthogonale rotaties.

**Givens-rotatie in $\mathbb{R}^2$:** Zoek een orthogonale matrix $G \in \mathbb{R}^{2 \times 2}$ die de tweede component van $\mathbf{x} = [x_1\ x_2]^T$ op nul brengt:

$$G^T \mathbf{x} = \begin{pmatrix}\sqrt{x_1^2 + x_2^2} \\ 0\end{pmatrix}, \quad G = \begin{pmatrix}c & -s \\ s & c\end{pmatrix}, \quad c = \frac{x_1}{\sqrt{x_1^2+x_2^2}},\ s = \frac{x_2}{\sqrt{x_1^2+x_2^2}}. \tag{2.26}$$

$c$ en $s$ zijn de cosinus en sinus van de hoek die $\mathbf{x}$ maakt met de positieve $x_1$-as. $G^T$ roteert $\mathbf{x}$ in wijzerzin naar de $x_1$-as.

**Gegeneraliseerde Givens-rotatie $G_{i,j} \in \mathbb{C}^{m \times m}$:** Identiteitsmatrix behalve op de $(i,i)$-, $(i,j)$-, $(j,i)$-, $(j,j)$-posities, waar $c^*, -s^*, s, c$ staan (met $|c|^2 + |s|^2 = 1$). Toepassen van $G_{i,j}$ op een vector maakt de $j$-de component nul terwijl de andere componenten ongewijzigd blijven (behalve de $i$-de).

**QR-factorisatie:** Maak element per element nullen onder de diagonaal door opeenvolgende Givens-rotaties te combineren. Voorbeeld voor een $4 \times 3$ matrix (nullen worden gecreëerd in de volgorde $G_{1,4}^{(3)}, G_{1,3}^{(2)}, G_{1,2}^{(1)}, \ldots$):

$$G_{3,4}^{(3)} \cdots G_{1,2}^{(1)} A = R \implies A = QR$$

**Eigenschappen van Givens-QR:**
- Berekent inherent een **volledige** QR-factorisatie
- **Achterwaarts stabiel**: $\kappa(G_{i,j}) = 1$, de orthogonaliteit van $Q$ is gegarandeerd
- Rekenkost (enkel $R$): $\mathcal{O}(3mn^2 - n^3)$; inclusief $Q$: $\mathcal{O}(6m^2n - n^3)$ — duurder dan Householder
- Voordeel: nuttig als $A$ al bijna driehoekig is (slechts een paar nullen nodig)

### 2.3.4 QR-factorisatie met Householder-transformaties

**Idee:** In plaats van één nul tegelijk (Givens) kunnen we alle nullen in een volledige kolom tegelijk creëren door gebruik te maken van een **reflector** (spiegeling).

**Reflector t.o.v. het orthogonaal complement van $\langle \mathbf{v} \rangle$:**

$$R_{\langle \mathbf{v} \rangle} = I - 2\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}. \tag{2.36}$$

Dit is een unitaire matrix die een vector $\mathbf{x}$ spiegelt t.o.v. het hypervlak loodrecht op $\mathbf{v}$.

**Householder-transformatie:** Zoek een reflector die de laatste $m-1$ componenten van een vector $\mathbf{x} \in \mathbb{C}^m$ op nul brengt:

$$F_m \mathbf{x} = \mp \|\mathbf{x}\|_2 \mathbf{e}_1.$$

De reflectierichting is $\mathbf{v} = \mathbf{x} \pm \|\mathbf{x}\|_2 \mathbf{e}_1$. Omwille van numerieke stabiliteit (gevaarlijke aftrekking vermijden) kies je:

$$\mathbf{v} = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2 \mathbf{e}_1, \quad F_m = I - 2\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}. \tag{2.39}$$

> ⚠️ **Belangrijk voor examen:** Weet de formule van de Householder-transformatie en begrijp de keuze van het teken (numerieke stabiliteit).

**Efficiënte toepassing:** Je hoeft $F_m$ nooit expliciet te berekenen! Het matrix-vector product $F_m \mathbf{q}$ is:

$$F_m \mathbf{q} = \mathbf{q} - 2\frac{\mathbf{v}^*\mathbf{q}}{\mathbf{v}^*\mathbf{v}}\mathbf{v}. \tag{2.40}$$

Dit kost slechts $\mathcal{O}(m)$ bewerkingen in plaats van $\mathcal{O}(m^2)$ voor een vol matrix-vector product.

**QR-factorisatie:** Pas kolom per kolom Householder-transformaties toe:

$$Q_{n,m} \cdots Q_{2,m} Q_{1,m} A = R \implies A = Q_{1,m} Q_{2,m} \cdots Q_{n,m} R = QR,$$

waarbij $Q_{j,m} = \text{diag}(I_{j-1}, F_{m-j+1})$.

**Rekenkost:** Enkel $R$: $\sim 2mn^2 - \frac{2}{3}n^3$; inclusief $Q$: $\sim 4m^2n - \frac{2}{3}n^3$. Significant goedkoper dan Givens voor volle matrices.

**Vergelijking van de drie methodes:**

| Eigenschap | Gram-Schmidt | Givens | Householder |
|---|---|---|---|
| Factorisatie type | Onvolledig $\hat{Q}\hat{R}$ | Volledig $QR$ | Volledig $QR$ |
| Numerieke stabiliteit | Matig (verbeterd met herortho.) | Stabiel | Stabiel |
| Rekenkost ($R$ + $Q$) | $\mathcal{O}(2mn^2)$ | $\mathcal{O}(6m^2n)$ | $\mathcal{O}(4m^2n)$ |
| Kolommen van $A$ op voorhand? | Nee | Ja | Ja |
| Beste voor | $m \gg n$, online/Krylov | Bijna-driehoekige $A$ | Algemene volle matrices |

---

## 2.4 Beste benadering van een vector in een deelruimte

### Probleemstelling

Gegeven een deelruimte $\mathcal{D} = \langle \mathbf{a}_1, \ldots, \mathbf{a}_n \rangle \subset \mathbb{C}^m$ (met $m > n$) en een vector $\mathbf{b} \in \mathbb{C}^m$. Stel $A = [\mathbf{a}_1 \cdots \mathbf{a}_n]$.

**Kleinste-kwadratenprobleem:** Vind $\hat{\mathbf{x}} \in \mathbb{C}^n$ zodat

$$\min_{\mathbf{x} \in \mathbb{C}^n} \|A\mathbf{x} - \mathbf{b}\|_2. \tag{2.42}$$

De vector $\hat{\mathbf{y}} = A\hat{\mathbf{x}}$ is dan de beste benadering van $\mathbf{b}$ in $\mathcal{D}$.

**Typische toepassing (curve fitting):** Gegeven meetpunten $(x_i, y_i)$ en basisfuncties $\{g_k\}$, zoek coëfficiënten $c_k$ zodat $f(x) = \sum_k c_k g_k(x)$ de data het best benadert. In matrixvorm: $A_{ik} = g_k(x_i)$, $b_i = y_i$, en het kleinste-kwadratenprobleem (2.42) geeft de beste fit.

### 2.4.1 Karakterisering: het normaalstelsel

> ⚠️ **Belangrijk voor examen:** Dit is de centrale stelling van het hoofdstuk. Ken de stellingen, de bewijzen, en de formules.

### Stelling 2.4.1 — Orthogonaliteitsvoorwaarde

**$\hat{\mathbf{y}}$ is de beste benadering van $\mathbf{b}$ in $\mathcal{D}$ als en alleen als $\mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}$.**

### Bewijs

We willen aantonen: als $\mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}$, dan geldt $\|\mathbf{b} - \mathbf{y}\|_2 \geq \|\mathbf{b} - \hat{\mathbf{y}}\|_2$ voor elke $\mathbf{y} \in \mathcal{D}$.

**Stap 1:** Schrijf $\mathbf{b} - \mathbf{y} = (\mathbf{b} - \hat{\mathbf{y}}) + (\hat{\mathbf{y}} - \mathbf{y})$.
→ $\mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}$ en $\hat{\mathbf{y}} - \mathbf{y} \in \mathcal{D}$, dus de twee termen zijn orthogonaal.

**Stap 2:** Pas de stelling van Pythagoras toe:

$$\|\mathbf{b} - \mathbf{y}\|_2^2 = \|\mathbf{b} - \hat{\mathbf{y}}\|_2^2 + \|\hat{\mathbf{y}} - \mathbf{y}\|_2^2 \geq \|\mathbf{b} - \hat{\mathbf{y}}\|_2^2. \tag{2.44}$$

**Conclusie:** De ongelijkheid geldt voor elke $\mathbf{y} \in \mathcal{D}$, dus $\hat{\mathbf{y}}$ is de beste benadering. $\square$

**Dit zegt ons:** De beste benadering is precies die vector in $\mathcal{D}$ zodat de restfout $\mathbf{b} - \hat{\mathbf{y}}$ loodrecht staat op $\mathcal{D}$. Dit is de geometrisch intuïtieve definitie van een orthogonale projectie.

### Stelling 2.4.2 — Orthogonale projectiestelling

**De beste benadering van $\mathbf{b}$ in $\mathcal{D}$ bestaat, is uniek, en wordt gegeven door $\hat{\mathbf{y}} = P_\mathcal{D}\mathbf{b}$.**

### Bewijs

**Stap 1:** Ontbind $\mathbf{b} = P_\mathcal{D}\mathbf{b} + (I - P_\mathcal{D})\mathbf{b}$.
→ $P_\mathcal{D}\mathbf{b} \in \mathcal{R}(P_\mathcal{D}) = \mathcal{D}$, en $(I-P_\mathcal{D})\mathbf{b} \in \mathcal{N}(P_\mathcal{D}) \perp \mathcal{D}$ (want $P_\mathcal{D}$ is orthogonaal).

**Stap 2:** Dus $\mathbf{b} - P_\mathcal{D}\mathbf{b} \perp \mathcal{D}$, en per Stelling 2.4.1 is $\hat{\mathbf{y}} = P_\mathcal{D}\mathbf{b}$ de beste benadering.

**Stap 3 (uniciteit):** Uit de ongelijkheid (2.44) volgt dat elke andere $\mathbf{y} \neq \hat{\mathbf{y}}$ een strikt grotere afstand heeft. $\square$

### Stelling 2.4.3 — Normaalstelsel

**$\hat{\mathbf{x}}$ is een oplossing van (2.42) als en alleen als**

$$A^*A\,\hat{\mathbf{x}} = A^*\mathbf{b}. \tag{2.46}$$

### Bewijs

**($\Rightarrow$) Noodzakelijk:**

Veronderstel $A\hat{\mathbf{x}} = \hat{\mathbf{y}} = P_\mathcal{D}\mathbf{b}$. Dan is $\mathbf{b} - A\hat{\mathbf{x}} \perp \mathcal{D}$, dus:

$$\mathbf{a}_i^*(\mathbf{b} - A\hat{\mathbf{x}}) = 0, \quad i = 1, \ldots, n.$$

In matrixvorm: $A^*(\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0}$, wat herschreven wordt als (2.46). ✓

**($\Leftarrow$) Voldoende:**

Als (2.46) geldt, dan $A^*(A\hat{\mathbf{x}} - \mathbf{b}) = \mathbf{0}$, dus $(A\hat{\mathbf{x}} - \mathbf{b}) \perp \mathcal{D}$. Per Stelling 2.4.1 is $A\hat{\mathbf{x}} = \hat{\mathbf{y}}$. $\square$

**Expliciete oplossing** (als $A$ van volle kolomrang is, d.w.z. $A^*A$ inverteerbaar):

$$\hat{\mathbf{x}} = (A^*A)^{-1}A^*\mathbf{b} =: A^+\mathbf{b}. \tag{2.47}$$

De matrix $A^+ = (A^*A)^{-1}A^*$ heet de **Moore-Penrose pseudo-inverse** van $A$.

**Bij orthogonale basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$** wordt het normaalstelsel diagonaal:

$$\hat{x}_j = \frac{\mathbf{a}_j^*\mathbf{b}}{\|\mathbf{a}_j\|_2^2}, \quad \hat{\mathbf{y}} = \sum_{j=1}^n \frac{\mathbf{a}_j^*\mathbf{b}}{\|\mathbf{a}_j\|_2^2}\,\mathbf{a}_j. \tag{2.48}$$

**Belangrijk voordeel van orthogonale bases:** Als je de deelruimte uitbreidt met een extra basisvector $\mathbf{a}_{n+1}$, veranderen de reeds berekende coëfficiënten $\hat{x}_1, \ldots, \hat{x}_n$ niet. Bij een scheve basis moet je het hele normaalstelsel opnieuw oplossen.

### 2.4.2 Oplossing via de QR-factorisatie

> ⚠️ **Belangrijk voor examen:** Dit is de voorkeursaanpak in de praktijk — beter geconditioneerd dan het normaalstelsel direct oplossen.

**Waarom niet gewoon het normaalstelsel oplossen?** Omdat $\kappa(A^*A) = \kappa(A)^2$. Als $A$ al een conditiegetal van $100$ heeft, is $A^*A$ een factor $10000$ slechter geconditioneerd.

**Aanpak via volledige QR-factorisatie $A = QR$:**

Schrijf $R = \begin{pmatrix}\hat{R} \\ 0\end{pmatrix}$ met $\hat{R} \in \mathbb{C}^{n \times n}$ bovendriehoeks, en partitioneer:

$$Q^*\mathbf{b} = \begin{pmatrix}\mathbf{b}_1 \\ \mathbf{b}_2\end{pmatrix}, \quad \mathbf{b}_1 \in \mathbb{C}^n,\ \mathbf{b}_2 \in \mathbb{C}^{m-n}.$$

Omdat unitaire transformaties de 2-norm bewaren:

$$\|A\mathbf{x} - \mathbf{b}\|_2^2 = \|Q^*(A\mathbf{x} - \mathbf{b})\|_2^2 = \left\|\begin{pmatrix}\mathbf{b}_1 - \hat{R}\mathbf{x} \\ \mathbf{b}_2\end{pmatrix}\right\|_2^2 = \|\mathbf{b}_1 - \hat{R}\mathbf{x}\|_2^2 + \|\mathbf{b}_2\|_2^2. \tag{2.49-2.50}$$

Alleen de eerste term kan geminimaliseerd worden. De kleinste-kwadratenoplossing $\hat{\mathbf{x}}$ voldoet dus aan het **bovendriehoekige stelsel**:

$$\hat{R}\,\hat{\mathbf{x}} = \mathbf{b}_1. \tag{2.51}$$

De minimale fout is $\min_\mathbf{x} \|A\mathbf{x} - \mathbf{b}\|_2 = \|\mathbf{b}_2\|_2$.

**Algoritme:**
1. Bereken de (onvolledige) QR-factorisatie: $A = \hat{Q}\hat{R}$
2. Bereken $\mathbf{b}_1 = \hat{Q}^*\mathbf{b}$ (de eerste $n$ componenten van $Q^*\mathbf{b}$)
3. Los het driehoekige stelsel $\hat{R}\hat{\mathbf{x}} = \mathbf{b}_1$ op via achterwaartse substitutie

**Interpretatie:** Stap 1 vervangt de basis $\{\mathbf{a}_k\}$ door een orthonormale basis $\{\mathbf{q}_k\}$. Stap 2 berekent de coëfficiënten in die orthonormale basis. Stap 3 transformeert terug naar de originele basis.

### 2.4.3 Conditie en stabiliteit

**Conditiegetallen van het kleinste-kwadratenprobleem:**

Definieer:
- $\theta = \arccos\left(\frac{\|\hat{\mathbf{y}}\|_2}{\|\mathbf{b}\|_2}\right)$: de hoek tussen $\mathbf{b}$ en de projectie $\hat{\mathbf{y}}$ (hoek met $\mathcal{D}$)
- $\eta = \frac{\|A\|_2 \|\hat{\mathbf{x}}\|_2}{\|A\hat{\mathbf{x}}\|_2}$: relatieve grootte van de oplossing

Dan gelden de conditiegetallen:

$$\kappa_b^{\hat{\mathbf{y}}} = \frac{1}{\cos\theta}, \qquad \kappa_b^{\hat{\mathbf{x}}} = \frac{\kappa(A)}{\eta\cos\theta},$$
$$\kappa_A^{\hat{\mathbf{y}}} \leq \frac{\kappa(A)}{\cos\theta}, \qquad \kappa_A^{\hat{\mathbf{x}}} \leq \frac{\kappa(A)}{\eta} + \frac{\kappa(A)^2 \tan\theta}{\eta}. \tag{2.52}$$

**Twee slechte situaties:**
- $\theta \approx \pi/2$ (d.w.z. $\mathbf{b}$ staat bijna loodrecht op $\mathcal{D}$, dus $\hat{\mathbf{y}} \approx 0$): de projectie is bijna nul en wordt enorm versterkt → slechte conditie.
- $\kappa(A)$ groot: directe voortplanting van fouten.

**Waarom normaalstelsel instabiel is:**

Als je het normaalstelsel $A^*A\hat{\mathbf{x}} = A^*\mathbf{b}$ direct oplost, verwacht je een fout:

$$\frac{\|\tilde{\mathbf{x}} - \hat{\mathbf{x}}\|_2}{\|\hat{\mathbf{x}}\|_2} = \mathcal{O}\big(\kappa(A^*A) \cdot \epsilon_{\text{mach}}\big) = \mathcal{O}\big(\kappa(A)^2 \cdot \epsilon_{\text{mach}}\big). \tag{2.53}$$

De QR-aanpak geeft slechts $\mathcal{O}(\kappa(A) \cdot \epsilon_{\text{mach}})$ — veel beter!

**Uitzondering:** Als de fit uitstekend is ($\theta$ klein), is de factor $\kappa(A)^2$ in (2.52) niet dominant, en kan het normaalstelsel acceptabel zijn.

> ⚠️ **Belangrijk voor examen:** Weet waarom de QR-aanpak beter geconditioneerd is dan het normaalstelsel ($\kappa(A)$ vs. $\kappa(A)^2$). Ken de conditiegetallen en begrijp de rol van $\theta$ en $\kappa(A)$.

---

## Samenvatting

| Concept | Formule / resultaat |
|---|---|
| Grammatrix | $G = A^*A$ |
| Coëfficiënten bij orthogonale basis | $c_j = \mathbf{a}_j^*\mathbf{v}/\|\mathbf{a}_j\|_2^2$ |
| Orthogonale projector | $P_\mathcal{D} = A(A^*A)^{-1}A^* = QQ^*$ (bij orthonormale basis) |
| Karakterisering beste benadering | $\mathbf{b} - \hat{\mathbf{y}} \perp \mathcal{D}$ |
| Normaalstelsel | $A^*A\hat{\mathbf{x}} = A^*\mathbf{b}$ |
| Beste oplossing (QR-aanpak) | $\hat{R}\hat{\mathbf{x}} = \mathbf{b}_1$ (bovendriehoeks) |
| Conditie normaalstelsel | $\kappa(A^*A) = \kappa(A)^2$ → instabiel |
| Conditie QR-aanpak | $\kappa(A)$ → stabiel |

**Verbinding met de rest van de cursus:** Hoofdstuk 3 veralgemeen deze ideeën naar abstracte functieruimten. Alle begrippen (normaalstelsel, orthogonale projectie, grammatrix) komen terug, maar nu voor oneindigdimensionale ruimten van functies.
