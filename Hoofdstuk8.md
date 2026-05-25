# Deel 8 – Grafen en eigenwaarden in datawetenschappen

## Overzicht

Dit hoofdstuk toont hoe **eigenwaarden en eigenvectoren** opdagen in drie concrete netwerkproblemen:

| Toepassing | Eigenwaardenkern |
|---|---|
| **PageRank** | Dominante eigenvector (eigenwaarde 1) van de linkmatrix |
| **Meest centrale knoop** | Diagonaal van de matrixexponentiaal $e^A$ via eigenwaardenontbinding |
| **Graafpartitionering** | Tweede-kleinste eigenvector (Fiedler-vector) van de graaf-Laplaciaan |

De rode draad: graafstructuren leiden op een natuurlijke manier tot matrices, en de antwoorden op praktische vragen over het netwerk zitten verscholen in de eigenwaarden en eigenvectoren van die matrices.

---

## 8.1 PageRank

### Het idee

Google's originele algoritme rankt webpagina's op basis van de **structuur van het hyperlinknetwerk**, niet alleen op de aanwezigheid van zoektermen. De redenering:
- Een pagina is belangrijk als *veel andere* pagina's naar haar verwijzen.
- Een pagina is belangrijker als de *verwijzende pagina's zelf* belangrijk zijn.
- Een pagina met *weinig links* verdeelt haar "belang" minder: haar verwijzing telt zwaarder.

### De linkmatrix $A$

Beschouw $N$ webpagina's. De **linkmatrix** $A \in \mathbb{R}^{N \times N}$ is:

$$A_{ij} = \begin{cases} \frac{1}{N_i} & \text{als pagina } P_i \text{ een link heeft naar } P_j \\ 0 & \text{anders,} \end{cases}$$

waarbij $N_i$ het aantal links op pagina $P_i$ is. Matrix $A$ is de **gewogen verbindingsmatrix** van de gerichte graaf van hyperlinks.

**Markov-interpretatie:** $A_{ij}$ is de kans dat een gebruiker die willekeurig op links klikt, van pagina $i$ naar pagina $j$ gaat. Matrix $A$ is de **transitiematrix** van een Markov-keten op het web.

### Van scores naar een eigenwaardenprobleem

We zoeken een scorevector $\boldsymbol{\pi} = [r(P_1)\ \cdots\ r(P_N)]^T$ zodat:

$$r(P_i) \approx \sum_{j \in B_i} \frac{r(P_j)}{N_j}, \quad i = 1, \ldots, N, \tag{8.2}$$

waarbij $B_i$ de verzameling pagina's is die naar $P_i$ verwijzen. In matrixvorm:

$$A^T \boldsymbol{\pi} \approx \boldsymbol{\pi}. \tag{8.3}$$

Dit zegt: de scorevector $\boldsymbol{\pi}$ is (bijna) een **eigenvector van $A^T$ bij eigenwaarde 1**.

### Twee aanpassingen voor wiskundige welgesteldheid

Het systeem (8.3) is mogelijk inconsistent of niet goed gedefinieerd als $A$ bijzondere structuur heeft. Twee aanpassingen verhelpen dit:

**Aanpassing 1 — Nul-rijen vervangen:**  
Pagina's zonder uitgaande links geven matrix $A$ nul-rijen, waardoor de scoreverdeling "wegloopt". Vervang elke nul-rij door $\frac{1}{N}\mathbf{1}^T$ (virtuele links naar alle pagina's). De resulterende matrix $\tilde{A}$ is **rij-stochastisch** (alle rijen sommeren naar 1, alle elementen $\geq 0$).

**Aanpassing 2 — Damping factor $\alpha$:**  
Voeg een convexe combinatie toe met een matrix die overal $\frac{1}{N}$ heeft:

$$\hat{A} = (1 - \alpha)\tilde{A} + \alpha \cdot \frac{1}{N}\mathbf{1}\mathbf{1}^T, \quad \alpha \in (0,1) \text{ klein}. \tag{\text{bv. } \alpha = 0.1}$$

Dit simuleert een gebruiker die af en toe willekeurig naar een pagina springt (in plaats van altijd op een link te klikken). Matrix $\hat{A}$ is nu ook **irreduceerbaar**: de graaf is sterk geconnecteerd (je kan vanuit elke knoop elke andere knoop bereiken).

### Sleutelstelling

### Stelling 8.1.1
**Als $\hat{A}$ een irreduceerbare rij-stochastische matrix is, dan:**
- **is de dominante eigenwaarde enkelvoudig en gelijk aan 1;**
- **kan de bijhorende eigenvector van $\hat{A}^T$ gekozen worden met alle elementen $\geq 0$.**

> ⚠️ **Belangrijk voor examen:** Weet deze stelling en begrijp waarom de twee aanpassingen aan $A$ nodig zijn om ze toe te kunnen passen.

Het PageRank-probleem wordt dan:

$$\hat{A}^T \boldsymbol{\pi} = \boldsymbol{\pi}. \tag{8.4}$$

Dit is een **eigenwaardenprobleem**: zoek de eigenvector bij eigenwaarde 1 (de dominante eigenwaarde) van $\hat{A}^T$. De elementen van $\boldsymbol{\pi}$ hebben allemaal hetzelfde teken en kunnen als scores geïnterpreteerd worden.

**Algoritme:** De methode van de machten (zie Hoofdstuk 9) berekent iteratief de eigenvector bij de dominante eigenwaarde. Elke iteratie is één matrix-vectorproduct met $\hat{A}^T$ — efficiënt voor de enorm grote maar ijle webgraaf.

### Voorbeeld: 6 pagina's

Matrix $A$ en $\hat{A}$ (met $\alpha = 0.1$) zoals beschreven in de cursus. De genormaliseerde dominante eigenvector (PageRank-vector) is:

$$\boldsymbol{\pi} \approx [0.037,\ 0.054,\ 0.042,\ 0.375,\ 0.206,\ 0.286]^T.$$

Pagina 4 heeft de hoogste score, ook al heeft ze misschien niet de meeste links — ze wordt verwezen vanuit zwaargewichtige pagina's.

**Markov-interpretatie:** $\boldsymbol{\pi}$ genormaliseerd met de $\ell^1$-norm is de **stationaire verdeling** van de Markov-keten: de kans om in elke pagina te zijn na oneindig lang willekeurig surfen.

---

## 8.2 Meest centrale knoop

### Probleemstelling

Gegeven een **ongerichte enkelvoudige graaf** (geen zelfverbindingen, maximaal één boog tussen twee knopen) met $N$ knopen en **binaire verbindingsmatrix** $A$ (symmetrisch, $A_{ij} \in \{0,1\}$). Welke knoop is het meest "centraal" in het netwerk?

### Van boogaantal naar lussen van alle lengtes

Een eerste maat voor centraliteit: het aantal bogen vanuit knoop $i$ = $\sum_j A_{ij}$. Merk op dat voor binaire, symmetrische $A$:

$$\sum_j A_{ij} = \sum_j A_{ij}A_{ji} = (A^2)_{ii}. \tag{8.5}$$

Dit is tegelijk het aantal **lussen van lengte 2** startend in knoop $i$ (heen en terug over een boog).

**Veralgemeening (Stelling 8.2.2):** Het aantal lussen van lengte $n$ startend in knoop $i$ is $(A^n)_{ii}$.

**Intuïtie voor grotere lengtes:** Een knoop verbonden met andere goed-verbonden knopen is centraler dan een knoop met evenveel directe buren maar die allemaal geïsoleerd zijn. Lussen van lengte 3, 4, ... meten dit indirect.

### De matrixexponentiaal als centraliteitsmaat

We combineren lussen van alle lengtes, met afnemend gewicht voor grotere lengtes (want verre verbindingen tellen minder):

$$\sum_{n=2}^\infty \frac{1}{n!}(A^n)_{ii} \quad \rightarrow \quad \left(e^A\right)_{ii} - 1 - (A)_{ii} = \left(e^A\right)_{ii} - 1,$$

want de diagonaal van $A$ is nul. Voor de rangschikking maakt de constante $-1$ niet uit, dus:

### Definitie 8.2.1 — Centraliteit
$$\text{centraliteit van knoop } i = \left(e^A\right)_{ii}.$$

De **meest centrale knoop** is de knoop waarvoor $(e^A)_{ii}$ maximaal is.

### Berekening via eigenwaardenontbinding

Als $A = X \Lambda X^{-1}$ (eigenwaardenontbinding), dan:

$$e^A = X\, e^\Lambda\, X^{-1} = X \begin{pmatrix} e^{\lambda_1} & & \\ & \ddots & \\ & & e^{\lambda_N} \end{pmatrix} X^{-1}.$$

Omdat $A$ symmetrisch is, is $X$ orthogonaal ($X^{-1} = X^T$) en heeft $A$ reële eigenwaarden. In de praktijk berekent men de volledige eigenwaardenontbinding van $A$ (via het QR-algoritme, zie Hoofdstuk 9) en assembleert men $e^A$ via bovenstaande formule.

> ⚠️ **Belangrijk voor examen:** Ken de definitie van centraliteit via $(e^A)_{ii}$ en weet hoe de matrixexponentiaal berekend wordt via de eigenwaardenontbinding.

### Voorbeeld

Voor de graaf uit Voorbeeld 8.3:
- $(A^2)_{ii} = [2,3,3,3,4,1,1,1]^T$ → knoop 5 heeft de meeste directe buren
- $(e^A)_{ii} \approx [3.03, 4.34, 3.99, 4.34, 3.91, 1.69, 1.69, 1.69]^T$ → knopen 2 en 4 zijn meest centraal

**Verklaring:** Knoop 5 heeft 4 directe buren maar die zijn minder goed verbonden. Knopen 2 en 4 hebben 3 directe buren maar zijn verbonden met andere goed-verbonden knopen. De lussen van grotere lengte geven de doorslag.

---

## 8.3 Graafpartitionering

### Probleemstelling

Gegeven een geconnecteerde enkelvoudige ongerichte graaf $G$ met $N$ knopen ($N$ even). Verdeel de knopen in twee gelijke groepen $K_1$ en $K_2$ zodat het **aantal bogen dat de grens oversteekt** minimaal is.

**Toepassing:** netwerken opdelen in clusters, parallel computing (taakverdeeling over processoren met minimale communicatie), afbeeldingssegmentatie.

### De graaf-Laplaciaan

Definieer de **graad** van knoop $i$ als $g(i) = \sum_j A_{ij}$ (aantal bogen uit $i$). De **graaf-Laplaciaan** is:

$$L = D - A, \quad D = \text{diag}(g(1), \ldots, g(N)),$$

ofwel element per element:

$$L_{ij} = \begin{cases} g(i) & \text{als } i = j \\ -1 & \text{als } i \neq j \text{ en knopen } i,j \text{ verbonden} \\ 0 & \text{anders.} \end{cases}$$

### Stelling 8.3.3 — Eigenschappen van de graaf-Laplaciaan

> ⚠️ **Belangrijk voor examen:** Ken deze drie eigenschappen en hun bewijsideeën.

1. **$L$ is symmetrisch positief semi-definiet.**
2. **$L$ heeft altijd eigenwaarde 0, met eigenvector $\mathbf{1} = [1\ \cdots\ 1]^T$.**  
   (Want $L\mathbf{1} = D\mathbf{1} - A\mathbf{1} = \mathbf{g} - \mathbf{g} = \mathbf{0}$.)
3. **Het aantal verbonden componenten van $G$ = de multipliciteit van eigenwaarde 0.**

Voor een geconnecteerde graaf is de eigenwaarde 0 enkelvoudig, en de eigenwaarden voldoen aan:

$$0 = \lambda_N < \lambda_{N-1} \leq \lambda_{N-2} \leq \cdots \leq \lambda_1.$$

De tweede-kleinste eigenwaarde $\lambda_{N-1}$ heet de **algebraïsche connectiviteit** van $G$. Hoe kleiner $\lambda_{N-1}$, hoe gemakkelijker de graaf te partitioneren is.

### Van partitieprobleem naar eigenwaardenprobleem

Stel de partitie voor door $\mathbf{x} \in \{-1, +1\}^N$:
$$x_i = \begin{cases} +1 & i \in K_1 \\ -1 & i \in K_2. \end{cases}$$

De eis $\#K_1 = \#K_2$ wordt: $\mathbf{1}^T\mathbf{x} = 0$.

### Sleutelgelijkheid

$$\frac{1}{4}\mathbf{x}^T L \mathbf{x} = \#B_{\mathbf{x}},$$

het aantal bogen dat de grens oversteekt. Dit volgt door directe uitwerking:

$$\mathbf{x}^T L \mathbf{x} = \sum_{i,j} L_{ij} x_i x_j = \sum_{(i,j)\in B} (-x_i x_j + x_i^2) = \sum_{(i,j)\in B} (-x_i x_j + 1).$$

Als $i$ en $j$ in verschillende partities zitten: $x_i x_j = -1$, bijdrage $= 2$. Als ze in dezelfde partitie zitten: $x_i x_j = +1$, bijdrage $= 0$.

Het **discrete optimalisatieprobleem** is dan:

$$\min_{\mathbf{x} \in \{-1,+1\}^N}\ \frac{1}{4}\mathbf{x}^T L \mathbf{x}, \quad \mathbf{1}^T \mathbf{x} = 0. \tag{8.7}$$

**Probleem:** het aantal mogelijke partities is $\frac{1}{2}\binom{N}{N/2}$ — exponentieel groot in $N$. Voor $N = 64$ zijn dat al $9 \times 10^{17}$ mogelijkheden. Brute-force is onmogelijk.

### Relaxatie naar een continue oplossing

Vervang de binaire beperking $\mathbf{x} \in \{-1,+1\}^N$ door de continue beperking $\|\mathbf{x}\|_2^2 = N$:

$$\min_{\mathbf{x} \in \mathbb{R}^N}\ \frac{1}{4}\mathbf{x}^T L \mathbf{x}, \quad \|\mathbf{x}\|_2^2 = N,\ \mathbf{1}^T\mathbf{x} = 0. \tag{8.8}$$

### Stelling 8.3.4 — Oplossing van de relaxatie

**Als $\lambda_{N-1} < \lambda_{N-2}$, dan is het minimum van (8.8) gelijk aan $\frac{N}{4}\lambda_{N-1}$, bereikt door $\mathbf{x} = \sqrt{N}\,\mathbf{v}_{N-1}$, de eigenvector bij $\lambda_{N-1}$.**

### Bewijs

**Stap 1:** Ontbind $L = V\,\text{diag}(\lambda_1, \ldots, \lambda_N)\,V^T$ met $V = [\mathbf{v}_1 \cdots \mathbf{v}_N]$ orthogonaal.

**Stap 2:** Merk op dat $\mathbf{v}_N = \frac{1}{\sqrt{N}}\mathbf{1}$ (de eigenvector bij $\lambda_N = 0$). De eis $\mathbf{1}^T\mathbf{x} = 0$ wordt daarmee $\mathbf{v}_N^T\mathbf{x} = 0$.

**Stap 3:** Substitueer $\mathbf{x} = V\mathbf{y}$. Omdat $V$ orthogonaal is:
- $\|\mathbf{x}\|_2^2 = \|\mathbf{y}\|_2^2 = N$
- $\mathbf{v}_N^T\mathbf{x} = \mathbf{v}_N^T V\mathbf{y} = y_N = 0$
- $\mathbf{x}^T L \mathbf{x} = \mathbf{y}^T \text{diag}(\lambda_1, \ldots, \lambda_N)\mathbf{y} = \sum_{i=1}^N \lambda_i y_i^2$

**Stap 4:** Het probleem wordt:

$$\min_{\mathbf{y} \in \mathbb{R}^N}\ \frac{1}{4}\sum_{i=1}^N \lambda_i y_i^2, \quad \|\mathbf{y}\|_2^2 = N,\ y_N = 0.$$

→ $y_N = 0$ is al opgelegd. De overige $y_i$ moeten een norm-eis voldoen. De doelfunctie wordt geminimaliseerd door **alle gewicht te leggen bij de kleinste $\lambda_i$** met $i < N$, d.i. $\lambda_{N-1}$.

**Stap 5:** Optimale $\mathbf{y} = [0\ \cdots\ 0\ \sqrt{N}\ 0]^T$ (component $N-1$ is $\sqrt{N}$, de rest nul). Minimum $= \frac{1}{4}\lambda_{N-1} \cdot N$.

Terug naar $\mathbf{x}$: $\mathbf{x} = V\mathbf{y} = \sqrt{N}\,\mathbf{v}_{N-1}$. $\square$

**Dit zegt ons:** De optimale continue partitie is volledig bepaald door de **Fiedler-vector** $\mathbf{v}_{N-1}$ — de eigenvector horende bij de tweede-kleinste eigenwaarde van $L$.

### Het algoritme

> ⚠️ **Belangrijk voor examen:** Ken het algoritme, de naam "Fiedler-vector" en de ondergrens voor het aantal te doorbreken bogen.

1. Bereken de **Fiedler-vector** $\mathbf{v}_{N-1}$: de eigenvector van $L$ bij $\lambda_{N-1}$.
2. Partitioneer op basis van het **teken** van de elementen van $\mathbf{v}_{N-1}$:
   $$K_1 = \{i : (\mathbf{v}_{N-1})_i \geq 0\}, \quad K_2 = \{i : (\mathbf{v}_{N-1})_i < 0\}.$$

**Ondergrens:** het minimaal aantal te doorbreken bogen is $\geq \frac{N}{4}\lambda_{N-1}$.

**Numerieke methode:** voor grote, ijle $L$ gebruikt men het **Lanczos-algoritme** (§9.1.3.5) om enkel de tweede-kleinste eigenwaarde en bijhorende eigenvector te berekenen — zonder de volledige eigenwaardenontbinding.

### Voorbeeld: $N = 10$ knopen

Kleinste eigenwaarden van $L$: $\lambda_{10} = 0$, $\lambda_9 = 0.451$, $\lambda_8 = 1.229$.

Fiedler-vector $\mathbf{v}_9 \approx [+0.36, +0.20, -0.21, -0.29, -0.47, -0.44, -0.06, +0.19, +0.36, +0.36]^T$.

Partitie op teken: $K_1 = \{1, 2, 8, 9, 10\}$, $K_2 = \{3, 4, 5, 6, 7\}$.

Ondergrens doorbroken bogen: $\lceil \frac{10}{4} \times 0.451 \rceil = \lceil 1.128 \rceil = 2$. De gevonden partitie doorbreekt inderdaad precies 2 bogen — de ondergrens is scherp.

---

## Samenvatting

| Toepassing | Matrix | Eigenwaardekern | Algoritme |
|---|---|---|---|
| **PageRank** | $\hat{A}^T$ (linkmatrix) | Dominante eigenvector ($\lambda = 1$) | Methode van de machten |
| **Meest centrale knoop** | $A$ (verbindingsmatrix) | Diagonaal van $e^A = X e^\Lambda X^{-1}$ | Volledige eigenwaardenontbinding |
| **Graafpartitionering** | $L = D - A$ (Laplaciaan) | Eigenvector bij $\lambda_{N-1}$ (Fiedler) | Lanczos-algoritme |

**Verbinding met de rest van de cursus:**
- De eigenwaardenalgoritmen (methode van de machten, Lanczos, QR-iteratie) die hier nodig zijn, worden uitgewerkt in **Hoofdstuk 9**.
- De structuur van de graaf-Laplaciaan ($L = D - A$, symmetrisch, positief semi-definiet) maakt het probleem rekenkundig haalbaar: men kan gebruik maken van efficiënte methodes voor symmetrische matrices.
