# Deel 12 – IJle representaties en benaderingen

## Overzicht

Dit hoofdstuk behandelt technieken om data (voorgesteld als een matrix) compact voor te stellen door gebruik te maken van onderliggende **ijle structuur** — de idee dat de meeste informatie in relatief weinig coëfficiënten zit.

| Sectie | Techniek | Kernvraag |
|---|---|---|
| 12.1 | Singuliere-waardenontbinding (SVD) | Hoe ontbind ik een willekeurige matrix in orthogonale factoren? |
| 12.2 | Lage-rangbenaderingen | Hoe benader ik een matrix $A$ het best door een matrix van lage rang? |
| 12.3 | Principal component analysis (PCA) | Hoe reduceer ik de dimensie van een dataset? |
| 12.4 | DCT en wavelets | Hoe benut ik frequentie- of resolutiestructuur voor compressie? |

**Rode draad:** Als de singuliere waarden van een matrix snel dalen, dan zit de meeste informatie in de eerste paar singuliere vectoren. Dit kan benut worden voor compressie, dimensiereductie, en efficiënte rekenmethodes.

---

## 12.1 De singuliere-waardenontbinding (SVD)

### Definitie

> ⚠️ **Belangrijk voor examen:** Ken de definitie, de eigenschappen, en het verband met eigenwaarden van $A^*A$.

De **singuliere-waardenontbinding** van $A \in \mathbb{C}^{m \times n}$ is een ontbinding van de vorm:

$$A = U \Sigma V^*, \tag{12.1}$$

met:
- $U \in \mathbb{C}^{m \times m}$ unitair: $U^* U = I_m$ (linker singuliere vectoren als kolommen)
- $V \in \mathbb{C}^{n \times n}$ unitair: $V^* V = I_n$ (rechter singuliere vectoren als kolommen)
- $\Sigma \in \mathbb{R}^{m \times n}$ met $\Sigma = \text{diag}(\sigma_1, \ldots, \sigma_p, 0, \ldots)$, $p = \min(m,n)$

De **singuliere waarden** $\sigma_i \geq 0$ zijn geordend zodat:

$$\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_p \geq 0.$$

De **rang** $r$ van $A$ is het aantal strikt positieve singuliere waarden: $\sigma_1 \geq \cdots \geq \sigma_r > 0 = \sigma_{r+1} = \cdots$.

**Intuïtie:** Elke lineaire afbeelding $A$ kan gezien worden als drie opeenvolgende operaties: roteer/reflecteer ($V^*$), schaal langs de assen ($\Sigma$), roteer/reflecteer opnieuw ($U$).

### Eigenschappen

Uit $A = U \Sigma V^*$ en $AV = U\Sigma$ volgt kolom per kolom:

$$A v_i = \sigma_i u_i. \tag{12.2}$$

Dit geeft de volgende structurele informatie:

| Ruimte | Basis |
|---|---|
| Kolomruimte $\mathcal{R}(A)$ | Eerste $r$ kolommen van $U$: $\{u_1, \ldots, u_r\}$ |
| Nulruimte $\mathcal{N}(A)$ | Laatste $n-r$ kolommen van $V$: $\{v_{r+1}, \ldots, v_n\}$ |
| Rijruimte van $A$ | Eerste $r$ kolommen van $V$: $\{v_1, \ldots, v_r\}$ |

**Rang-1 ontbinding:**

$$A = \sum_{i=1}^r \sigma_i u_i v_i^* = \sigma_1 u_1 v_1^* + \sigma_2 u_2 v_2^* + \cdots + \sigma_r u_r v_r^*.$$

Elke term $\sigma_i u_i v_i^*$ is een matrix van rang 1. Dit is de bouwsteen voor lage-rangbenaderingen.

**Verband met eigenwaarden van $A^*A$:**

Uit $A^*AV = V\Sigma^T \Sigma$ (via $A = U\Sigma V^*$) volgt dat $v_i$ een **eigenvector** is van $A^*A$ met eigenwaarde $\sigma_i^2$. Voor $i > r$: eigenwaarde is nul.

**Normen:**

$$\|A\|_2 = \sigma_1, \qquad \|A\|_F = \sqrt{\sum_{i=1}^r \sigma_i^2}.$$

Als $A \in \mathbb{R}^{m \times n}$, zijn $U$ en $V$ orthogonale matrices (reëel).

### Berekening

De berekening van de SVD verloopt in twee stappen, analoog aan de QR-$\to$-Hessenberg-$\to$-QR aanpak voor eigenwaarden (Sectie 9.2):

**Stap 1 — Bidiagonalisatie:** Vermenigvuldig links en rechts met orthogonale matrices om $A$ te herleiden tot een (bovenste) bi-diagonale matrix:

$$Q_l A Q_r = \begin{pmatrix} A_0 \\ 0 \end{pmatrix}, \quad A_0 = \begin{pmatrix} b_{11} & b_{12} & & \\ & b_{22} & b_{23} & \\ & & \ddots & \ddots \\ & & & b_{nn} \end{pmatrix}.$$

Anders dan bij de Hessenberg-reductie kunnen $Q_l$ en $Q_r$ **onafhankelijk** gekozen worden, waardoor een bi-diagonale (ipv Hessenberg) structuur altijd haalbaar is.

**Stap 2 — SVD van $A_0$:** Bereken de SVD van de bi-diagonale matrix $A_0 = U_0 \Sigma_0 V_0^*$ iteratief. De volledige SVD van $A$ is dan:

$$A = \underbrace{Q_l \begin{pmatrix} U_0 \\ & I \end{pmatrix}}_{U} \Sigma_0 \underbrace{V_0^* Q_r^*}_{V^*}.$$

**Stabiele aanpak:** Ipv eigenwaarden van $A_0^* A_0$ te berekenen (numeriek onstabiel wegens kwadrateren van fouten), berekent men eigenwaarden van de Hermitische matrix:

$$H = \begin{pmatrix} 0 & A_0 \\ A_0^* & 0 \end{pmatrix},$$

waarvan de eigenwaarden precies $\pm\sigma_1, \ldots, \pm\sigma_n$ zijn. Zowel linker als rechter singuliere vectoren volgen uit de bijhorende eigenvectoren.

---

## 12.2 Lage-rangbenaderingen

**Probleemstelling:** Gegeven $A \in \mathbb{C}^{m \times n}$ van rang $r$, zoek de beste benadering $B$ van rang $\leq k$, met $k < r$.

Definieer: $\mathcal{M}_k(m,n) := \{B \in \mathbb{C}^{m \times n} : \text{rang}(B) \leq k\}$.

**Waarom geen standaardtheorie?** $\mathcal{M}_k(m,n)$ is **geen deelruimte** — ze is niet gesloten onder optelling:

$$B_1 = \begin{pmatrix}1 & 0\\ 0 & 0\end{pmatrix}, \quad B_2 = \begin{pmatrix}0 & 0\\ 0 & 1\end{pmatrix} \in \mathcal{M}_1(2,2), \quad \text{maar } B_1 + B_2 = I_2 \notin \mathcal{M}_1(2,2).$$

De theorie van de beste benadering in deelruimten (Hoofdstukken 2–3) is dus niet direct toepasbaar.

**Motivatie voor lage-rangbenaderingen:**
- Een rang-$k$ matrix $B \in \mathbb{C}^{m \times n}$ kan worden geschreven als $B = XY^*$ met $X \in \mathbb{C}^{m \times k}$, $Y \in \mathbb{C}^{n \times k}$: geheugen $O((m+n)k)$ ipv $O(mn)$.
- Veel fysische systemen zijn inherent laag-dimensionaal: kleine $k$ volstaat al voor een goede benadering.

### Methode 1 — QR-factorisatie met kolomverwisselingen

**Idee:** Pas de QR-factorisatie aan door bij elke stap de kolom met grootste norm naar voor te plaatsen. Dit geeft:

$$AP = QR, \qquad \text{of} \qquad A = QRP^T, \tag{12.4}$$

met $Q$ unitair, $P$ een permutatiematrix, en $R$ bovendriehoeks met:

$$|r_{11}| \geq |r_{22}| \geq \cdots \geq |r_{nn}|, \qquad |r_{ii}| \geq |r_{ij}|, \quad i \leq j. \tag{12.5}$$

De eerste rijen van $R$ dragen de meest informatie. De rang-$k$ benadering is:

$$A_k^{qr} = Q_k R_k P^T,$$

met $Q_k$ de eerste $k$ kolommen van $Q$ en $R_k$ de eerste $k$ rijen van $R$.

**Voordeel:** Significant goedkoper te berekenen dan de SVD.  
**Nadeel:** Niet optimaal — er bestaat geen garantie dat dit de beste rang-$k$ benadering is.

### Methode 2 — Getrunceerde SVD

Gebruik de rang-1 ontbinding $A = \sum_{i=1}^r \sigma_i u_i v_i^*$ en knip af na $k$ termen:

$$A_k^{svo} = \sum_{i=1}^k \sigma_i u_i v_i^* = U_k \Sigma_k V_k^*, \tag{12.6}$$

met $U_k = [u_1 \cdots u_k]$, $V_k = [v_1 \cdots v_k]$, $\Sigma_k = \text{diag}(\sigma_1, \ldots, \sigma_k)$.

**Gemaakte fout:**

$$\|A_k^{svo} - A\|_2 = \sigma_{k+1}, \qquad \|A_k^{svo} - A\|_F = \sqrt{\sigma_{k+1}^2 + \cdots + \sigma_r^2}. \tag{12.7}$$

**Compressie:** De drie factoren van $A_k^{svo}$ bevatten $(m + 1 + n) \cdot k$ getallen, ipv $mn$.

### Stelling 12.2.1 — Beste rang-$k$ benadering (Eckart–Young–Mirsky)

> ⚠️ **Belangrijk voor examen:** Dit is het kernresultaat van het hoofdstuk. Ken het bewijs voor de spectraalnorm.

$$\|A_k^{svo} - A\|_2 = \min_{B \in \mathcal{M}_k(m,n)} \|B - A\|_2$$

$$\|A_k^{svo} - A\|_F = \min_{B \in \mathcal{M}_k(m,n)} \|B - A\|_F.$$

**Dit zegt ons:** De getrunceerde SVD levert de **beste** rang-$k$ benadering in zowel de spectraalnorm als de Frobeniusnorm. Geen andere rang-$k$ matrix benadert $A$ beter.

### Bewijs (spectraalnorm)

**Doel:** Tonen dat voor elke $B \in \mathcal{M}_k(m,n)$ geldt $\|B - A\|_2 \geq \sigma_{k+1}$.

**Stap 1:** Schrijf $B = XY^*$ met $X \in \mathbb{C}^{m \times k}$, $Y \in \mathbb{C}^{n \times k}$.

**Stap 2:** Zoek een vector $w \in \mathbb{C}^n$, $\|w\|_2 = 1$, zodat tegelijkertijd:
- $w \in \mathcal{R}([v_1 \cdots v_{k+1}])$ (de ruimte opgespannen door de eerste $k+1$ rechter singuliere vectoren, dimensie $k+1$)
- $w \in \mathcal{N}(Y^*)$ (de nulruimte van $Y^*$, dimensie $\geq n - k$)

Zo'n $w$ bestaat omdat de som van de dimensies $(k+1) + (n-k) = n+1 > n$, dus de twee ruimten hebben een niet-triviale doorsnede.

**Stap 3:** Schrijf $w = \sum_{i=1}^{k+1} c_i v_i$ (met $\sum c_i^2 = 1$ want $\{v_i\}$ orthonormaal en $\|w\|_2 = 1$).

**Stap 4:** Omdat $w \in \mathcal{N}(Y^*)$, geldt $Bw = XY^*w = 0$, dus:

$$\|B - A\|_2 \geq \frac{\|(B-A)w\|_2}{\|w\|_2} = \|Aw\|_2.$$

**Stap 5:** Bereken $\|Aw\|_2$:

$$Aw = \sum_{i=1}^{k+1} c_i A v_i = \sum_{i=1}^{k+1} c_i \sigma_i u_i.$$

Omdat $\{u_i\}$ orthonormaal:

$$\|Aw\|_2^2 = \sum_{i=1}^{k+1} c_i^2 \sigma_i^2 \geq \sigma_{k+1}^2 \sum_{i=1}^{k+1} c_i^2 = \sigma_{k+1}^2.$$

**Conclusie:** $\|B - A\|_2 \geq \sigma_{k+1} = \|A_k^{svo} - A\|_2$, dus $A_k^{svo}$ is optimaal. ∎

### Bewijs (Frobeniusnorm, schets)

Via de eigenschap dat voor elke $B \in \mathcal{M}_k(m,n)$ de $i$-de singuliere waarde van $(B-A)$ voldoet aan $\sigma_i(B-A) \geq \sigma_{k+i}(A)$, volgt:

$$\|B - A\|_F^2 = \sum_i \sigma_i(B-A)^2 \geq \sum_{i=1}^{p-k} \sigma_{k+i}(A)^2 = \|A_k^{svo} - A\|_F^2. \quad \square$$

### Vergelijking van de twee methodes

| | QR met kolomverwisselingen ($A_k^{qr}$) | Getrunceerde SVD ($A_k^{svo}$) |
|---|---|---|
| **Kwaliteit** | Goed, maar niet optimaal | Optimaal (Stelling 12.2.1) |
| **Rekenkost** | Significant goedkoper | Duurder |
| **Garantie** | $\|A_k^{qr} - A\| \geq \|A_k^{svo} - A\|$ | Minimale fout gegarandeerd |
| **Gebruik** | Snelle benadering bij grote matrices | Wanneer optimale kwaliteit vereist is |

### Voorbeeld — beeldcompressie

Een zwart-wit beeld van $512 \times 512$ pixels is een matrix $A \in \mathbb{R}^{512 \times 512}$ (grijswaarden 0–255). De benadering $A_{20}^{svo}$ (rang 20) vereist slechts $512 \cdot 20 + 20 + 20 \cdot 512 = 20\,500$ getallen i.p.v. $262\,144$ — een factor $\approx 13$ compressie, met visueel aanvaardbare kwaliteit.

---

## 12.3 Principal Component Analysis (PCA)

### Opzet

Gegeven een **datamatrix** $X \in \mathbb{R}^{m \times n}$ met $m \geq n$:
- Elke **rij** is een experiment/observatie: $m$ observaties van $n$ variabelen.
- Elke **kolom** is een variabele (sensor, feature, meting).

**Veronderstelling 12.1:** Het steekproefgemiddelde van elke kolom is nul (eventueel te bereiken door verschuiving).

De $m$ observaties vormen een **puntenwolk** in $\mathbb{R}^n$. PCA zoekt de richtingen in die ruimte waarlangs de puntenwolk het meest spreidt.

### Principal components via de SVD

**Eerste hoofdrichting $r_1$** (maximale variantie):

$$r_1 = \arg\max_{\|r\|_2 = 1} \sum_{i=1}^m (r^T x^{(i)})^2 = \arg\max_{\|r\|_2 = 1} \|Xr\|_2^2 = \arg\max_{\|r\|_2 = 1} r^T X^T X r.$$

Gebruik makend van de SVD $X = U \Sigma V^T$ en de substitutie $r = Vy$:

$$\max_{\|y\|_2 = 1} y^T \Sigma^T \Sigma y,$$

met oplossing $\hat{y} = e_1$ en dus $r_1 = v_1$ — de rechter singuliere vector bij $\sigma_1$.

**Tweede hoofdrichting $r_2$** (maximale variantie, orthogonaal tot $r_1$):

$$r_2 = \arg\max_{\|r\|_2 = 1,\, r \perp v_1} r^T X^T X r \implies r_2 = v_2.$$

**Algemeen:** De $j$-de hoofdrichting is $r_j = v_j$, de $j$-de rechter singuliere vector van $X$.

> ⚠️ **Belangrijk voor examen:** De principal components zijn precies de rechter singuliere vectoren van de datamatrix $X$, gerangschikt naar afnemende singuliere waarde.

### Coördinaten in de nieuwe basis

Uit de SVD $X^T = \sum_{i=1}^n \sigma_i v_i u_i^T$ volgt dat de $j$-de datapunt $x^{(j)}$ kan geschreven worden als:

$$x^{(j)} = \sum_{i=1}^n \underbrace{(\sigma_i u_{ij})}_{\text{coördinaat}} v_i. \tag{12.10}$$

De coördinaten in de basis $\{v_1, \ldots, v_n\}$ zijn de linker singuliere vectoren gewogen met de singuliere waarden.

### Dimensiereductie

Door slechts de eerste $k$ termen te behouden, correspondeert de gereduceerde dataset precies met de getrunceerde SVD $X_k^{svo}$:

$$X \approx X_k^{svo} = \sum_{i=1}^k \sigma_i u_i v_i^T.$$

**Toepassing:** Als we een functie $y = f(x_1, \ldots, x_n)$ willen schatten uit meetdata, maar de puntenwolk leeft bijna volledig in een $k$-dimensionale deelruimte ($k \ll n$), dan volstaat het om slechts $k$ variabelen te beschouwen: $y \approx F(w_1, \ldots, w_k)$ met $w_j = v_j^T x$. Dit reduceert de complexiteit drastisch en de hoeveelheid benodigde trainingsdata.

---

## 12.4 Andere ijle representaties en benaderingen

Naast de SVD bestaan er transformaties gebaseerd op frequentie- of resolutiestructuur die voor bepaalde typen data (beelden, signalen) tot efficiëntere ijle representaties leiden.

### 12.4.1 Fourier-benaderingen — de Discrete Cosinustransformatie (DCT)

De **tweedimensionale DCT** ontbindt de matrix $A \in \mathbb{R}^{m \times n}$ in periodieke basisfuncties:

$$A_{ij} = \sum_{k=0}^{m-1} \sum_{l=0}^{n-1} c_k d_l B_{kl} \cos\!\left(\frac{(i+\tfrac{1}{2})k\pi}{m}\right) \cos\!\left(\frac{(j+\tfrac{1}{2})l\pi}{n}\right), \tag{12.11}$$

met $c_k = \tfrac{1}{\sqrt{m}}$ voor $k=0$ en $\sqrt{\tfrac{2}{m}}$ anders (analoog voor $d_l$).

De **coëfficiënten** $B_{kl}$ worden berekend als de gewogen innerproducten van $A$ met de cosinusbasisf uncties.

**Sleutelobservatie:** De transformatie van $A$ naar $B$ brengt **spaarsheid** aan het licht — de meeste coëfficiënten zijn klein in vergelijking met de dominante paar. Door enkel de $K$ grootste coëfficiënten te behouden en de inverse DCT toe te passen, bekomt men een goede benadering van $A$.

**Praktische eigenschappen:**
- De DCT (en haar inverse) zijn rekenkundig efficiënt via de **Fast Fourier Transform (FFT)** als $m$ en $n$ machten van twee zijn.
- De originele **JPEG-standaard** past de DCT toe op $8 \times 8$ blokken van het beeld.
- Voor gelijke compressieverhouding is de DCT efficiënter dan lage-rangbenaderingen voor typische fotografie.

### 12.4.2 Wavelets

Een **wavelettransformatie** geeft een ontbinding op basis van **resolutie**: alle basisfuncties zijn verschoven en/of geschaalde kopieën van één **mother wavelet** met compact drager.

**Haar-wavelettransformatie (eenvoudigste geval):**

Veronderstel $m = n$ en $m$ een macht van 2. De transformatie verloopt stapsgewijs:

**Stap 1:**
- Middel de grijswaarden in blokken van $2 \times 2$ pixels → benadering met resolutie $m/2 \times m/2$ ($m^2/4$ coëfficiënten)
- Sla de **verschilwaarden** op die reconstructie mogelijk maken ($3m^2/4$ coëfficiënten)

**Stap 2:** Herhaal stap 1 op de $m/2 \times m/2$ uitgemiddelde matrix → resolutie $m/4 \times m/4$.

**Itereer** tot volledig uitgemiddeld. Het resultaat is een hiërarchische structuur waarbij grovere resoluties hersteld worden uit steeds minder coëfficiënten.

**Ijle representatie:** Typisch zijn de meeste verschilcoëfficiënten klein (vlakke gebieden in een beeld). Door kleine coëfficiënten op nul te zetten en de inverse transformatie toe te passen, bekomt men een gecomprimeerde benadering.

**Standaarden:** De **JPEG2000-standaard** gebruikt de Daubechies-wavelettransformatie, die betere visuele kwaliteit geeft bij hoge compressie dan de JPEG-standaard (DCT).

---

## Samenvatting en vergelijking

| Methode | Basisidee | Sterktes | Zwaktes |
|---|---|---|---|
| **SVD → getrunceerde SVD** | Rang-1 ontbinding, knip af na $k$ termen | Optimaal (Eckart-Young-Mirsky); PCA; universeel | Duur om te berekenen |
| **QR met kolomverwisselingen** | Pivot op grootste kolom, knip af na $k$ rijen/kolommen | Goedkoper dan SVD | Suboptimaal |
| **DCT** | Ontbinding in cosinusbasisfuncties | Efficiënt voor beelden (JPEG); FFT-gebaseerd | Niet optimaal voor niet-periodieke data |
| **Wavelets** | Multi-resolutie ontbinding | Goede lokalisatie in tijd én frequentie (JPEG2000) | Complexer te begrijpen en implementeren |

**Verbanden:**
- PCA = getrunceerde SVD van de gecentreerde datamatrix $X$.
- De beste rang-$k$ benadering in spectrale en Frobenius-norm wordt steeds gegeven door de getrunceerde SVD (Stelling 12.2.1 — Eckart–Young–Mirsky).
- DCT en wavelets zijn beter dan SVD voor beeldcompressie omdat ze de **lokale** structuur van een beeld benutten (blokken, randen), terwijl de SVD de **globale** structuur benut.
