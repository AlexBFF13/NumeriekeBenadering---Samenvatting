# Deel 9 – Eigenwaardenalgoritmes

## Overzicht

Dit hoofdstuk behandelt twee klassen van problemen:

| Klasse | Vraag | Methodes |
|---|---|---|
| **Selecte eigenwaarden** | Vind de $k$ meest extreme eigenwaarden ($k \ll m$) | Methode van de machten, deelruimte-iteratie, Arnoldi/Lanczos |
| **Alle eigenwaarden** | Vind het volledige spectrum van een matrix | Householder → Hessenberg, impliciet verschoven QR-algoritme |

**Fundamentele barrière (Stelling 9.0.1, Abel):** Er bestaat geen gesloten formule voor de wortels van een algemeen veelterm van graad ≥ 5. Eigenwaarden berekenen is dus **inherent iteratief**.

De rode draad: grote matrices worden nooit volledig in driehoeksvorm gebracht in één stap. Alle methodes zijn iteratief — ze verfijnen een beginschatting tot convergentie.

---

## 9.1 Berekening van selecte eigenwaarden

### 9.1.1 Methode van de machten

**Intuïtie:** Stel een vector $\mathbf{x}^{(0)}$ voor als lineaire combinatie van eigenvectoren. Herhaal de vermenigvuldiging met $A$. De component in de richting van de dominante eigenvector groeit het snelst — uiteindelijk domineert die component volledig.

**Basisiteratie:** voor $A \in \mathbb{C}^{m \times m}$, startvector $\mathbf{x}^{(0)}$:

$$\mathbf{p}^{(k)} = A\mathbf{x}^{(k-1)}, \qquad \mathbf{x}^{(k)} = \frac{\mathbf{p}^{(k)}}{\|\mathbf{p}^{(k)}\|}, \quad k = 1, 2, 3, \ldots \tag{9.4}$$

**Convergentieanalyse:**

Schrijf $\mathbf{x}^{(0)} = \sum_{i=1}^m c_i \mathbf{x}_i$ (eigenvectorbasis, $c_1 \neq 0$). Dan:

$$\frac{1}{c_1 \lambda_1^k} A^k \mathbf{x}^{(0)} = \mathbf{x}_1 + \sum_{i=2}^m \frac{c_i}{c_1} \left(\frac{\lambda_i}{\lambda_1}\right)^k \mathbf{x}_i \xrightarrow{k \to \infty} \mathbf{x}_1,$$

mits $|\lambda_1| > |\lambda_2|$. De convergentie is **lineair** met convergentiefactor:

$$\rho = \left|\frac{\lambda_2}{\lambda_1}\right|.$$

Hoe kleiner de kloof tussen $|\lambda_1|$ en $|\lambda_2|$, hoe trager de convergentie.

**Eigenwaardeschatting via Rayleigh-quotiënt:**

$$r(\mathbf{x}) = \frac{\mathbf{x}^* A \mathbf{x}}{\mathbf{x}^* \mathbf{x}} \xrightarrow{k \to \infty} \lambda_1.$$

> ⚠️ **Belangrijk voor examen:** Ken de basisiteratie, de convergentiefactor $|\lambda_2/\lambda_1|$, en de varianten hieronder.

**Varianten:**

| Variant | Aanpassing | Convergeert naar |
| --- | --- | --- |
| **Inverse iteratie** | Vervang $A$ door $A^{-1}$ (los $A\mathbf{p}^{(k)} = \mathbf{x}^{(k-1)}$ op) | Eigenvector bij kleinste $\lvert\lambda\rvert$ |
| **Inverse iteratie met verschuiving $\sigma$** | Vervang $A$ door $(A - \sigma I)^{-1}$ | Eigenvector bij eigenwaarde dichtstbij $\sigma$ |
| **Rayleigh-quotiënt iteratie** | Pas $\sigma$ elke stap aan: $\sigma = r(\mathbf{x}^{(k-1)})$ | Asymptotisch **kubische** convergentie |

### 9.1.2 Deelruimte-iteratie

**Intuïtie:** In plaats van op één vector te itereren, itereer je op een $n$-dimensionale deelruimte. De deelruimte trekt aan naar de ruimte opgespannen door de $n$ dominante eigenvectoren.

**Basisiteratie:** start met $Q_0 \in \mathbb{C}^{m \times n}$ (orthonormale kolommen), dan:

$$P_k = A Q_{k-1}, \qquad Q_k R_k = P_k \quad (\text{onvolledige QR-factorisatie}), \quad k = 1, 2, \ldots \tag{9.5}$$

De normalisatie van de methode van de machten (schaardelen) wordt vervangen door QR-factorisatie: die herstelt zowel de lengte als de orthogonaliteit van de kolommen.

**Convergentie:** als $|\lambda_n| > |\lambda_{n+1}|$:

$$\lim_{k \to \infty} \langle \mathbf{q}_1^{(k)}, \ldots, \mathbf{q}_n^{(k)} \rangle = \langle \mathbf{x}_1, \ldots, \mathbf{x}_n \rangle, \tag{9.6}$$

met convergentiefactor $|\lambda_{n+1}/\lambda_n|$.

**Verband met methode van de machten:** de eerste kolom van $Q_k$ in deelruimte-iteratie gedraagt zich als de methode van de machten voor één vector. Het geneste karakter van QR (§2.3.1) zorgt dat tegelijk ook de eerste $j$ kolommen convergeren naar $\langle \mathbf{x}_1, \ldots, \mathbf{x}_j \rangle$.

**Eigenwaardeschatting:** als $Q_k \approx [\mathbf{x}_1 \cdots \mathbf{x}_n]$, dan is $Q_k^* A Q_k$ een $n \times n$ matrix waarvan de eigenwaarden de $n$ dominante eigenwaarden van $A$ benaderen. Los dit kleine eigenwaardenprobleem op met het algoritme van §9.2.

### 9.1.3 De methode van Arnoldi

**Kernidee:** de methode van de machten genereert in $k$ stappen vectoren $\{\mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots\}$, maar gebruikt alleen de laatste. Arnoldi benut *alle* gegenereerde vectoren door benaderingen van eigenvectoren te zoeken in de volledige **Krylov-ruimte**:

$$\mathcal{K}_k(A, \mathbf{b}) = \langle \mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots, A^{k-1}\mathbf{b} \rangle. \tag{9.7}$$

In de $k$-de stap worden **$k$ eigenwaardbenaderingen** bekomen — veel meer informatie per iteratie dan de methode van de machten.

#### 9.1.3.1 Orthonormale basis van de Krylov-ruimte

De vectoren $\{\mathbf{b}, A\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}$ zijn een slechte basis: voor grote $k$ zijn ze bijna allemaal evenredig met de dominante eigenvector $\mathbf{x}_1$ (gevaarlijke aftrekking). Oplossing: bouw via **Gram-Schmidt** een geneste orthonormale basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_k\}$ van $\mathcal{K}_k(A, \mathbf{b})$.

**Arnoldi's algoritme (Algoritme 6):**
```
q₁ = b / ||b||₂
for j = 1, 2, 3, ...:
    vⱼ = A qⱼ
    for i = 1 to j:
        hᵢⱼ = qᵢ* vⱼ
        vⱼ = vⱼ - hᵢⱼ qᵢ        ← Gram-Schmidt: verwijder component langs qᵢ
    h_{j+1,j} = ||vⱼ||₂
    q_{j+1} = vⱼ / h_{j+1,j}
```

**Recursievergelijking van Arnoldi:** de uitkomst van $k$ iteraties kan samengevat worden als:

$$A Q_k = Q_k H_k + h_{k+1,k}\, \mathbf{q}_{k+1} \mathbf{e}_k^T, \tag{9.8}$$

waarbij:
- $Q_k = [\mathbf{q}_1 \cdots \mathbf{q}_k] \in \mathbb{C}^{m \times k}$ (orthonormale kolommen)
- $H_k \in \mathbb{C}^{k \times k}$: **Hessenberg-matrix** (bovendriehoeks + één subdiagonaal)
- $h_{k+1,k}$: de laatste subdiagonaalcoëfficiënt

> ⚠️ **Belangrijk voor examen:** Ken de recursievergelijking (9.8) en zijn betekenis. Weet dat $Q_k^* A Q_k = H_k$ (volgt direct uit (9.8) na links vermenigvuldigen met $Q_k^*$).

#### 9.1.3.2 Ritz-waarden: eigenwaarden benaderen

We zoeken een benadering $(\theta, \mathbf{q})$ van een eigenwaardepaar met $\mathbf{q} \in \mathcal{K}_k$, d.w.z. $\mathbf{q} = Q_k \mathbf{y}$. Eis dat het residu **orthogonaal** is op $\mathcal{K}_k$:

$$Q_k^*(A - \theta I)Q_k \mathbf{y} = \mathbf{0} \implies H_k \mathbf{y} = \theta \mathbf{y}.$$

De **Ritz-waarden** zijn de eigenwaarden van de kleine $k \times k$ matrix $H_k$. In stap $k$ lost men dit kleine eigenwaardenprobleem op (met het algoritme van §9.2) en verkrijgt $k$ benaderingen van eigenwaarden van $A$.

**Voordeel:** in stap $k$ berekent men een $k \times k$ eigenwaardenprobleem in plaats van een $m \times m$ probleem, met $k \ll m$.

#### 9.1.3.3 Prematuur afbreken

Als $h_{k+1,k} = 0$, kan het algoritme niet verder. Maar dit is eigenlijk goed nieuws: dan geldt $A Q_k = Q_k H_k$, d.w.z. $\mathcal{K}_k(A,\mathbf{b})$ is een **$A$-invariante deelruimte**. Voor diagonaliseerbare $A$ betekent dit dat de Ritz-waarden exact de eigenwaarden zijn, en $\mathbf{b}$ heeft geen componenten buiten $\langle \mathbf{x}_1, \ldots, \mathbf{x}_k \rangle$.

#### 9.1.3.4 Convergentie-eigenschappen

**Wat convergeert eerst?** De **extreme eigenwaarden** — die welke kwalitatief goed gescheiden zijn aan de rand van het spectrum. Dit zijn niet noodzakelijk de grootste in absolute waarde, maar de meest geïsoleerde.

> ⚠️ **Belangrijk voor examen:** Ken de twee stellingen hieronder.

### Stelling 9.1.2 — Invariantie voor verschuivingen

**Als Arnoldi toegepast op $(A, \mathbf{b})$ Ritz-waarden $\{r_1^{(k)}, \ldots, r_k^{(k)}\}$ geeft, dan geeft Arnoldi toegepast op $(A - \sigma I, \mathbf{b})$ Ritz-waarden $\{r_1^{(k)} - \sigma, \ldots, r_k^{(k)} - \sigma\}$.**

**Bewijs:** Uit (9.8): $(A - \sigma I)Q_k = Q_k(H_k - \sigma I) + h_{k+1,k}\mathbf{q}_{k+1}\mathbf{e}_k^T$, wat een geldige Arnoldi-recursie is voor $(A-\sigma I, \mathbf{b})$. Bovendien geldt $\mathcal{K}_k(A-\sigma I, \mathbf{b}) = \mathcal{K}_k(A, \mathbf{b})$ (want de Krylov-ruimte hangt niet af van $\sigma$). $\square$

**Gevolg voor convergentie:** Arnoldi differentieert niet tussen de kleinste en de grootste eigenwaarde (in tegenstelling tot de methode van de machten). Zowel $\lambda_{\min}$ als $\lambda_{\max}$ worden vroeg teruggevonden — dit is een groot voordeel.

### Stelling 9.1.3 — Relatie met veeltermbenadering

**De karakteristieke veelterm van $H_k$ is de unieke oplossing van:**

$$\min_{p \in \mathcal{M}_k} \|p(A)\mathbf{b}\|_2, \tag{9.10}$$

**waarbij $\mathcal{M}_k$ de verzameling is van alle monische veeltermen van graad $k$.**

**Bewijs (schets):**

**Stap 1:** Elke $p \in \mathcal{M}_k$ geeft $p(A)\mathbf{b} = A^k\mathbf{b} - Q_k\mathbf{y}$ voor een unieke $\mathbf{y} \in \mathbb{C}^k$ (want $A^k\mathbf{b}$ en $Q_k$ spannen samen een ruimte op die $\{A^k\mathbf{b} - Q_k\mathbf{y} : \mathbf{y} \in \mathbb{C}^k\}$ omvat).

**Stap 2:** Het minimalisatieprobleem (9.10) is dus equivalent met $\min_{\mathbf{y}} \|A^k\mathbf{b} - Q_k\mathbf{y}\|_2$ — een kleinste-kwadratenprobleem.

**Stap 3:** De karakteristieke veelterm $\hat{p}$ van $H_k$ voldoet aan $\hat{p}(H_k) = 0$ (Cayley-Hamilton). Hieruit volgt via de Arnoldi-recursie dat $Q_k^*(A^k\mathbf{b} - Q_k\hat{\mathbf{y}}) = 0$: het residu is orthogonaal op $\mathcal{K}_k$. Per de karakterisering van de kleinste-kwadratenoplossing (§2.4.1) is dit de optimale $\mathbf{y}$. $\square$

**Interpretatie:** De Ritz-waarden (nulpunten van de karakteristieke veelterm van $H_k$) zijn zo gekozen dat $|p(\lambda_i)|$ klein is voor zoveel mogelijk eigenwaarden $\lambda_i$ van $A$. De geïsoleerde extreme eigenwaarden worden als eerste "opgepikt" als afzonderlijke nulpunten.

**Herschrijving via eigenwaardenontbinding (9.12):** Gebruik makend van (9.3), d.w.z. $p(A) = X\,\text{diag}(p(\lambda_1), \ldots, p(\lambda_m)) X^{-1}$, kunnen we het minimiseringsrobleem (9.10) herschrijven als:

$$\min_{p \in \mathcal{M}_k} \left\| X \begin{pmatrix} p(\lambda_1) \\ & \ddots \\ & & p(\lambda_m) \end{pmatrix} X^{-1} \mathbf{b} \right\|_2. \tag{9.12}$$

Dit maakt expliciet waarom **geclusterde of geïsoleerde eigenwaarden** de convergentie bepalen. Als één eigenwaarde $\lambda_j$ sterk geïsoleerd is van de rest, kan een nulpunt van $p$ dicht bij $\lambda_j$ gelegd worden zonder de andere termen te verstoren — en convergeert de overeenkomstige Ritz-waarde snel.

#### 9.1.3.5 Lanczos' algoritme voor Hermitiaanse matrices

Als $A = A^*$ (Hermitiaans), dan zijn ook de matrices $H_k$ Hermitiaans (want $H_k = Q_k^* A Q_k$). Een Hermitiaanse Hessenberg-matrix is **tridiagonaal**.

**Gevolg:** in het Gram-Schmidt-proces hoeft men $A\mathbf{q}_j$ slechts te orthogonaliseren t.o.v. $\mathbf{q}_{j-1}$ en $\mathbf{q}_j$ (niet t.o.v. alle eerdere vectoren) — het **Lanczos-algoritme** (Algoritme 7):

```
q₁ = b / ||b||₂
v₁ = Aq₁;  h₁₁ = q₁*v₁;  v₁ = v₁ - h₁₁q₁
h₂₁ = ||v₁||₂;  q₂ = v₁/h₂₁
for j = 2, 3, ...:
    vⱼ = Aqⱼ
    h_{j-1,j} = h_{j,j-1}           ← symmetrie!
    hⱼⱼ = qⱼ* vⱼ
    vⱼ = vⱼ - h_{j-1,j} q_{j-1} - hⱼⱼ qⱼ   ← slechts 2 termen!
    h_{j+1,j} = ||vⱼ||₂;  q_{j+1} = vⱼ/h_{j+1,j}
```

**Voordeel:** per stap slechts 2 inwendige producten en 2 vectorbewerkingen (in plaats van $j$ bij Arnoldi) → veel goedkoper.

**Toepassing:** Lanczos is het standaardalgoritme voor grote symmetrische/Hermitiaanse matrices (bijv. de graaf-Laplaciaan $L$ bij graafpartitionering uit §8.3).

---

## 9.2 Berekenen van alle eigenwaarden

Het standaardalgoritme voor het volledige spectrum: het **impliciet verschoven QR-algoritme**. Directe toepassing op kleine/middelgrote matrices ($m \sim \mathcal{O}(10^3)$); ook gebruikt voor het kleine eigenwaardenprobleem bij Arnoldi/deelruimte-iteratie.

Twee stappen:
1. **§9.2.1** Reduceer $A$ naar een gelijkvormige Hessenberg-matrix (eenmalig, direct)
2. **§9.2.2** Iteratief de eigenwaarden afsplitsen via QR-iteratie

### 9.2.1 Gelijkvormige Hessenberg-reductie

**Doel:** Breng $A$ via een tweezijdige gelijkvormigheidstransformatie naar een **Hessenberg-matrix** $H$ (bovendriehoeks + één subdiagonaal):

$$H = Q^* A Q, \quad Q \text{ unitair.} \tag{9.13}$$

**Waarom Hessenberg en niet driehoeks?** Bij een tweezijdige transformatie ($Q^* A Q$) kunnen we nullen creëren met de linker-multiplicatie, maar de rechter-multiplicatie gooit ze deels terug weg. We kunnen wel garanderen dat alles *onder* de eerste subdiagonaal nul wordt — dat geeft precies de Hessenberg-structuur.

**Methode:** Householder-transformaties, kolom per kolom (Algoritme 8):

- **Pre-multiplicatie** met $Q_{j+1,m}$: maakt nullen in kolom $j$ onder de subdiagonaal.
- **Post-multiplicatie** met $Q_{j+1,m}$: herstelt gelijkvormigheid (eigenwaarden blijven bewaard), maar introduceert geen nieuwe niet-nullen omdat $Q_{j+1,m}$ enkel kolommen $j+1, \ldots, m$ combineert en die al nul zijn in de reeds verwerkte rijen.

**Rekenkost:** $\sim \frac{10}{3}m^3$ bewerkingen (enkel $H$); $\sim \frac{16}{3}m^3$ (ook $Q$ expliciet).

**Verband met Arnoldi:** $m$ stappen van Arnoldi geven ook $AQ_m = Q_mH_m$, dus $Q_m^*AQ_m = H_m$ — dezelfde Hessenberg-reductie. Maar Algoritme 8 (Householder) is numeriek stabieler.

> ⚠️ **Belangrijk voor examen:** Begrijp waarom we Hessenberg en niet driehoeks kunnen bereiken met tweezijdige transformaties. Ken de structuur van het algoritme.

### 9.2.2 Iteratief bepalen van de eigenwaarden

#### Theoretische basis: deelruimte-iteratie met $p(A)$

Kies $s$ **verschuivingen** $\sigma_1, \ldots, \sigma_s \in \mathbb{C}$ (typisch $s = 1$ of $s = 2$) en definieer:

$$p(z) = \prod_{i=1}^s (z - \sigma_i).$$

Pas deelruimte-iteratie toe op de **volledige ruimte** ($n = m$) met matrix $p(A)$ in plaats van $A$. Dit verschuift de eigenwaarden via het spectral mapping theorem:

$$p(A) = X \begin{pmatrix} p(\lambda_1) & & \\ & \ddots & \\ & & p(\lambda_m) \end{pmatrix} X^{-1}.$$

Als de verschuivingen goede benaderingen zijn van $s$ eigenwaarden, dan zijn $|p(\lambda_i)|$ klein voor die eigenwaarden — hun eigenvectoren domineren niet meer, maar de andere wel. Convergentiefactor:

$$\rho = \left|\frac{p(\lambda_{m-s+1})}{p(\lambda_{m-s})}\right|.$$

Hoe beter de verschuivingen, hoe kleiner $\rho$ en hoe sneller de convergentie.

#### Eén iteratiestap van het impliciet verschoven QR-algoritme

**Het probleem:** deelruimte-iteratie met $p(A)$ op de volledige ruimte is rekenduur. Het impliciet verschoven QR-algoritme realiseert dit effect impliciet, door te exploiteren dat $A^{(k-1)}$ al in Hessenberg-vorm is.

**Procedure (één stap):**

**Stap 1:** Bereken $\mathbf{x} = p(A^{(k-1)})\mathbf{e}_1$ (geschaald zodat $\|\mathbf{x}\|_2 = 1$). Door de Hessenberg-structuur zijn enkel de eerste $s+1$ elementen van $\mathbf{x}$ niet-nul.

**Stap 2:** Bepaal een Householder-transformatie $Q_1$ die $\tilde{\mathbf{x}}$ (de eerste $s+1$ elementen van $\mathbf{x}$) afbeeldt op $\mathbf{e}_1$:
$$Q_1 = \begin{pmatrix} Q_{1,s+1} & 0 \\ 0 & I_{m-s-1} \end{pmatrix}.$$

**Stap 3:** Bereken $Q_1 A^{(k-1)} Q_1$. Dit verstoort de Hessenberg-structuur: er ontstaat een "uitstulping" (bulge).

**Stap 4:** Herstel de Hessenberg-vorm door opeenvolgende Householder-transformaties $Q_2 = Q_{2,s+2}Q_{3,s+3}\cdots Q_{m-1,m}$ die de uitstulping "langs de diagonaal jagen" tot ze uit de matrix verdwijnt:

$$A^{(k)} = Q_2^* Q_1^* A^{(k-1)} Q_1 Q_2 = Q^{(k)*} A^{(k-1)} Q^{(k)}.$$

**Sleuteleigenschap van $Q^{(k)}$:** de eerste kolom van $Q^{(k)}$ is gelijk aan $p(A^{(k-1)})\mathbf{e}_1$.

**Resultaat:**
- $A^{(k)}$ is **gelijkvormig** met $A^{(k-1)}$ (zelfde eigenwaarden)
- $A^{(k)}$ is ook in **Hessenberg-vorm**
- Via Stelling 9.2.5 is dit equivalent met een stap van deelruimte-iteratie met $p(A)$

> ⚠️ **Belangrijk voor examen:** Begrijp de logica van het algoritme: Hessenberg-reductie eenmalig, daarna iteratief de bulge doorkrijgen. Ken de eigenschappen van $A^{(k)}$.

#### Eigenwaarden afsplitsen (Stelling 9.2.4)

Als de $s$ verschuivingen goede benaderingen zijn van eigenwaarden $\lambda_{m-s+1}, \ldots, \lambda_m$, convergeert $|a_{m-s+1, m-s}^{(k)}| \to 0$. Zodra dit element (of een ander subdiagonaalelement) verwaarloosbaar klein is, kan de matrix gesplitst worden:

$$Q^*AQ \approx \begin{pmatrix} S_{11} & S_{12} \\ 0 & S_{22} \end{pmatrix},$$

waarbij $S_{22} \in \mathbb{C}^{s \times s}$ de afgesplitste eigenwaarden bevat. De eigenwaarden van $S_{22}$ dienen als nieuwe verschuivingen voor de volgende iteratiestap.

**Recursief algoritme (Algoritme 9):**
```
Gegeven Hessenberg-matrix A⁽⁰⁾:
for k = 1, 2, ...:
    Bereken A⁽ᵏ⁾ via één stap van het impliciet verschoven QR-algoritme
    Als een subdiagonaalelement kleiner is dan tolerantie:
        Splits de matrix op en pas het algoritme toe op het kleinere blok
```

**Keuze van verschuivingen:**
- $s = 1$: de verschuiving is het element $a_{mm}^{(k-1)}$ (onderaan rechts). Werkt goed voor reële eigenwaarden.
- $s = 2$: de verschuivingen zijn de eigenwaarden van de $2 \times 2$ deelmatrix rechts onderaan. Werkt voor reële matrices met complexe eigenwaardenparen (convergeert altijd).

**Convergentiesnelheid:** de subdiagonaalelementen convergeren typisch **kwadratisch** naar nul (zoals het geval $3.1876 \times 10^{-1} \to 2.7305 \times 10^{-1} \to \cdots \to 1.9477 \times 10^{-19}$ in het voorbeeld toont).

#### Stelling 9.2.5 — Verband met deelruimte-iteratie

> ⚠️ **Belangrijk voor examen:** Ken de stelling en het bewijs.

**Stelling:** Zij $A^{(k)}, A^{(k-1)} \in \mathbb{C}^{m \times m}$ Hessenberg-matrices en $Q^{(k)} \in \mathbb{C}^{m \times m}$ unitair die voldoen aan $A^{(k)} = Q^{(k)*} A^{(k-1)} Q^{(k)}$, waarbij de eerste kolom van $Q^{(k)}$ gelijk is aan $p^{(k)}(A^{(k-1)}) \mathbf{e}_1$. Veronderstel dat alle sub-diagonaalelementen van $A^{(k)}$ verschillend van nul zijn. Dan geldt:

$$\langle \mathbf{q}_1^{(k)}, \ldots, \mathbf{q}_j^{(k)} \rangle = p^{(k)}(A^{(k-1)}) \langle \mathbf{e}_1, \ldots, \mathbf{e}_j \rangle, \quad j \in \{1, \ldots, m\}.$$

### Bewijs van Stelling 9.2.5

**Doel:** Tonen dat de deelruimte opgespannen door de eerste $j$ kolommen van $Q^{(k)}$ gelijk is aan $p^{(k)}(A^{(k-1)})$ toegepast op $\langle \mathbf{e}_1, \ldots, \mathbf{e}_j \rangle$.

**Stap 1 — Herformuleer (9.18) als een Arnoldi-recursie.**

Schrijf $A^{(k)} = Q^{(k)*} A^{(k-1)} Q^{(k)}$ om als:

$$A^{(k-1)} Q^{(k)} = Q^{(k)} A^{(k)}.$$

Dit is precies de Arnoldi-recursievergelijking (9.8) voor de matrix $A^{(k-1)}$ en de startvector $\mathbf{q}_1^{(k)}$, omdat $A^{(k)}$ een Hessenberg-matrix is. Bijgevolg geldt voor elke $j$:

$$\langle \mathbf{q}_1^{(k)}, \ldots, \mathbf{q}_j^{(k)} \rangle = \mathcal{K}_j\!\left(A^{(k-1)},\, \mathbf{q}_1^{(k)}\right).$$

**Stap 2 — Substitueer de eerste kolom van $Q^{(k)}$.**

Omdat $\mathbf{q}_1^{(k)} = p^{(k)}(A^{(k-1)}) \mathbf{e}_1$:

$$\mathcal{K}_j\!\left(A^{(k-1)},\, p^{(k)}(A^{(k-1)}) \mathbf{e}_1\right) = p^{(k)}(A^{(k-1)})\, \mathcal{K}_j\!\left(A^{(k-1)},\, \mathbf{e}_1\right).$$

→ Dit volgt uit de commutatierelatie $A^{(k-1)} p^{(k)}(A^{(k-1)}) = p^{(k)}(A^{(k-1)}) A^{(k-1)}$: beide zijn veeltermen in dezelfde matrix.

**Stap 3 — Gebruik de Hessenberg-structuur van $A^{(k-1)}$.**

Omdat $A^{(k-1)}$ een Hessenberg-matrix is met alle sub-diagonaalelementen $\neq 0$, geldt:

$$\mathcal{K}_j\!\left(A^{(k-1)},\, \mathbf{e}_1\right) = \langle \mathbf{e}_1, \ldots, \mathbf{e}_j \rangle.$$

→ Immers, $A^{(k-1)} \mathbf{e}_1$ heeft een niet-nul tweede component (sub-diagonaal $\neq 0$), $A^{(k-1)2} \mathbf{e}_1$ heeft een niet-nul derde component, enz.

**Conclusie:** Samenvoeging van de drie stappen geeft het gestelde. $\square$

### Cumulatief effect van $k$ iteraties

Definieer de **cumulatieve transformatie** $\tilde{Q}^{(k)} = Q^{(1)} \cdots Q^{(k)}$, zodat:

$$A^{(k)} = \tilde{Q}^{(k)*} A^{(0)} \tilde{Q}^{(k)}.$$

Schrijf $\tilde{Q}^{(k)} = [\tilde{\mathbf{q}}_1^{(k)} \cdots \tilde{\mathbf{q}}_m^{(k)}]$. Dan geldt voor elke $j$:

$$\langle \tilde{\mathbf{q}}_1^{(k)}, \ldots, \tilde{\mathbf{q}}_j^{(k)} \rangle = \tilde{Q}^{(k-1)} \langle \mathbf{q}_1^{(k)}, \ldots, \mathbf{q}_j^{(k)} \rangle = \tilde{Q}^{(k-1)}\, p^{(k)}(A^{(k-1)})\, \langle \mathbf{e}_1, \ldots, \mathbf{e}_j \rangle.$$

Gebruik nu dat $A^{(k-1)} = \tilde{Q}^{(k-1)*} A^{(0)} \tilde{Q}^{(k-1)}$, zodat:

$$p^{(k)}\!\left(A^{(k-1)}\right) = \tilde{Q}^{(k-1)*}\, p^{(k)}\!\left(A^{(0)}\right)\, \tilde{Q}^{(k-1)}.$$

Substitueren en herhalen van hetzelfde argument voor $k-1, k-2, \ldots, 1$ geeft uiteindelijk:

$$\langle \tilde{\mathbf{q}}_1^{(k)}, \ldots, \tilde{\mathbf{q}}_j^{(k)} \rangle = p^{(k)}\!\left(A^{(0)}\right) \cdots p^{(1)}\!\left(A^{(0)}\right) \langle \mathbf{e}_1, \ldots, \mathbf{e}_j \rangle, \quad j = 1, \ldots, m. \tag{9.16'}$$

**Dit zegt ons:** $k$ stappen van het impliciet verschoven QR-algoritme zijn wiskundig equivalent met deelruimte-iteratie op de volledige ruimte $\mathbb{C}^m$, toegepast met de samengestelde matrix $p^{(k)}(A^{(0)}) \cdots p^{(1)}(A^{(0)})$. Hoe beter de verschuivingen eigenwaarden benaderen, hoe kleiner $|p^{(j)}(\lambda_i)|$ voor de te isoleren eigenwaarden, en hoe sneller de convergentie.

---

## Samenvatting en vergelijking van methodes

| Methode | Zoekt | Convergentie | Beste voor |
| --- | --- | --- | --- |
| **Methode van de machten** | 1 dominante eigenvector | Lineair, $\left\|\lambda_2/\lambda_1\right\|$ | Eén eigenwaarde, matrix-vector product goedkoop |
| **Inverse iteratie (+shift)** | Eigenwaarde dichtstbij $\sigma$ | Lineair (kubisch bij Rayleigh) | Verfijnen van een schatting |
| **Deelruimte-iteratie** | $n$ dominante eigenvectoren | Lineair, $\left\|\lambda_{n+1}/\lambda_n\right\|$ | Meerdere eigenwaarden |
| **Arnoldi** | $k$ extreme eigenwaarden | Kwalitatief lineair, extreme eerst | Grote ijle matrices, $k \ll m$ |
| **Lanczos** | Idem, voor Hermitiaanse $A$ | Idem, maar goedkoper per stap | Grote symmetrische matrices |
| **Impliciet QR** | Alle $m$ eigenwaarden | Kwadratisch (per blok) | Kleine/middelgrote volle matrices |

**Rekenkost impliciet QR-algoritme:** gedraagt zich in de praktijk als $\mathcal{O}(m^3)$ — dominiert door de eenmalige Hessenberg-reductie ($\frac{10}{3}m^3$).

### Verbinding met de rest van de cursus

- **Hoofdstuk 8:** PageRank (methode van de machten), meest centrale knoop (volledige eigenwaardenontbinding), graafpartitionering (Lanczos voor $\lambda_{N-1}$ van $L$)
- **Sectie 2.3:** Gram-Schmidt en QR-factorisatie zijn de bouwstenen van Arnoldi en de deelruimte-iteratie
- **Stelling 9.1.3** verbindt Arnoldi met de veeltermbenadering uit de rest van de cursus: de Ritz-waarden minimaliseren een polynoombenaderingsprobleem
