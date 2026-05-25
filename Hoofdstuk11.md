# Deel 11 – Optimalisatie-algoritmes

## Overzicht

Dit hoofdstuk behandelt numerieke methodes voor **onbeperkte optimalisatie**:

$$\min_{x \in \mathbb{R}^n} f(x),$$

waarbij $f : \mathbb{R}^n \to \mathbb{R}$ tweemaal continu differentieerbaar is. We beperken ons tot **eerste-orde methodes** die enkel gebruik maken van functie-evaluaties en gradënten (geen tweede afgeleiden).

| Sectie | Methode | Toepassing |
|---|---|---|
| 11.1 | Inleidende begrippen | Gradiënt, Hessiaan, min-/maxcondities |
| 11.2 | Steilste afdaling | Algemene niet-lineaire optimalisatie |
| 11.3 | Toegevoegde gradiënten (CG) | Lineaire stelsels $Ax = b$ met $A > 0$; ook algemene optimalisatie |
| 11.4 | Gauss-Newton | Niet-lineaire kleinste-kwadratenproblemen |
| 11.5 | Stochastische gradiënt | Training van neurale netwerken |

**Rode draad:** iteratieve methodes die in elke stap een **daalrichting** berekenen en dan een **stapgrootte** bepalen. De keuze van richting en stapgrootte bepaalt de convergentiesnelheid.

---

## 11.1 Inleidende begrippen

### Gradiënt en richtingsafgeleide

Voor $f : \mathbb{R}^n \to \mathbb{R}$ definiëren we de **gradiënt**:

$$\nabla f(x) := \left[\frac{\partial f}{\partial x_1}(x) \;\cdots\; \frac{\partial f}{\partial x_n}(x)\right]^T \in \mathbb{R}^n.$$

Voor $n = 1$ is de gradiënt gewoon de eerste afgeleide.

Voor een vectorfunctie $K : \mathbb{R}^n \to \mathbb{R}^m$ is de **Jacobiaan**:

$$J_K(x) := \left[\frac{\partial k_i}{\partial x_j}(x)\right]_{i,j=1}^{m,n} = \begin{bmatrix} \nabla k_1(x)^T \\ \vdots \\ \nabla k_m(x)^T \end{bmatrix} \in \mathbb{R}^{m \times n}.$$

Beschouw de rechte lijn $x(\alpha) := x + \alpha s$ met $s$ een eenheidsvector ($\|s\|_2 = 1$). De **richtingsafgeleide** van $f$ in richting $s$ in punt $x$ is:

$$\frac{dF(\alpha)}{d\alpha}\bigg|_{\alpha=0} = \sum_{i=1}^n \frac{\partial f}{\partial x_i}(x) \cdot s_i = s^T \nabla f(x) = (\nabla f(x))^T s. \tag{11.1}$$

### Eigenschap 11.1 — Steilste richting

> ⚠️ **Belangrijk voor examen:** Dit rechtvaardigt de keuze $p = -\nabla f(x)$ als daalrichting.

Van alle eenheidsvectoren $s$ is de **helling het grootst** (meest positief) in de richting van $+\nabla f(x)$ en het **kleinst** (meest negatief) in de richting van $-\nabla f(x)$.

**Bewijs:** Uit (11.1) volgt dat de richtingsafgeleide gelijk is aan $s^T \nabla f(x)$. Via de Cauchy-Schwarz ongelijkheid geldt:

$$s^T \nabla f(x) \leq \|s\|_2 \cdot \|\nabla f(x)\|_2 = \|\nabla f(x)\|_2,$$

met gelijkheid precies wanneer $s = \frac{\nabla f(x)}{\|\nabla f(x)\|_2}$. De minimale waarde wordt bereikt voor $s = -\frac{\nabla f(x)}{\|\nabla f(x)\|_2}$.

**Extra:** de gradiënt $\nabla f(x)$ staat **loodrecht** op de niveauverzameling $\{y \in \mathbb{R}^n : f(y) = f(x)\}$ door $x$.

### Hessiaan

De **Hessiaan** van $f$ is de matrix van tweede afgeleiden:

$$\nabla^2 f(x) := \left[\frac{\partial^2 f}{\partial x_i \partial x_j}(x)\right]_{i,j=1}^n \in \mathbb{R}^{n \times n}.$$

Deze matrix is **symmetrisch** (o.w.v. de stelling van Schwarz). De tweede richtingsafgeleide in richting $s$ is:

$$\frac{d^2 F(\alpha)}{d\alpha^2} = s^T \nabla^2 f(x) s,$$

en dit is de **kromming** van $f$ in punt $x$ in de richting $s$.

### Stellingen over minima

> ⚠️ **Belangrijk voor examen:** Ken de nodige en voldoende voorwaarden, inclusief het bewijs van de voldoende voorwaarde.

#### Stelling 11.1.1 — Nodige voorwaarde voor een minimum

Als $f$ een lokaal minimum heeft in $x^*$, dan geldt:

$$\nabla f(x^*) = 0.$$

**Opmerking:** $\nabla f(x^*) = 0$ is niet voldoende — het kan ook een maximum of zadelpunt zijn.

#### Stelling 11.1.2 — Voldoende voorwaarde voor een lokaal minimum

Als $\nabla f(x^*) = 0$ én de Hessiaan $\nabla^2 f(x^*)$ is **positief definiet**, dan heeft $f$ een geïsoleerd lokaal minimum in $x^*$.

### Bewijs (Stelling 11.1.2)

**Doel:** Tonen dat $x^*$ een geïsoleerd lokaal minimum is.

**Stap 1:** Uit $\nabla f(x^*) = 0$ volgt via (11.1) dat de richtingsafgeleide in $x^*$ nul is in **elke** richting $s$.
→ De eerste-orde term in de Taylorontwikkeling draagt niets bij.

**Stap 2:** Uit $\nabla^2 f(x^*) > 0$ (positief definiet) volgt dat $s^T \nabla^2 f(x^*) s > 0$ voor alle $s \neq 0$.
→ In **elke** richting is de kromming strikt positief: $f$ buigt overal omhoog.

**Conclusie:** Vanuit $x^*$ stijgt $f$ in elke richting, dus $x^*$ is een geïsoleerd lokaal minimum.

**Dit zegt ons:** Een kritiek punt waar de Hessiaan positief definiet is, is altijd een (strikt) lokaal minimum. Dit is de veelvariantse versie van "$f'(x^*) = 0$ en $f''(x^*) > 0 \Rightarrow$ minimum".

---

## 11.2 Methode van de steilste afdaling

### Afdalingsmethoden — algemene structuur

Een **iteratieve afdalingsmethode** verbetert $x^{(k-1)}$ in twee stappen:

1. Bereken een **daalrichting** $p^{(k)}$ zodanig dat $p^{(k)T} \nabla f(x^{(k-1)}) < 0$.
2. Bepaal een **stapgrootte** $\alpha^{(k)} > 0$ zodanig dat $x^{(k)} = x^{(k-1)} + \alpha^{(k)} p^{(k)}$ een lagere functiewaarde geeft.

De voorwaarde $p^{(k)T} \nabla f(x^{(k-1)}) < 0$ zegt dat de richtingsafgeleide negatief is: we gaan inderdaad omlaag (althans voor kleine stapgrootte).

Bij de **methode van de steilste afdaling** kiest men:

$$p^{(k)} = -\nabla f(x^{(k-1)}),$$

de richting waar de helling het meest negatief is (zie Eigenschap 11.1).

### Lijnzoekmethoden (bepaling van $\alpha^{(k)}$)

#### Exacte lijnzoekmethode

$$\alpha^{(k)} = \arg\min_{\alpha > 0} f(x^{(k-1)} + \alpha p^{(k)}).$$

Dit is optimaal langs de gekozen richting, maar vereist het oplossen van een scalair optimalisatieprobleem per iteratie — duur. Bovendien leidt de combinatie "steilste afdaling + exacte lijnzoeking" **niet** altijd tot snelle convergentie (de zigzag-beweging is beroemd traag).

#### Backtracking met de Armijo-voorwaarde

> ⚠️ **Belangrijk voor examen:** Ken de Armijo-voorwaarde en de intuïtie erachter.

In de praktijk gebruikt men **backtracking**: vertrek van $\alpha = 1$ en halveer (of vermenigvuldig met $q < 1$) tot de **Armijo-voorwaarde** voldaan is:

$$f\!\left(x^{(k-1)} + \alpha p^{(k)}\right) \leq f\!\left(x^{(k-1)}\right) + c \cdot \alpha \cdot \nabla f(x^{(k-1)})^T p^{(k)}, \tag{11.2}$$

met parameter $c \in (0, 1)$ (typisch klein, bv. $c = 10^{-4}$).

**Intuïtie:** De raaklijn aan $F(\alpha) = f(x^{(k-1)} + \alpha p^{(k)})$ in $\alpha = 0$ heeft helling $\nabla f(x^{(k-1)})^T p^{(k)} < 0$. Het rechterlid van (11.2) is die raaklijn, maar met helling vermenigvuldigd met $c < 1$ — een lichtjes minder steile neerwaartse lijn. De Armijo-voorwaarde eist dat $F(\alpha)$ **onder** die lijn ligt.

**Voordelen van Armijo:**
- Garandeert **monotoniciteit**: $f(x^{(k)}) < f(x^{(k-1)})$.
- Stopt zodra de voorwaarde voldaan is → vermijdt te kleine stappen (risico op stagnatie).

---

## 11.3 Algoritme van de toegevoegde gradiënten

> ⚠️ **Belangrijk voor examen:** Dit is het centrale algoritme van het hoofdstuk. Ken het algoritme, de eigenschappen (Krylov-ruimte, orthogonaliteit, A-norm minimaliteit) en het verband met veeltermbenadering.

### 11.3.1 Principe en basisiteratie

**Probleem met steilste afdaling:** De zoekrichting $p^{(k)} = -\nabla f(x^{(k)})$ houdt geen rekening met de vorige stappen. Dit leidt tot de beruchte "zigzag"-convergentie: opeenvolgende zoekrichtingen staan loodrecht op elkaar (bij exacte lijnzoeking), waardoor er veel iteraties verloren gaan.

**Idee van de toegevoegde gradiënten:** neem als zoekrichting een **lineaire combinatie** van de huidige gradiënt en de vorige zoekrichting:

$$p^{(k)} = -\nabla f(x^{(k)}) + \beta^{(k)} p^{(k-1)},$$

waarbij $\beta^{(k)}$ zo gekozen wordt dat de convergentie versneld wordt. Dit is slechts één extra berekening ten opzichte van steilste afdaling, maar de winst kan spectaculair zijn.

**Algoritme 10 — Toegevoegde gradiënten (basis)**

```
Vereist: x^(0)
p^(0) = -∇f(x^(0))
voor k = 1, 2, ...:
    α^(k) = argmin_α f(x^(k-1) + α·p^(k-1))     [exacte lijnzoeking]
    x^(k) = x^(k-1) + α^(k)·p^(k-1)
    β^(k) = ‖∇f(x^(k))‖² / ‖∇f(x^(k-1))‖²
    p^(k) = -∇f(x^(k)) + β^(k)·p^(k-1)
```

**Opmerking:** Door de exacte lijnzoeking is $p^{(k-1)}$ loodrecht op $\nabla f(x^{(k)})$ (de richting van steilste afdaling is loodrecht op de vorige zoekrichting). Dit is de sleutel tot de verbeterde convergentie.

### 11.3.2 Toepassing op lineaire stelsels

De methode van de toegevoegde gradiënten wordt vaak gebruikt voor het oplossen van **grote, ijle, symmetrisch positief definiete** (SPD) stelsels:

$$Ax = b, \quad A \in \mathbb{R}^{m \times m},\; A > 0. \tag{11.3}$$

**Verband met optimalisatie:** Als $A > 0$, dan is de oplossing $x^* = A^{-1}b$ het **unieke minimum** van de kwadratische doelfunctie:

$$f(x) := \frac{1}{2} x^T A x - x^T b. \tag{11.4}$$

**Verificatie:** $\nabla f(x) = Ax - b$ (zodat $\nabla f(x^*) = 0 \Leftrightarrow Ax^* = b$) en $\nabla^2 f(x) = A > 0$ (Stelling 11.1.2 is voldaan).

**Uitwerking van de stappen voor kwadratische $f$:**

Bij een kwadratische doelfunctie kunnen de stappen van Algoritme 10 analytisch uitgewerkt worden. We introduceren het **residu**:

$$r^{(k)} = b - Ax^{(k)}.$$

Merk op dat $r^{(k)} = -\nabla f(x^{(k)})$ (de negatieve gradiënt = het residu).

**Stap exacte lijnzoeking** (min langs $p^{(k-1)}$): de richtingsafgeleide nul stellen geeft:

$$\alpha^{(k)} = \frac{p^{(k-1)T} r^{(k-1)}}{p^{(k-1)T} A p^{(k-1)}}. \tag{11.5}$$

**Update van $x$ en $r$:**

$$x^{(k)} = x^{(k-1)} + \alpha^{(k)} p^{(k-1)}, \qquad r^{(k)} = r^{(k-1)} - \alpha^{(k)} A p^{(k-1)}.$$

**Parameter $\beta^{(k)}$** (met $\nabla f(x^{(k)}) = -r^{(k)}$):

$$\beta^{(k)} = \frac{r^{(k)T} r^{(k)}}{r^{(k-1)T} r^{(k-1)}}.$$

Met startwaarde $x^{(0)} = 0$ (zodat $r^{(0)} = b$ en $p^{(0)} = b$) komen we tot Algoritme 11 en de standaardvorm Algoritme 12.

**Algoritme 12 — Toegevoegde gradiënten voor stelsels (standaard)**

> ⚠️ **Belangrijk voor examen:** Ken dit algoritme in pseudocode en verstaan elke stap.

```
Vereist: A > 0, b, kmax, ε > 0
p^(0) = r^(0) = b,  x^(0) = 0
voor k = 1, 2, ..., kmax:
    α^(k) = (r^(k-1)ᵀ r^(k-1)) / (p^(k-1)ᵀ A p^(k-1))
    x^(k) = x^(k-1) + α^(k) · p^(k-1)
    r^(k) = r^(k-1) - α^(k) · A p^(k-1)
    als ‖r^(k)‖₂ < ε: geef x^(k) terug
    β^(k) = (r^(k)ᵀ r^(k)) / (r^(k-1)ᵀ r^(k-1))
    p^(k) = r^(k) + β^(k) · p^(k-1)
```

**Implementatie-opmerking:** per iteratie is slechts **één matrix-vectorproduct** $Ap^{(k-1)}$ nodig. Dit maakt de methode ideaal voor grote, ijle matrices — enkel het product $Av$ moet beschikbaar zijn, niet de matrix $A$ zelf.

### Eigenschappen van de iteratie

#### Eigenschap 11.2 — Krylov-ruimte

> ⚠️ **Belangrijk voor examen:** Toon aan dat de methode een Krylov-ruimte opbouwt.

Als $r^{(k-1)} \neq 0$, dan geldt na $k-1$ stappen van Algoritme 11:

$$\mathcal{K}_k(A, b) = \langle b, Ab, \ldots, A^{k-1}b \rangle = \langle x^{(1)}, \ldots, x^{(k)} \rangle = \langle p^{(0)}, \ldots, p^{(k-1)} \rangle = \langle r^{(0)}, \ldots, r^{(k-1)} \rangle. \tag{11.6}$$

**Interpretatie:** Na $k$ iteraties leeft de benaderende oplossing $x^{(k)}$ in de Krylov-ruimte $\mathcal{K}_k(A,b)$. Dit is dezelfde ruimte als bij het Arnoldi-algoritme (Sectie 9.1.3).

#### Eigenschap 11.3 — Orthogonaliteit van residu's en A-toegevoegdheid van zoekrichtingen

Als $r^{(k-1)} \neq 0$, dan geldt:

$$r^{(k)T} r^{(j)} = 0, \quad j = 0, \ldots, k-1, \tag{11.7}$$

$$p^{(k)T} A p^{(j)} = 0, \quad j = 0, \ldots, k-1. \tag{11.8}$$

**Interpretatie:**
- De residu's zijn **onderling loodrecht** (mutueel orthogonaal).
- De zoekrichtingen zijn **A-toegevoegd** (of A-geconjugeerd): $p^{(k)T} A p^{(j)} = 0$ voor $k \neq j$. Vandaar de naam "methode van de toegevoegde gradiënten" (conjugate gradients).

### Bewijs van Eigenschap 11.3 (door inductie)

**Doel:** Tonen dat (11.7) en (11.8) gelden voor alle $k \geq 1$, gegeven dat ze gelden voor kleinere indices.

**Basisgeval** $k = 1$: na te gaan door directe berekening.

**Inductiestap:** We veronderstellen dat (11.7)-(11.8) gelden voor indices $< k$.

**Bewijs van (11.7):**

Neem het inwendig product van de update $r^{(k)} = r^{(k-1)} - \alpha^{(k)} A p^{(k-1)}$ (lijn 5 van Alg. 12) met $r^{(j)}$:

$$r^{(k)T} r^{(j)} = r^{(k-1)T} r^{(j)} - \alpha^{(k)} r^{(j)T} A p^{(k-1)}. \tag{11.9}$$

*Geval $j < k-1$:* Uit (11.6) volgt $r^{(j)} \in \mathcal{K}_{j+1}(A,b) = \langle p^{(0)}, \ldots, p^{(j)} \rangle$. Beide termen in (11.9) zijn nul door de inductiehypothese.

*Geval $j = k-1$:* Het rechterlid van (11.9) is nul mits $\alpha^{(k)} = \frac{r^{(k-1)T} r^{(k-1)}}{r^{(k-1)T} A p^{(k-1)}}$. We tonen aan dat dit inderdaad de juiste $\alpha^{(k)}$ is.

Vermits $r^{(k-1)} = p^{(k-1)} - \beta^{(k-1)} p^{(k-2)}$ en $p^{(k-2)} \in \mathcal{K}_{k-1}(A,b) = \langle r^{(0)}, \ldots, r^{(k-2)} \rangle$, volgt uit de inductiehypothese:

$$r^{(k-1)T} r^{(k-1)} = r^{(k-1)T} p^{(k-1)}. \tag{11.11}$$

Op gelijkaardige wijze kan aangetoond worden dat $r^{(k-1)T} A p^{(k-1)} = p^{(k-1)T} A p^{(k-1)}$, zodat de formule voor $\alpha^{(k)}$ consistent is met (11.5). ∎

**Bewijs van (11.8):**

Neem het inwendig product van de update $p^{(k)} = r^{(k)} + \beta^{(k)} p^{(k-1)}$ met $Ap^{(j)}$:

$$p^{(k)T} A p^{(j)} = r^{(k)T} A p^{(j)} + \beta^{(k)} p^{(k-1)T} A p^{(j)}. \tag{11.12}$$

*Geval $j < k-1$:* Uit (11.6) volgt $Ap^{(j)} \in \mathcal{K}_{j+2}(A,b) = \langle r^{(0)}, \ldots, r^{(j+1)} \rangle$. Beide termen zijn nul door (11.7) en de inductiehypothese.

*Geval $j = k-1$:* Nul mits $\beta^{(k)} = -\frac{p^{(k-1)T} A r^{(k)}}{p^{(k-1)T} A p^{(k-1)}}$, wat men kan verifiëren is consistent met de definitie van $\beta^{(k)}$ in het algoritme. ∎

### Stelling 11.3.3 — Monotoniciteit van de fout (A-norm minimaliteit)

> ⚠️ **Belangrijk voor examen:** Dit is het kernesultaat over de optimale eigenschappen van CG.

Definieer de **fout** $e^{(k)} = x^{(k)} - x^*$ en de **A-norm** $\|x\|_A = \sqrt{x^T A x}$.

Als $r^{(k-1)} \neq 0$, dan geldt:

$$\|e^{(k)}\|_A^2 = \min_{x \in \mathcal{K}_k(A,b)} \|x - x^*\|_A^2, \tag{11.14}$$

en dus:

$$\|e^{(k)}\|_A \leq \|e^{(k-1)}\|_A \quad \text{(monotoon dalend)}. \tag{11.15}$$

Bovendien is $e^{(m)} = 0$ voor een zekere $m \leq \text{dim}$.

**Dit zegt ons:** $x^{(k)}$ is het element van de Krylov-ruimte $\mathcal{K}_k(A,b)$ dat de fout (gemeten met de A-norm) minimaliseert. De methode van de toegevoegde gradiënten is dus een **Krylov-methode** — net als Arnoldi, maar dan voor stelsels in plaats van eigenwaarden.

### Bewijs van Stelling 11.3.3

**Doel:** Tonen dat het minimum in (11.14) bereikt wordt voor $x = x^{(k)}$.

Elke $x \in \mathcal{K}_k(A,b)$ schrijven we als $x = x^{(k)} + \delta x$ met ook $\delta x \in \mathcal{K}_k(A,b)$. Dan:

$$\|x - x^*\|_A^2 = \|\delta x + e^{(k)}\|_A^2 = \|e^{(k)}\|_A^2 + 2 e^{(k)T} A \delta x + \|\delta x\|_A^2.$$

Nu is $A e^{(k)} = A(x^{(k)} - x^*) = -r^{(k)}$. Uit (11.7) volgt dat $r^{(k)} \perp \langle r^{(0)}, \ldots, r^{(k-1)} \rangle = \mathcal{K}_k(A,b)$, dus $e^{(k)T} A \delta x = -r^{(k)T} \delta x = 0$ voor alle $\delta x \in \mathcal{K}_k(A,b)$.

De gekruiste term is altijd nul, terwijl $\|\delta x\|_A^2 \geq 0$. Het minimum wordt bereikt voor $\delta x = 0$, d.w.z. $x = x^{(k)}$. ∎

### Stelling 11.3.4 — Verband met veeltermbenadering

> ⚠️ **Belangrijk voor examen:** Dit verband verklaart waarom geclusterde eigenwaarden snelle convergentie geven.

Als $r^{(k-1)} \neq 0$, dan geldt:

$$\|e^{(k)}\|_A^2 = \min_{p \in \mathcal{N}_k} \|p(A) e^{(0)}\|_A^2, \tag{11.17}$$

waarbij $\mathcal{N}_k$ de verzameling is van alle veeltermen van graad $\leq k$ die voldoen aan $p(0) = 1$.

### Bewijs van Stelling 11.3.4 (schets)

Elk element van $\mathcal{K}_k(A,b)$ heeft de vorm $c_1 b + c_2 Ab + \cdots + c_k A^{k-1}b$. Omdat $b = Ax^* = -Ae^{(0)}$, kan men schrijven:

$$x - x^* = e^{(0)} - c_1 A e^{(0)} - \cdots - c_k A^k e^{(0)} = \underbrace{(I - c_1 A - \cdots - c_k A^k)}_{\hat{p}(A)} e^{(0)} = \hat{p}(A) e^{(0)},$$

voor een veelterm $\hat{p} \in \mathcal{N}_k$ (want $\hat{p}(0) = 1$). Het gestelde volgt dan uit (11.14). ∎

**Implicaties voor convergentie:** Via de eigenwaardontbinding $A = X \Lambda X^{-1}$ (met $\Lambda = \text{diag}(\lambda_1, \ldots, \lambda_m)$) kan men (11.17) herschrijven als:

$$\|e^{(k)}\|_A^2 = \min_{p \in \mathcal{N}_k} \left\| \begin{pmatrix} p(\lambda_1) \\ \vdots \\ p(\lambda_m) \end{pmatrix} \right\|_A^2.$$

Hoe goed de convergentie na $k$ iteraties is, hangt af van hoe goed het spectrum van $A$ **afgedekt** kan worden door een veelterm $p \in \mathcal{N}_k$:

- Als $A$ slechts $k$ **verschillende eigenwaarden** heeft, convergeert de methode exact in $k$ iteraties (de veelterm met nulpunten in de eigenwaarden voldoet aan $p(0) = 1$).
- Als de eigenwaarden **geclusterd** zijn in weinig groepen, is ook snelle convergentie te verwachten.

### Preconditioning

Om de convergentie te versnellen kan men het stelsel (11.3) **voorconditoneren**. Kies een inverteerbare matrix $D$ en beschouw het equivalente stelsel:

$$(D^{-1} A D^{-T}) y = D^{-1} b, \quad \text{met } D^{-T} y = x.$$

Als $D^{-1} A D^{-T}$ beter geclusterde eigenwaarden heeft dan $A$, convergeert CG sneller. Matrix $D$ heet een **preconditioner**.

**Voorbeeld (diagonaalpreconditioner):** Als de diagonaalelementen van $A$ relatief groot zijn t.o.v. de niet-diagonaalelementen, kies $D_{ii} = \sqrt{A_{ii}}$. Dit schaalt alle diagonaalelementen naar 1 en concentreert de eigenwaarden rond 1.

---

## 11.4 De methode van Gauss-Newton

**Probleem:** Niet-lineaire kleinste-kwadratenproblemen van de vorm:

$$\min_c \|F(c)\|_2^2, \tag{11.20}$$

met $F : \mathbb{R}^n \to \mathbb{R}^m$ en $m > n$. Typisch voorbeeld: het fitten van een niet-lineair model $g(t_i; c_1, \ldots, c_n)$ aan meetdata $(t_i, y_i)$, waarbij $F_i(c) = y_i - g(t_i; c)$.

**Idee:** Lineariseer $F$ rondom de huidige schatting $c^{(k-1)}$:

$$F(c) \approx F(c^{(k-1)}) + J_F(c^{(k-1)}) \cdot (c - c^{(k-1)}),$$

met $J_F$ de Jacobiaan van $F$. Los het lineaire kleinste-kwadraatprobleem op:

$$\min_p \left\| F(c^{(k-1)}) + J_F(c^{(k-1)}) p \right\|_2^2, \tag{11.21}$$

en stel $c^{(k)} = c^{(k-1)} + p^{(k)}$ met $p^{(k)}$ de oplossing van (11.21).

**Dit is** een gewoon lineair kleinste-kwadratenprobleem met $A = J_F(c^{(k-1)})$ en $b = -F(c^{(k-1)})$, oplosbaar via de methodologie van Sectie 2.4.2 (QR-factorisatie).

**Opmerking:** Als $F$ lineair is, $F(c) = Ac - b$, dan is Gauss-Newton identiek aan het direct oplossen van $Ac = b$. Een convergentie-analyse valt buiten de scope van de cursus.

---

## 11.5 Optimalisatie bij neurale netwerktraining

De kostfunctie (10.14) voor het trainen van een neuraal netwerk met parameters $c \in \mathbb{R}^n$ en $N$ trainingspunten is:

$$\text{Kost}(c) = \frac{1}{N} \sum_{i=1}^N \frac{1}{2} \|y(x^{\{i\}}) - a^{[L]}(x^{\{i\}})\|_2^2.$$

### Gradiëntafdaling

De basisiteratie is:

$$c^{(k+1)} = c^{(k)} - \alpha \nabla \text{Kost}(c^{(k)}), \tag{11.22}$$

met $\alpha > 0$ de **learning rate** (leersnelheid). Men itereert tot een stopcriterium bereikt is.

**Berekening van de gradiënt — backpropagation:** De gradiënt $\nabla \text{Kost}(c)$ kan efficiënt berekend worden via **backpropagation**: een achterwaartse doorloop van het netwerk, die neerkomt op systematisch gebruik van de **kettingregel** voor afgeleiden.

### Stochastische gradiëntmethode (SGD)

> ⚠️ **Belangrijk voor examen:** Ken de varianten van SGD en de gebruikte terminologie.

**Probleem met (11.22):** Bij grote $N$ en $n$ vereist het berekenen van $\nabla \text{Kost}(c^{(k)}) = \frac{1}{N} \sum_i \nabla C_{x^{\{i\}}}(c^{(k)})$ een volledige pass over alle trainingsdata — dit is duur.

**Oplossing:** Vervang het gemiddelde over alle $N$ punten door de gradiënt bij één willekeurig gekozen trainingspunt:

1. Kies een willekeurig geheel getal $i \in \{1, \ldots, N\}$.
2. Update: $c^{(k+1)} \leftarrow c^{(k)} - \alpha \nabla C_{x^{\{i\}}}(c^{(k)})$. $(11.25)$

**Belangrijk:** SGD garandeert **niet** dat de kostfunctie bij elke stap daalt (door de willekeur). Vandaar de naam "stochastische gradiënt" (niet "stochastische steilste afdaling").

### Epochs en sampling zonder teruglegging

In plaats van sampling met teruglegging kan men **zonder teruglegging** werken: doorloop alle $N$ trainingspunten in willekeurige volgorde. Eén volledige pass heet een **epoch**:

1. Permuteer de indices $\{1, \ldots, N\}$ naar $\{k_1, k_2, \ldots, k_N\}$.
2. Voor $i = 1$ tot $N$: $c^{(k+1)} = c^{(k)} - \alpha \nabla C_{x^{\{k_i\}}}(c^{(k)})$. $(11.26)$

### Mini-batch SGD

Een compromis tussen volledige gradiënt en SGD met één punt: gebruik een **mini-batch** van $m$ willekeurig gekozen trainingspunten:

1. Kies $m$ indices $k_1, \ldots, k_m$ uniform en willekeurig uit $\{1, \ldots, N\}$.
2. Update: $c^{(k+1)} = c^{(k)} - \frac{\alpha}{m} \sum_{i=1}^m \nabla C_{x^{\{k_i\}}}(c^{(k)})$. $(11.27)$

De verzameling $\{x^{\{k_i\}}\}_{i=1}^m$ heet de **mini-batch**. De batchgrootte $m$ en de learning rate $\alpha$ zijn **hyperparameters** van het algoritme.

**Keuze van hyperparameters:** In de praktijk worden $m$, $\alpha$, en de randomisatiestrategie vaak bepaald door de vereisten van de rekenkrachtige hardware (GPU's zijn geoptimaliseerd voor bepaalde batchgroottes). Men kan $\alpha$ ook dynamisch aanpassen naarmate het trainingsproces vordert.

---

## Samenvatting

| Methode | Zoekrichting $p^{(k)}$ | Stapgrootte | Toepassing |
|---|---|---|---|
| Steilste afdaling | $-\nabla f(x^{(k)})$ | Armijo of exact | Algemene optimalisatie |
| Toegevoegde gradiënten | $-\nabla f(x^{(k)}) + \beta^{(k)} p^{(k-1)}$ | Exact | SPD stelsels; algemene opt. |
| Gauss-Newton | $\arg\min_p \|F(c^{(k)}) + J_F p\|^2$ | 1 (of met lijnzoeking) | Niet-lineaire KKM |
| SGD | $-\nabla C_{x^{\{i\}}}(c)$ | Learning rate $\alpha$ | Neurale netwerken |

**Verband Hoofdstuk 9 ↔ Hoofdstuk 11:** De methode van de toegevoegde gradiënten voor stelsels is een Krylov-methode, net als het Arnoldi-algoritme voor eigenwaarden. In beide gevallen wordt de benaderende oplossing gezocht in een Krylov-ruimte, maar de optimalisatiecriterium verschilt: Arnoldi minimaliseert het residu loodrecht op de ruimte, CG minimaliseert de fout in de A-norm.
