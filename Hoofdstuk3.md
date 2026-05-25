# Deel 3 – Benaderen van functies

## Overzicht en rode draad

In Hoofdstuk 2 werkten we in $\mathbb{R}^n$ of $\mathbb{C}^n$: vectoren met $n$ componenten. Nu willen we hetzelfde doen — benaderingen zoeken, projecties berekenen, Gram-Schmidt uitvoeren — maar in **abstracte ruimten**, zoals ruimten van functies.

De aanpak is als volgt: we bouwen een hiërarchie van steeds rijkere ruimten.

```
Metrische ruimte   ⊃   Genormeerde ruimte   ⊃   Unitaire ruimte
  (alleen afstand)       (norm → afstand)        (scalair product → norm → afstand)
```

Bij elke stap leggen we meer structuur op, en bij elke stap kunnen we meer resultaten bewijzen. Uiteindelijk in unitaire ruimten kunnen we de volledige theorie van Hoofdstuk 2 veralgemenen.

**Rode draad:** het zoeken van een **beste benadering** van een element in een deelverzameling of deelruimte.

---

## 3.1 Metrische ruimte en afstand

### Het begrip afstand

We kennen afstand intuïtief uit de meetkunde: $|x - y|$ op de getallenrechte, $\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2}$ in het vlak. Om dit bruikbaar te maken voor functieruimten, vatten we de *essentiële eigenschappen* van afstand samen in een axiomatische definitie.

> ⚠️ **Belangrijk voor examen:** De definitie van een afstandsfunctie (metriek) en haar vier axioma's.

### Definitie

Een **afstandsfunctie** (of metriek) op een verzameling $A$ is een functie $\delta: A \times A \to \mathbb{R}$ met de volgende eigenschappen voor alle $x, y, z \in A$:

$$
\begin{aligned}
1.\quad & \delta(x, y) \geq 0 \\
2.\quad & \delta(x, y) = 0 \iff x = y \\
3.\quad & \delta(x, y) = \delta(y, x) \quad \text{(symmetrie)} \\
4.\quad & \delta(x, y) \leq \delta(z, x) + \delta(z, y) \quad \text{(driehoeksongelijkheid)}
\end{aligned}
$$

De verzameling $A$ uitgerust met $\delta$ noemt men een **metrische ruimte**, genoteerd $(A, \delta)$.

**Intuïtie:** De eerste twee eigenschappen zeggen dat afstand altijd niet-negatief is en nul is precies wanneer twee punten samenvallen. De driehoeksongelijkheid zegt dat een directe weg nooit langer is dan een omweg via een derde punt.

---

### Het begrip beste benadering

**Intuïtie:** Gegeven een punt $x$ en een deelverzameling $D$, willen we het punt in $D$ vinden dat het dichtst bij $x$ ligt.

### Definitie

Zij $D$ een deelverzameling van een metrische ruimte $(A, \delta)$. Een element $d \in D$ is een **beste benadering** van $x \in A$ als

$$
\delta(x, D) := \inf\{\delta(x, y) : y \in D\} = \delta(x, d).
$$

Twee **waarschuwingen:**
- Een beste benadering **bestaat niet altijd** (bijv. $D = (0, 1)$, $x = 2$: de afstand $1$ wordt nooit bereikt).
- Een beste benadering is **niet altijd uniek** (bijv. alle punten op de eenheidscirkel zijn even ver van de oorsprong).

---

### Voorbeelden van afstandsfuncties

#### Continue functieruimten

Op de ruimte $C[a,b]$ van continue reële functies op $[a,b]$ zijn meerdere afstanden mogelijk.

**Kwadratische afstand ($L^2$-afstand):**

$$
\delta(f, g) = \left( \int_a^b w(x)\, |f(x) - g(x)|^2\, dx \right)^{1/2}
$$

met $w(x) > 0$ een gewichtsfunctie. Speciaal geval van de **$L^p$-afstand**:

$$
\delta_p(f, g) = \left( \int_a^b w(x)\, |f(x) - g(x)|^p\, dx \right)^{1/p}, \quad p \geq 1.
$$

**Chebyshev-afstand (maximumafstand):**

$$
\delta_\infty(f, g) = \max_{x \in [a,b]} |f(x) - g(x)|.
$$

**Belangrijk subtiel punt over $L^p[a,b]$:** De integraal

$$
\delta(f, g) = \left( \int_a^b |f(x) - g(x)|^p\, dx \right)^{1/p}
$$

definieert **geen** afstand als we met gewone functies werken: twee functies die in één punt verschillen hebben afstand nul maar zijn toch ongelijk. Oplossing: we identificeren functies die *bijna overal* (b.o.) gelijk zijn, d.w.z. die slechts op een verzameling van maat nul verschillen. Zo onstaat de ruimte $L^p[a,b]$ van *functieklassen*, uitgerust met de Lebesgue-integraal.

#### Discrete ruimten

Op een eindig stel punten $x_1, \ldots, x_n$ is $\delta(f, g) = \left(\sum_{i=1}^n |f_i - g_i|^p\right)^{1/p}$ een afstandsfunctie. De ruimte $\ell^p$ van $p$-sommeerbare rijen $\{x_i\}_{i \geq 1}$ wordt gemetriseerd door

$$
\delta(x, y) = \left( \sum_{i=1}^\infty |x_i - y_i|^p \right)^{1/p}.
$$

#### Niet-klassieke afstanden

- **Triviale afstand:** $\delta(x,y) = 1$ als $x \neq y$. Nuttig als tegenvoorbeeld, geen praktische toepassing.
- **Hamming-afstand:** het aantal posities waarop twee binaire tuples verschillen.
- **Silverman-afstand:** $\delta(A, B) = \#\{(A \cup B) \setminus (A \cap B)\}$ voor eindige verzamelingen.

---

## 3.2 Genormeerde ruimte en lengte

### Vectorruimte

Een **vectorruimte** $V$ over een veld $F$ ($= \mathbb{R}$ of $\mathbb{C}$) is een verzameling met optelling en scalaire vermenigvuldiging die de standaard 10 axioma's voldoet (associativiteit, commutatief, nulelement, invers, distributiviteit, etc.). Denk aan $\mathbb{R}^n$, $C[a,b]$, de ruimte van polynomen.

**Van metrisch naar genormeerd:** We leggen nu *meer* structuur op: de afstand moet compatibel zijn met de vectorstructuur.

> ⚠️ **Belangrijk voor examen:** De definitie van norm en haar vier axioma's.

### Definitie

Een **norm** op een vectorruimte $V$ is een functionaal $\|\cdot\|: V \to \mathbb{R}$ met de volgende eigenschappen:

$$
\begin{aligned}
1.\quad & \|x\| \geq 0 \\
2.\quad & \|x\| = 0 \iff x = 0 \\
3.\quad & \|ax\| = |a|\, \|x\| \quad \text{(homogeniteit)} \\
4.\quad & \|x + y\| \leq \|x\| + \|y\| \quad \text{(driehoeksongelijkheid)}
\end{aligned}
$$

Een vectorruimte met een norm heet een **genormeerde ruimte**.

Een nuttige afgeleide ongelijkheid:

$$
\big|\|x\| - \|y\|\big| \leq \|x - y\|.
$$

---

### Verband met metrische ruimten

> ⚠️ **Belangrijk voor examen:** Elke norm induceert een afstand, maar niet elke afstand komt van een norm.

### Stelling 3.2.1 — Norm induceert afstand

**Stelling:** Als $V$ een genormeerde ruimte is, dan is $V$ ook metrisch via $\delta(x, y) = \|x - y\|$.

**Intuïtie:** De afstand tussen $x$ en $y$ is de lengte van het verschilvector $x - y$.

### Bewijs

We moeten de vier eigenschappen van een afstandsfunctie verifiëren. Eigenschappen 1–3 volgen direct uit de normeigenschappen. Enkel de driehoeksongelijkheid vraagt werk.

**Te bewijzen:** $\|x - y\| \leq \|x - z\| + \|y - z\|$.

**Stap 1:** Stel $\alpha = x - z$ en $\beta = z - y$. Dan is $x - y = \alpha + \beta$.

**Stap 2:** De driehoeksongelijkheid van de norm geeft:
$$
\|x - y\| = \|\alpha + \beta\| \leq \|\alpha\| + \|\beta\| = \|x - z\| + \|z - y\|.
$$

**Conclusie:** De functie $\delta(x,y) = \|x - y\|$ voldoet aan alle vier axioma's. $\square$

**Dit zegt ons:** Elke genormeerde ruimte is automatisch ook een metrische ruimte. De omgekeerde implicatie is niet geldig: een metrische vectorruimte is slechts genormeerbaar als haar afstandsfunctie **translatie-invariant** ($\delta(x,y) = \delta(x+z, y+z)$) én **homogeen** ($\delta(ax, ay) = |a|\, \delta(x,y)$) is.

Afstanden die *niet* van een norm komen: chordale afstand, triviale afstand, Silverman-afstand.

---

### Voorbeelden van normen

**In $\mathbb{R}^n$ of $\mathbb{C}^n$:**

$$
\|x\|_\infty = \max_i |x_i|, \quad \|x\|_1 = \sum_i |x_i|, \quad \|x\|_p = \left(\sum_i |x_i|^p\right)^{1/p} \text{ voor } p \geq 1.
$$

**In $C[a,b]$ (continue functies):**

$$
\|f\|_\infty = \max_{x \in [a,b]} |f(x)|, \qquad \|f\|_p = \left(\int_a^b w(x)|f(x)|^p\, dx\right)^{1/p}.
$$

De driehoeksongelijkheid voor de $\infty$-norm in $C[a,b]$ bewijzen we:

Voor alle $x \in [a,b]$: $|f(x)| \leq \sup |f|$ en $|g(x)| \leq \sup |g|$. Dus
$$
|f(x) + g(x)| \leq |f(x)| + |g(x)| \leq \sup|f| + \sup|g|.
$$
Dit geldt voor alle $x$, dus $\sup|f+g| \leq \sup|f| + \sup|g|$.

De $L^p$-norm voldoet aan de driehoeksongelijkheid dankzij de **ongelijkheid van Minkowski**.

**Equivalentie van normen:** In een *eindigdimensionale* vectorruimte zijn alle normen equivalent: convergentie in de ene norm impliceert convergentie in elke andere norm.

---

### Beste benadering in een deelruimte: convexiteit

**Kernvraag:** Wanneer bestaat er een unieke beste benadering?

Het antwoord hangt af van de *geometrie* van de eenheidsbol.

### Definitie

Een verzameling $C$ is **convex** als voor alle $x_1, x_2 \in C$ en $\lambda \in (0,1)$:
$$
\lambda x_1 + (1-\lambda) x_2 \in C.
$$

$C$ is **strikt convex** als het inwendige lijnstuk volledig in het *inwendige* van $C$ ligt.

**Feit:** In een genormeerde ruimte is elke gesloten bol convex. Bewijs: voor $\|x_1 - a\| \leq r$ en $\|x_2 - a\| \leq r$ geldt
$$
\|\lambda x_1 + (1-\lambda) x_2 - a\| \leq \lambda\|x_1-a\| + (1-\lambda)\|x_2-a\| \leq r.
$$

### Definitie — Strikte norm

Een genormeerde ruimte is **strikt genormeerd** als haar eenheidsbol strikt convex is, d.w.z.:
$$
x \neq y, \quad \|x\| = \|y\| = 1 \implies \|x + y\| < 2.
$$

**Niet-strikte normen:** $\|\cdot\|_1$ en $\|\cdot\|_\infty$ in $\mathbb{R}^n$ zijn **niet** strikt: er bestaan koorden op het boloppervlak. De $\|\cdot\|_p$ voor $1 < p < \infty$ is wel strikt.

---

### Stelling 3.2.3 — Existentie en uniciteit van beste benadering

> ⚠️ **Belangrijk voor examen:** Stelling én bewijs (uniciteit via strikte convexiteit).

**Stelling:** Zij $D$ een **eindigdimensionale deelruimte** van een **strikt genormeerde ruimte** $V$, en zij $v \in V$. Dan **bestaat** de beste benadering van $v$ in $D$ en is ze **uniek**.

**Intuïtie:** Existentie volgt doordat eindigdimensionale deelruimten "gesloten" zijn (rijen in $D$ hebben grenswaarden in $D$). Uniciteit volgt doordat de strikte norm ervoor zorgt dat het middelpunt van twee beste benaderingen altijd beter is — een contradictie.

### Bewijs

**Existentie:**

Zij $d = \inf\{\|v - w\| : w \in D\}$. We tonen aan dat dit infimum bereikt wordt.

**Stap 1:** Er bestaat een rij $\{w_k\}_{k \geq 1}$ in $D$ met $\|v - w_k\| \to d$.

**Stap 2:** De rij $\{w_k\}$ is uniform begrensd, want:
$$
\|w_k\| \leq \|w_k - v\| + \|v\| \leq \|w_1 - v\| + \|v\|.
$$
→ de normen $\|w_k\|$ zijn begrensd.

**Stap 3:** Schrijf $w_k = \sum_{i=1}^n c_i^{(k)} a_i$ in de basis $\{a_1, \ldots, a_n\}$ van $D$. De coördinaatvectoren $(c_1^{(k)}, \ldots, c_n^{(k)})$ zijn begrensd in $\mathbb{R}^n$.

**Stap 4:** Volgens de stelling van Bolzano-Weierstrass heeft elke begrensde rij in $\mathbb{R}^n$ een convergente deelrij. Zonder beperking van de algemeenheid convergeert de rij naar $(\hat{c}_1, \ldots, \hat{c}_n)$.

**Stap 5:** Definieer $\hat{w} = \sum_{i=1}^n \hat{c}_i a_i \in D$. Dan:
$$
\|v - \hat{w}\| \leq \|v - w_k\| + \|w_k - \hat{w}\| \xrightarrow{k \to \infty} d + 0 = d.
$$
Dus $\|v - \hat{w}\| = d$ en $\hat{w}$ is een beste benadering.

**Uniciteit (bewijs uit het ongerijmde):**

Veronderstel dat $v_1, v_2 \in D$ met $v_1 \neq v_2$ beiden een beste benadering zijn:
$$
\|v - v_1\| = \|v - v_2\| = d.
$$

**Stap 1:** Definieer $e_i = \frac{1}{d}(v - v_i)$ voor $i = 1, 2$. Dan geldt $\|e_1\| = \|e_2\| = 1$.

**Stap 2:** Omdat de eenheidsbol strikt convex is ($v_1 \neq v_2$ impliceert $e_1 \neq e_2$):
$$
\left\|\frac{e_1 + e_2}{2}\right\| < 1.
$$

**Stap 3:** Dit geeft:
$$
\left\|v - \frac{v_1 + v_2}{2}\right\| = \left\|\frac{(v-v_1) + (v-v_2)}{2}\right\| = d \left\|\frac{e_1 + e_2}{2}\right\| < d.
$$

**Stap 4:** Maar $\frac{v_1 + v_2}{2} \in D$ (deelruimten zijn gesloten voor lineaire combinaties), en deze vector benadert $v$ beter dan $d$ — een contradictie met de definitie van $d$ als infimum.

**Conclusie:** De beste benadering is uniek. $\square$

**Dit zegt ons:** Zodra onze norm "strikt convex" is (geen rechte stukken op de eenheidsbol), bestaat er altijd precies één beste benadering in elke eindigdimensionale deelruimte.

---

## 3.3 Unitaire ruimte en orthogonaliteit

### Het begrip scalair product

**Intuïtie:** Een scalair product voegt nog meer meetkundige structuur toe: hoeken, loodrechtheid, en een voorkeursnorm. Het generaliseert het klassieke inwendig product $x \cdot y = \sum_i x_i y_i$ uit $\mathbb{R}^n$.

> ⚠️ **Belangrijk voor examen:** De axioma's van het scalair product (complexe versie).

### Definitie

Een vectorruimte $V$ over $\mathbb{C}$ is **unitair** als er een **scalair product** $(x, y): V \times V \to \mathbb{C}$ gedefinieerd is met de volgende eigenschappen:

$$
\begin{aligned}
1.\quad & (x, ay) = a(x,y) \quad \text{voor alle } a \in \mathbb{C} \quad \text{(lineair in 2e argument)} \\
2.\quad & (x + y, z) = (x, z) + (y, z) \quad \text{(additief in 1e argument)} \\
3.\quad & (x, y) = \overline{(y, x)} \quad \text{(Hermitische symmetrie)} \\
4.\quad & (x, x) > 0 \quad \text{als } x \neq 0 \quad \text{(positief definiet)}
\end{aligned}
$$

Voor een ruimte over $\mathbb{R}$ is eigenschap 3 gewone symmetrie $(x,y) = (y,x)$.

**Opmerking:** Uit 1 en 3 volgt $(ax, y) = \bar{a}(x,y)$ voor $a \in \mathbb{C}$: het scalair product is **conjugaat-lineair** in het eerste argument.

### Cauchy-Schwarz ongelijkheid

> ⚠️ **Belangrijk voor examen:** Kennen en kunnen gebruiken.

$$
|(x, y)| \leq \sqrt{(x,x)}\sqrt{(y,y)}.
$$

Dit wordt later, na de introductie van de geïnduceerde norm, geschreven als:

$$
|(x, y)| \leq \|x\|\, \|y\|.
$$

---

### Verband met genormeerde ruimten

De hiërarchie wordt nu completer:

- Elke unitaire ruimte is genormeerd.
- Niet elke genormeerde ruimte is unitair (alleen als de parallellogramgelijkheid geldt).
- Unitaire ruimten hebben een *strikte* norm, dus de existentie-uniciteitsresultaten van §3.2 gelden automatisch.

### Stelling 3.3.4 — Unitaire ruimte is genormeerd

> ⚠️ **Belangrijk voor examen:** Stelling + bewijs van driehoeksongelijkheid via Cauchy-Schwarz.

**Stelling:** In een unitaire ruimte definieert

$$
\|x\| = \sqrt{(x, x)}
$$

een norm, de **geïnduceerde** of **natuurlijke norm**.

### Bewijs

Eigenschappen 1–3 zijn eenvoudig. We bewijzen de driehoeksongelijkheid.

**Voor $\|x+y\| = 0$:** triviaal.

**Voor $\|x+y\| > 0$:** Expansie van $(x+y, x+y)$:

$$
\|x+y\|^2 = (x+y, x+y) = (x, x+y) + (y, x+y).
$$

**Stap 1:** Cauchy-Schwarz op elke term:
$$
|(x, x+y)| \leq \|x\|\, \|x+y\|, \qquad |(y, x+y)| \leq \|y\|\, \|x+y\|.
$$

**Stap 2:** Dus:
$$
\|x+y\|^2 \leq \|x\|\,\|x+y\| + \|y\|\,\|x+y\|.
$$

**Stap 3:** Deel door $\|x+y\|$:
$$
\|x+y\| \leq \|x\| + \|y\|. \quad \square
$$

**Dit zegt ons:** Het scalair product voorziet ons gratis van een norm en via die norm van een afstand. We hebben zo de volledige hiërarchie: scalair product $\Rightarrow$ norm $\Rightarrow$ afstand.

---

### Stelling 3.3.5 — Parallellogramgelijkheid

> ⚠️ **Belangrijk voor examen:** Weten wanneer een norm van een scalair product afkomt.

**Stelling:** Een genormeerde ruimte is unitair (met scalair product dat de norm induceert) **als en slechts als** de norm voldoet aan de **parallellogramgelijkheid**:

$$
\|x + y\|^2 + \|x - y\|^2 = 2(\|x\|^2 + \|y\|^2).
$$

**Bewijs (noodzakelijkheid):**

$$
\|x+y\|^2 + \|x-y\|^2 = (x+y,x+y) + (x-y,x-y)
= 2(x,x) + 2(y,y) = 2\|x\|^2 + 2\|y\|^2. \quad \square
$$

**Meetkundige interpretatie:** In een parallellogram is de som van de kwadraten van de diagonalen gelijk aan de som van de kwadraten van de vier zijden.

**Corollarium — welke normen zijn *geen* scalair product:**

- De $\|\cdot\|_p$-normen voor $p \neq 2$ geven **geen** scalair product (de parallellogramgelijkheid geldt niet).
- De maximumnorm in $C[a,b]$ geeft **geen** scalair product.
- De $L^p[a,b]$-normen voor $p \neq 2$ geven **geen** scalair product.

De $L^2$-norm met gewicht wél: het overeenkomstige scalair product is

$$
(f, g) = \int_a^b w(x)\, f(x)\, \overline{g(x)}\, dx.
$$

---

### Stelling 3.3.6 — Eenheidsbol in unitaire ruimte is strikt convex

> ⚠️ **Belangrijk voor examen:** Dit is de sleutelstelling die uniciteit van beste benadering garandeert in unitaire ruimten.

**Stelling:** De eenheidsbol in een unitaire ruimte is strikt convex.

### Bewijs

We moeten tonen dat voor $x \neq y$ met $\|x\| = \|y\| = 1$: $\|x + y\| < 2$.

**Stap 1:** Pas de parallellogramgelijkheid toe:
$$
\|x+y\|^2 = -\|x-y\|^2 + 2(\|x\|^2 + \|y\|^2) = -\|x-y\|^2 + 4.
$$

**Stap 2:** Omdat $x \neq y$, geldt $\|x - y\| > 0$, dus:
$$
\|x+y\|^2 < 4 \implies \|x+y\| < 2. \quad \square
$$

**Dit zegt ons:** In combinatie met Stelling 3.2.3 weten we nu dat in elke unitaire ruimte, voor elke eindigdimensionale deelruimte $D$, de beste benadering (m.b.t. de geïnduceerde norm) altijd **bestaat en uniek** is.

---

### Orthogonaliteit en de stelling van Pythagoras

### Definitie

Twee vectoren $x, y \in V$ zijn **orthogonaal** ($x \perp y$) als $(x, y) = 0$.

Een vector $x$ is orthogonaal tot een deelverzameling $D$ ($x \perp D$) als $x \perp d$ voor alle $d \in D$.

### Stelling van Pythagoras (veralgemening)

**Stelling:** Als $x \perp y$ in een unitaire ruimte $V$, dan geldt:

$$
\|x + y\|^2 = \|x\|^2 + \|y\|^2.
$$

**Bewijs:** Triviaal via expansie:
$$
\|x+y\|^2 = (x+y, x+y) = (x,x) + \underbrace{(x,y)}_{=0} + \underbrace{(y,x)}_{=0} + (y,y) = \|x\|^2 + \|y\|^2. \quad \square
$$

---

### Voorbeelden van unitaire ruimten

1. **$\mathbb{R}^n$ of $\mathbb{C}^n$:** scalair product $(x,y) = \sum_i \bar{x}_i y_i$ (of $x^* y$ in matrixnotatie). Geïnduceerde norm = $\|\cdot\|_2$.

2. **$\ell^2$:** kwadratisch sommeerbare rijen, scalair product $(x,y) = \sum_{i=1}^\infty \bar{x}_i y_i$.

3. **$C([a,b], w)$ of $L^2[a,b]$:** scalair product $(f,g) = \int_a^b w(x) f(x)\overline{g(x)}\, dx$.

**Niet-unitaire genormeerde ruimten:** $\ell^p$ voor $p \neq 2$, $L^p[a,b]$ voor $p \neq 2$, $C[a,b]$ met de maximumnorm.

---

## 3.4 Benaderen in een unitaire ruimte

**Context:** We beschouwen beste benaderingen in een *eindigdimensionale deelruimte* van een mogelijks *oneindigdimensionale* unitaire ruimte. Denk: benader een willekeurige continue functie door een veelterm van graad $\leq n$.

Alles wat we in Hoofdstuk 2 deden voor $\mathbb{C}^m$ veralgement hier naar abstracte unitaire ruimten, met scalaire producten in plaats van inwendige producten.

---

### 3.4.1 De grammatrix

### Definitie

De **grammatrix** van een stel vectoren $\{a_1, \ldots, a_n\}$ in een unitaire ruimte is:

$$
G(a_1, \ldots, a_n) = \begin{pmatrix}
(a_1, a_1) & (a_1, a_2) & \cdots & (a_1, a_n) \\
(a_2, a_1) & (a_2, a_2) & \cdots & (a_2, a_n) \\
\vdots & & \ddots & \vdots \\
(a_n, a_1) & (a_n, a_2) & \cdots & (a_n, a_n)
\end{pmatrix}.
$$

**Opmerking:** In $\mathbb{C}^m$ reduceert dit tot $G = A^*A$ (de matrix waarvan de kolommen $a_1, \ldots, a_n$ zijn).

> ⚠️ **Belangrijk voor examen:** Stelling over HPD ↔ lineaire onafhankelijkheid.

### Stelling 3.4.8 — Grammatrix en lineaire onafhankelijkheid

**Stelling:** De grammatrix $G$ is **Hermitisch positief definiet** (HPD) als en slechts als $\{a_1, \ldots, a_n\}$ **lineair onafhankelijk** zijn.

### Bewijs

**$G$ is Hermitisch:** Uit $(a_i, a_j) = \overline{(a_j, a_i)}$ volgt $G = G^*$.

**HPD $\Leftarrow$ lineair onafhankelijk:**

Neem een willekeurige $\xi = [v_1, \ldots, v_n]^T \in \mathbb{C}^n$ met $\xi \neq 0$. Dan:

$$
\xi^* G \xi = \sum_{i,j} \bar{v}_i (a_i, a_j) v_j = \left(\sum_i v_i a_i,\, \sum_j v_j a_j\right) = \left\|\sum_i v_i a_i\right\|^2.
$$

→ We gebruikten de **lineariteit van het scalair product** in beide argumenten om de dubbele som samen te voegen.

Als $\{a_i\}$ lineair onafhankelijk zijn en $\xi \neq 0$, dan is $\sum_i v_i a_i \neq 0$, dus $\xi^* G \xi > 0$.

**Niet HPD $\Leftarrow$ lineair afhankelijk:**

Als $\{a_i\}$ lineair afhankelijk zijn, bestaat er $\hat{\xi} \neq 0$ met $\sum_i \hat{v}_i a_i = 0$. Dan is $\hat{\xi}^* G \hat{\xi} = 0$, dus $G$ is singulier, niet HPD. $\square$

**Gevolg:** $G$ is inverteerbaar $\iff$ $\{a_1, \ldots, a_n\}$ zijn lineair onafhankelijk.

---

### 3.4.2 Orthogonale projector

### Definitie

Een **projector** $P: V \to V$ is een lineaire operator die **idempotent** is: $P(Pv) = Pv$.

Elke vector $v$ ontbindt uniek als $v = Pv + (I-P)v$ met $Pv \in R(P)$ en $(I-P)v \in N(P)$.

Een projector is **orthogonaal** als $R(P) \perp N(P)$.

### Eigenschap 3.1 — Karakterisering van orthogonale projector

> ⚠️ **Belangrijk voor examen:** Karakterisering via zelfgeadjungeerdheid.

**Stelling:** $P$ is een orthogonale projector $\iff$ $(Pw, v) = (w, Pv)$ voor alle $v, w \in V$.

**Bewijs (noodzakelijk):** Als $R(P) \perp N(P)$, dan voor alle $v, w$:
$$
(Pw, (I-P)v) = 0, \quad ((I-P)w, Pv) = 0.
$$
Uitschrijven: $(Pw, v) - (Pw, Pv) = 0$ en $(w, Pv) - (Pw, Pv) = 0$. Dus $(Pw, v) = (Pw, Pv) = (w, Pv)$.

**Bewijs (voldoende):** Neem $x = Pu \in R(P)$ en $y \in N(P)$. Dan:
$$
(x, y) = (Pu, y) = (u, Py) = (u, 0) = 0. \quad \square
$$

---

### Constructie van de orthogonale projector op $D$

Zij $D = \text{span}\{a_1, \ldots, a_n\}$ een deelruimte van $V$.

**Doel:** bereken $P_D v$ voor een gegeven $v \in V$.

**Aanpak:** we willen $v = v_1 + v_2$ met $v_1 \in D$ en $v_2 \perp D$.

**Stap 1:** Schrijf $v_1 = \sum_{i=1}^n c_i a_i$.

**Stap 2:** De eis $v_2 = v - v_1 \perp D$ geeft voor elke $k$:
$$
(a_k, v - v_1) = 0 \implies \sum_{i=1}^n c_i (a_k, a_i) = (a_k, v), \quad k = 1, \ldots, n.
$$

Dit is het **normaalstelsel**:

$$
G \begin{pmatrix} c_1 \\ \vdots \\ c_n \end{pmatrix} = \begin{pmatrix} (a_1, v) \\ \vdots \\ (a_n, v) \end{pmatrix},
$$

met als oplossing $c = G^{-1} b$ (uniek, want $G$ is inverteerbaar voor een basis).

**Stap 3:** De orthogonale projector is:

$$
P_D v = \sum_{i=1}^n c_i\, a_i.
$$

**Speciale gevallen:**

- **Orthogonale basis** ($a_i \perp a_j$ voor $i \neq j$): $G$ is diagonaal, dus
$$
c_i = \frac{(a_i, v)}{\|a_i\|^2}, \qquad P_D v = \sum_{i=1}^n \frac{(a_i, v)}{\|a_i\|^2} a_i.
$$

- **Orthonormale basis** ($\|a_i\| = 1$ bovenop orthogonaliteit):
$$
\boxed{P_D v = \sum_{i=1}^n (a_i, v)\, a_i.}
$$

> ⚠️ **Belangrijk voor examen:** Dit is de centrale formule voor orthogonale projectie in een unitaire ruimte — veralgemening van formule (2.16) uit Hoofdstuk 2.

---

### Stelling 3.4.9 — Orthogonale projectiestelling

> ⚠️ **Belangrijk voor examen:** Dit is het hoofdresultaat van het hoofdstuk.

**Stelling:** De **beste benadering** van $v \in V$ in de deelruimte $D$ van een unitaire ruimte, m.b.t. de **geïnduceerde norm**, wordt gegeven door

$$
\hat{y} = P_D v.
$$

**Intuïtie:** De loodrechte projectie van $v$ op $D$ is de dichtstbijzijnde punt in $D$. Het residu $v - \hat{y}$ is loodrecht op $D$ — precies de orthogonaliteitsconditie.

Dit is een directe veralgemening van Stelling 2.4.2 en 2.4.1 uit Hoofdstuk 2.

---

### 3.4.3 Orthogonalisatieprocedures: Gram-Schmidt in unitaire ruimten

**Vertrekpunt:** Een basis $\{a_1, \ldots, a_n\}$ van een deelruimte $D \subset V$.

**Doel:** Vervang door een orthonormale basis $\{q_1, \ldots, q_n\}$ zodat $\text{span}\{q_1, \ldots, q_k\} = \text{span}\{a_1, \ldots, a_k\}$ voor alle $k$ (geneste eigenschap).

**Structuur:** De relatie tussen de twee basissen is

$$
a_j = r_{1j} q_1 + r_{2j} q_2 + \cdots + r_{jj} q_j, \qquad j = 1, \ldots, n.
$$

Dit is de veralgemening van de (onvolledige) QR-factorisatie $A = \hat{Q}\hat{R}$.

### Algoritme 5 — Klassiek Gram-Schmidt in unitaire ruimte

```
voor j = 1 tot n:
    v_j = a_j
    voor i = 1 tot j-1:
        r_ij = (q_i, a_j)          ← scalair product i.p.v. inwendig product
        v_j = v_j - r_ij * q_i
    r_jj = ||v_j||                 ← geïnduceerde norm
    q_j = v_j / r_jj
```

> ⚠️ **Belangrijk voor examen:** De analogie met Algoritme 1 uit Hoofdstuk 2: inwendige producten $\to$ scalaire producten, Euclidische norm $\to$ geïnduceerde norm.

**Meetkundige interpretatie van stap $j$:**

$$
v_j = a_j - P_{\langle q_1, \ldots, q_{j-1} \rangle} a_j.
$$

We trekken de orthogonale projectie van $a_j$ op de reeds gebouwde deelruimte $\langle q_1, \ldots, q_{j-1} \rangle$ af. Het resultaat $v_j$ is orthogonaal op alle voorgaande $q_i$. Na normalisatie geeft dit $q_j$.

**Opmerking:** Gewijzigd Gram-Schmidt (Algoritme 2) en herorthogonalisatie blijven ook geldig in de abstracte setting.

**Beperking:** De orthogonalisatiemethoden gebaseerd op *unitaire transformaties* (Householder, Givens) uit §2.3.4 veralgemenen **niet** naar abstracte unitaire ruimten. De reden: die methoden werken op *rijen* van de matrix, wat geen analoog heeft in de abstracte setting waar de "kolommen" functies zijn.

---

### Voorbeeld 3.2 — Gram-Schmidt in $L^2[-1,1]$

**Ruimte:** $L^2[-1,1]$ met scalair product $(f,g) = \int_{-1}^1 f(x) g(x)\, dx$.

**Startbasis:** $\{1, x, x^2, x^3\}$ voor $P_3[-1,1]$.

We berekenen een orthonormale basis $\{P_0, P_1, P_2, P_3\}$.

**Iteratie 1:**

$$
r_{11} = \left(\int_{-1}^1 1\, dt\right)^{1/2} = \sqrt{2}, \qquad P_0(x) = \frac{1}{\sqrt{2}}.
$$

**Iteratie 2:**

$$
r_{12} = (P_0, x) = \frac{1}{\sqrt{2}}\int_{-1}^1 t\, dt = 0 \quad \text{(oneven integraal)}.
$$

$$
v_2 = x - 0 = x, \qquad r_{22} = \left(\int_{-1}^1 t^2\, dt\right)^{1/2} = \sqrt{\frac{2}{3}}, \qquad P_1(x) = \sqrt{\frac{3}{2}}\, x.
$$

**Iteratie 3:**

$$
r_{13} = (P_0, x^2) = \frac{1}{\sqrt{2}}\int_{-1}^1 t^2\, dt = \frac{1}{\sqrt{2}} \cdot \frac{2}{3} = \frac{\sqrt{2}}{3}.
$$

$$
r_{23} = (P_1, x^2) = \sqrt{\frac{3}{2}} \int_{-1}^1 t^3\, dt = 0 \quad \text{(oneven)}.
$$

$$
v_3 = x^2 - \frac{\sqrt{2}}{3} \cdot \frac{1}{\sqrt{2}} = x^2 - \frac{1}{3}.
$$

$$
r_{33} = \left(\int_{-1}^1 \left(t^2 - \frac{1}{3}\right)^2 dt\right)^{1/2} = \sqrt{\frac{8}{45}}, \qquad P_2(x) = \sqrt{\frac{5}{2}} \cdot \frac{1}{\sqrt{2}}\left(3x^2 - 1\right) = \sqrt{\frac{5}{8}}(3x^2 - 1).
$$

**Iteratie 4:** Analoog, resultaat:

$$
P_3(x) = \sqrt{\frac{7}{2}} \cdot \frac{\sqrt{5}}{2}\left(x^3 - \frac{3}{5}x\right).
$$

**Opmerking:** De functies $P_0, P_1, P_2, P_3$ zijn genormeerde versies van de **Legendre-polynomen** — een fundamentele familie in wiskundige analyse en numerieke integratie.

---

### 3.4.4 Beste benadering: aanbevolen aanpak

Gegeven een niet-orthogonale basis $\{a_1, \ldots, a_n\}$, is de aangewezen aanpak voor het berekenen van de beste benadering $\hat{y}$ van $v$ in $D$:

**Stap 1:** Bereken een orthonormale basis $\{q_1, \ldots, q_n\}$ via Gram-Schmidt (Algoritme 5).

**Stap 2:** Bereken de coëfficiënten in de orthonormale basis:
$$
d_i = (q_i, v), \qquad i = 1, \ldots, n.
$$

**Stap 3:** De beste benadering is:
$$
\hat{y} = \sum_{i=1}^n d_i q_i.
$$

**Stap 4 (optioneel):** Als men de coëfficiënten $c = [c_1, \ldots, c_n]^T$ in de oorspronkelijke basis $\{a_i\}$ wil (zodat $\hat{y} = \sum_i c_i a_i$), los dan het driehoeksstelsel $\hat{R} c = d$ op, met $d = [d_1, \ldots, d_n]^T$ en $\hat{R}$ de bovendriehoeksmatrix van Gram-Schmidt-coëfficiënten.

**Waarschuwing:** Het direct oplossen van het normaalstelsel $Gc = b$ is numeriek gevaarlijk door het hoge conditiegetal van de grammatrix $G$ (zie Tabel 3.1: voor monomiale basis en $n=18$ is $\kappa(G) \approx 10^{13}$).

### Voorbeeld 3.4

In $L^2[-1,1]$ met orthonormale basis $\{P_0, P_1, P_2, P_3\}$ uit Voorbeeld 3.2 is de beste benadering van $f$ door een veelterm van graad $\leq 3$:

$$
g(x) = \sum_{i=1}^4 d_i\, P_{i-1}(x), \qquad d_i = (P_{i-1}, f) = \int_{-1}^1 P_{i-1}(x) f(x)\, dx.
$$

In de monomiale basis: $g(x) = \sum_{i=1}^4 c_i x^{i-1}$ met $c = \hat{R}^{-1} d$.

---

## Samenvatting en verbanden

| Structuur | Wat extra | Afstand | Norm | Scalair product |
|---|---|---|---|---|
| Metrische ruimte | alleen afstand | ✓ | — | — |
| Genormeerde ruimte | norm → afstand | ✓ | ✓ | — |
| Unitaire ruimte | scalair product → norm → afstand | ✓ | ✓ | ✓ |

**Hiërarchie van resultaten:**

- In een *metrische ruimte*: beste benadering bestaat niet altijd, niet altijd uniek.
- In een *strikt genormeerde ruimte* + eindigdim. deelruimte: **bestaat en is uniek** (Stelling 3.2.3).
- In een *unitaire ruimte*: beste benadering = **orthogonale projectie** $P_D v$, berekend via het normaalstelsel of Gram-Schmidt (Stelling 3.4.9).

**Veralgemening van Hoofdstuk 2:**

| Hoofdstuk 2 (in $\mathbb{C}^m$) | Hoofdstuk 3 (in unitaire ruimte) |
|---|---|
| Inwendig product $x^* y$ | Scalair product $(x, y)$ |
| Euclidische norm $\|\cdot\|_2$ | Geïnduceerde norm $\|\cdot\| = \sqrt{(\cdot,\cdot)}$ |
| Grammatrix $G = A^*A$ | Grammatrix $G_{ij} = (a_i, a_j)$ |
| Projectie $P_D v = A(A^*A)^{-1}A^*v$ | Projectie via normaalstelsel (3.47)–(3.49) |
| QR-factorisatie (Algoritme 1) | Gram-Schmidt in unitaire ruimte (Algoritme 5) |
