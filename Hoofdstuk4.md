# Deel 4 – Benadering door middel van veeltermen

## Overzicht en rode draad

In Hoofdstuk 3 legden we de abstracte theorie vast: metrische ruimten, genormeerde ruimten, unitaire ruimten en orthogonale projectie. Nu passen we die theorie toe op een concrete klasse van benaderingen: **veeltermen**.

Veeltermen zijn ideale benaderingsfuncties — ze zijn gemakkelijk te evalueren, te differentiëren, te integreren. De centrale vraag is: gegeven een functie $f$, wat is de "beste" benadering door een veelterm van graad $\leq n$?

De rode draad:
1. Definieer een goed benaderingscriterium → **kleinste-kwadraten** in een gewogen $L^2$-ruimte.
2. Gebruik een orthogonale veeltermbasis om de berekeningen efficiënt en numeriek stabiel te maken.
3. Bestudeer de eigenschappen van orthogonale veeltermen (recursiebetrekking, nulpunten).
4. Bekijk de klassieke families: **Legendre** en **Chebyshev**.

---

## 4.1 Inleiding

### Het kleinste-kwadratencriterium

Een goed benaderingscriterium moet aan vier eisen voldoen:

1. **Praktische betekenis:** de fout moet meetbaar zijn en een zinvolle interpretatie hebben.
2. **Ordening:** men moet twee benaderingen met elkaar kunnen vergelijken — welke is beter?
3. **Realistische eis:** de beste benadering moet bereikbaar zijn; je kunt niet eisen dat de fout overal tegelijk minimaal is.
4. **Berekenbaarheid:** er moet een algoritme bestaan om de beste benadering te vinden.

Het **kleinste-kwadratencriterium** voldoet aan alle vier:

$$
\int_a^b w(x)\, |f(x) - y(x)|^2\, dx \quad \text{minimaliseren.}
$$

Dit is de natuurlijke afstand in een gewogen $L^2$-ruimte met scalair product

$$
(f, g) = \int_a^b w(x)\, f(x)\, g(x)\, dx,
$$

met $w(x) > 0$ een gewichtsfunctie en $\int_a^b w(x)\, dx < \infty$.

**Bijzonder voordeel:** omdat dit criterium uit een scalair product komt, zijn existentie en uniciteit van de beste benadering gegarandeerd (Stelling 3.2.3 + Stelling 3.3.6).

**Momenten van de gewichtsfunctie:** Om te garanderen dat alle veeltermen tot de gewogen $L^2$-ruimte behoren, eisen we dat de **momenten** eindig zijn:

$$
m_k = \int_a^b w(x)\, x^k\, dx < \infty \quad \text{voor } k = 0, 1, 2, \ldots
$$

---

### De monomiale basis: een numerieke valkuil

De meest voor de hand liggende keuze voor een basis van $P_n[a,b]$ is de **monomiale basis**:

$$
\phi_k(x) = x^k, \quad k = 0, \ldots, n.
$$

De coëfficiënten van de beste benadering $y_n(x) = \sum_{k=0}^n a_k x^k$ voldoen dan aan het normaalstelsel $Gc = b$ met grammatrix

$$
G_{jk} = (\phi_j, \phi_k) = \int_a^b w(x)\, x^{j+k}\, dx.
$$

> ⚠️ **Belangrijk voor examen:** De grammatrix voor de monomiale basis op $[0,1]$ met $w \equiv 1$ is de **Hilbert-matrix**:
> $$G_{jk} = \frac{1}{j+k+1}.$$
> Dit is een berucht slecht geconditioneerde matrix.

### Voorbeeld 4.1 — Hilbert-matrix

Voor $n = 4$: $\kappa(G) \approx 47661$. Voor $n = 10$: $\kappa(G) \approx 10^{13}$.

Dit betekent dat men bij het oplossen van het normaalstelsel tot 13 decimale cijfers aan nauwkeurigheid kan verliezen. De inverse van de $5 \times 5$ Hilbert-matrix bevat al waarden tot $179200$ — een teken dat de matrix extreem gevoelig is voor afrondingsfouten.

**Conclusie:** Gebruik de monomiale basis nooit rechtstreeks voor het oplossen van het normaalstelsel. Gebruik altijd een orthogonale basis.

---

### De orthogonale veeltermbasis: de correcte aanpak

Als we een orthogonale basis $\{\phi_0, \phi_1, \ldots, \phi_n\}$ gebruiken, wordt de grammatrix **diagonaal** en reduceert de projectieformule tot (zie §3.4):

$$
\boxed{y_n(x) = \sum_{k=0}^n a_k \phi_k(x) \quad \text{met} \quad a_k = \frac{(\phi_k, f)}{(\phi_k, \phi_k)} = \frac{\int_a^b w(x)\, f(x)\, \phi_k(x)\, dx}{\int_a^b w(x)\, \phi_k^2(x)\, dx}.}
$$

> ⚠️ **Belangrijk voor examen:** Dit is de centrale formule voor de kleinste-kwadratenbenadering in een orthogonale basis. Kennen en kunnen toepassen.

**Normalisatie-vrijheid:** Wanneer $\{\phi_k\}$ orthogonaal zijn, kan men de schaalfactor van elke $\phi_k$ vrij kiezen. Gebruikelijke keuzen:
- $\|\phi_k\| = 1$: orthonormale basis (Gram-Schmidt resultaat).
- $\phi_k(b) = 1$: normalisatie in een eindpunt.
- **Monische** normalisatie: de coëfficiënt van $x^k$ in $\phi_k(x)$ is 1, dus $\phi_k(x) = x^k + \text{lagere machten}$.

---

## 4.2 Eigenschappen van orthogonale veeltermen

We beschouwen een rij $\phi_0(x), \phi_1(x), \phi_2(x), \ldots$ waarbij $\phi_k$ een orthogonale veelterm van **strikte graad** $k$ is.

### Stelling 4.2.1 — Basis voor $P_n[a,b]$

**Stelling:** De verzameling $\{\phi_0, \phi_1, \ldots, \phi_n\}$ is een basis van $P_n[a,b]$.

**Bewijs:** De grammatrix van een orthogonaal stel is inverteerbaar (Stelling 3.4.8), dus de $n+1$ veeltermen zijn lineair onafhankelijk in een ruimte van dimensie $n+1$. Bijgevolg vormen ze een basis. $\square$

### Stelling 4.2.2 — Orthogonaliteit tot lagere-graadsveeltermen

**Stelling:** Een veelterm $\phi_k$ uit een orthogonale rij is orthogonaal tot **elke** veelterm van graad $< k$.

**Bewijs:** Elke veelterm van graad $< k$ is een lineaire combinatie van $\phi_0, \ldots, \phi_{k-1}$, en $\phi_k$ is orthogonaal tot elk van die. $\square$

**Concreet:** $\phi_k(x) \perp 1, x, x^2, \ldots, x^{k-1}$.

---

### Stelling 4.2.3 — Drietermsrecursiebetrekking

> ⚠️ **Belangrijk voor examen:** Dit is een van de belangrijkste stellingen van het hoofdstuk. Bewijs kennen.

**Stelling:** Elke rij orthogonale veeltermen van opeenvolgende graad voldoet aan een **drietermsrecursiebetrekking**:

$$
\phi_0(x) = \lambda_0
$$
$$
\phi_1(x) = \lambda_1 \left(x - \alpha_1\right) \phi_0(x)
$$
$$
\phi_k(x) = \lambda_k \left[(x - \alpha_k)\, \phi_{k-1}(x) - \beta_k\, \phi_{k-2}(x)\right] \quad \text{voor } k \geq 2,
$$

met

$$
\alpha_k = \frac{(x\phi_{k-1},\, \phi_{k-1})}{(\phi_{k-1},\, \phi_{k-1})}, \qquad \beta_k = \frac{(x\phi_{k-1},\, \phi_{k-2})}{(\phi_{k-2},\, \phi_{k-2})},
$$

en $\lambda_k$ bepaald door de normalisatievoorwaarde.

**Intuïtie:** In plaats van Gram-Schmidt opnieuw te berekenen (wat $O(n^2)$ scalaire producten kost), kunnen we elke nieuwe orthogonale veelterm berekenen uit enkel de twee vorige. Dit kost slechts $O(n)$ scalaire producten — een enorme besparing.

### Bewijs

**Stap 1:** De veelterm $x\phi_{k-1}(x)$ heeft graad $k$ en kan worden ontbonden in de basis $\{\phi_0, \ldots, \phi_k\}$:

$$
x\phi_{k-1}(x) = b_0 \phi_0(x) + b_1 \phi_1(x) + \cdots + b_k \phi_k(x),
$$

met coëfficiënten (via formule 4.7):

$$
b_l = \frac{(x\phi_{k-1},\, \phi_l)}{(\phi_l,\, \phi_l)} = \frac{(\phi_{k-1},\, x\phi_l)}{(\phi_l,\, \phi_l)} \quad \text{voor } l \leq k.
$$

→ We gebruiken dat het scalair product symmetrisch is: $(x\phi_{k-1}, \phi_l) = (\phi_{k-1}, x\phi_l)$.

**Stap 2:** De veelterm $x\phi_l(x)$ heeft graad $l + 1$. Door Stelling 4.2.2 geldt dat $b_l = 0$ wanneer $k - 1 > l + 1$, d.w.z. wanneer $l < k - 2$.

→ **Sleutelinzicht:** $\phi_{k-1}$ is orthogonaal tot alle veeltermen van graad $< k-1$. Maar $x\phi_l$ heeft graad $l+1$. Als $l+1 < k-1$, d.w.z. $l < k-2$, dan is $x\phi_l$ van te lage graad om een niet-nul inproduct te geven.

**Stap 3:** Er blijven slechts drie termen over:

$$
x\phi_{k-1}(x) = b_{k-2}\phi_{k-2}(x) + b_{k-1}\phi_{k-1}(x) + b_k\phi_k(x).
$$

**Stap 4:** Herschrijf naar $\phi_k$:

$$
\phi_k(x) = \frac{1}{b_k}\left[(x - b_{k-1})\phi_{k-1}(x) - b_{k-2}\phi_{k-2}(x)\right].
$$

Identificeer $\alpha_k = b_{k-1}$, $\beta_k = b_{k-2}/b_{k-1,\text{norm}}$ en $\lambda_k = 1/b_k$ overeenkomstig de normalisatiekeuze. $\square$

**Dit zegt ons:** De volledige rij $\phi_0, \phi_1, \ldots, \phi_n$ kan worden opgebouwd met $O(n)$ scalaire producten in plaats van $O(n^2)$, doordat elke nieuwe veelterm enkel afhangt van de twee vorige.

**Expliciete coëfficiënten** (continue geval):

$$
\alpha_k = \frac{\int_a^b w(x)\, x\, \phi_{k-1}^2(x)\, dx}{\int_a^b w(x)\, \phi_{k-1}^2(x)\, dx}, \qquad \beta_k = \frac{\int_a^b w(x)\, x\, \phi_{k-1}(x)\, \phi_{k-2}(x)\, dx}{\int_a^b w(x)\, \phi_{k-2}^2(x)\, dx}.
$$

---

### Voorbeeld 4.2 — Monische orthogonale veeltermen op $[-1,1]$

We bouwen de monische orthogonale veeltermen $P_k$ op voor $w(x) \equiv 1$ op $[-1,1]$.

Door de **symmetrie** van het interval: $\alpha_k = 0$ voor alle $k$ (de integralen over oneven functies zijn nul).

**Startwaarden:** $P_0(x) = 1$, $P_1(x) = x$.

**$k = 2$:**

$$
\beta_2 = \frac{\int_{-1}^1 x \cdot P_1(x) \cdot P_0(x)\, dx}{\int_{-1}^1 P_0^2(x)\, dx} = \frac{\int_{-1}^1 x^2\, dx}{\int_{-1}^1 1\, dx} = \frac{2/3}{2} = \frac{1}{3}.
$$

$$
P_2(x) = x \cdot P_1(x) - \frac{1}{3} P_0(x) = x^2 - \frac{1}{3}.
$$

**$k = 3$:**

$$
\beta_3 = \frac{\int_{-1}^1 x \cdot P_2(x) \cdot P_1(x)\, dx}{\int_{-1}^1 P_1^2(x)\, dx} = \frac{\int_{-1}^1 x^2(x^2 - 1/3)\, dx}{2/3} = \frac{4/15}{2/3} = \frac{2}{5}.
$$

$$
P_3(x) = x^3 - \frac{3}{5}x.
$$

**$k = 4$:** $\beta_4 = 9/35$, geeft $P_4(x) = x^4 - \frac{6}{7}x^2 + \frac{3}{35}$.

Dit zijn, op evenredigheidsfactoren na, de **Legendre-veeltermen**.

---

### Stelling 4.2.4 — Nulpunten van orthogonale veeltermen

> ⚠️ **Belangrijk voor examen:** Stelling en bewijs kennen.

**Stelling:** De veelterm $\phi_k(x)$ heeft precies **$k$ enkelvoudige reële nulpunten**, alle gelegen in het **open interval** $(a, b)$.

**Intuïtie:** Orthogonale veeltermen oscilleren sterk in $[a,b]$. Een veelterm van graad $k$ die orthogonaal is tot alle constanten moet $k$ keer van teken wisselen in $(a,b)$.

### Bewijs

Voor $k = 0$: $\phi_0$ is een constante, geen nulpunten — triviaal.

Voor $k > 0$: Veronderstel dat $\phi_k(x)$ slechts $m < k$ tekenwisselingen heeft in $(a,b)$, in de punten $x_1 < x_2 < \cdots < x_m$.

**Stap 1:** Definieer $q(x) = (x - x_1)(x - x_2) \cdots (x - x_m)$. Dit is een veelterm van graad $m < k$.

**Stap 2:** De functie $\phi_k(x) \cdot q(x)$ wisselt **nooit** van teken in $[a,b]$: elke tekenwissel van $\phi_k$ wordt gecompenseerd door de factor $(x - x_i)$ die ook van teken wisselt.

**Stap 3:** Dus $\phi_k(x) q(x)$ is overal positief of overal negatief op $[a,b]$, en:

$$
\int_a^b w(x)\, \phi_k(x)\, q(x)\, dx \neq 0.
$$

**Stap 4:** Maar $q(x)$ is een veelterm van graad $m < k$, en $\phi_k \perp$ alle veeltermen van lagere graad (Stelling 4.2.2). Dus de integraal **moet** nul zijn — een contradictie.

**Conclusie:** $m < k$ is onmogelijk, dus $\phi_k$ heeft minstens $k$ tekenwisselingen en dus minstens $k$ reële nulpunten in $(a,b)$. Omdat $\phi_k$ van graad $k$ is, zijn het er precies $k$, alle enkelvoudig. $\square$

**Gevolgen:**
- De randpunten $a$ en $b$ zijn **nooit** nulpunten van $\phi_k$.
- De normalisatievoorwaarden $\phi_k(a) = 1$ of $\phi_k(b) = 1$ zijn dus altijd zinvol.
- Buiten $[a,b]$ is $\phi_k$ monotoon (geen verdere nulpunten).

---

### Stelling 4.2.5 — Nulpunten als eigenwaarden van een tridiagonale matrix

> ⚠️ **Belangrijk voor examen:** Dit is de praktische methode om de nulpunten te berekenen.

**Stelling:** De $n$ nulpunten van $\phi_n(x)$ zijn de **eigenwaarden** van de $n \times n$ tridiagonale matrix

$$
A = \text{tridiag}(\beta_{k-1}; \alpha_k; \tilde{\beta}_k), \quad \tilde{\beta}_k = 1/\lambda_k,
$$

waarbij de exacte vorm afhangt van de normalisatiekeuze.

### Bewijs

**Stap 1:** Herschrijf de drietermsrecursiebetrekking:

$$
\frac{1}{\lambda_k}\phi_{k-2}(x) + \alpha_k \phi_{k-1}(x) + \tilde{\beta}_k \phi_k(x) = x\phi_{k-1}(x).
$$

**Stap 2:** Definieer de vector $\Phi(x) = [\phi_0(x), \phi_1(x), \ldots, \phi_{n-1}(x)]^T$. Dan schrijft het stelsel van recursiebetrekkingen (voor $k = 1, \ldots, n-1$) zich als:

$$
A\Phi(x) = x\Phi(x) - \tilde{\beta}_n \phi_n(x)\, e_n,
$$

met $e_n$ de $n$-de standaardbasisvector.

**Stap 3:** Evalueer in een nulpunt $x_k$ van $\phi_n$: dan is $\phi_n(x_k) = 0$, en:

$$
A\Phi(x_k) = x_k \Phi(x_k).
$$

Dus $x_k$ is een eigenwaarde van $A$ met eigenvector $\Phi(x_k)$.

**Stap 4:** Er zijn $n$ nulpunten en $A$ heeft hoogstens $n$ eigenwaarden. Alle $n$ nulpunten zijn eigenwaarden en omgekeerd. $\square$

**Praktisch belang:** Dit verbindt orthogonale veeltermen met eigenwaardeproblemen. De nulpunten — die we nodig hebben voor Gauss-kwadratuur en optimale interpolatie — kunnen efficiënt berekend worden via eigenwaardealgoritmen voor tridiagonale matrices.

---

## 4.3 Continue kleinste-kwadratenbenadering

### Foutenkromme en haar nulpunten

De kleinste-kwadratenbenadering van graad $n$ in een orthogonale basis is:

$$
y_n(x) = \sum_{k=0}^n a_k \phi_k(x), \qquad a_k = \frac{(\phi_k, f)}{(\phi_k, \phi_k)}.
$$

Het **residu** (fout) is:

$$
r_n(x) = f(x) - y_n(x).
$$

**Heuristiek voor de fout:** Onder milde voorwaarden dalen de coëfficiënten $|a_k|$ snel naar nul. Een schatting van het residu is:

$$
r_n(x) \approx -a_{n+1}\phi_{n+1}(x).
$$

Omdat $\phi_{n+1}$ precies $n+1$ nulpunten heeft in $(a,b)$, verwachten we dat ook het residu nulpunten heeft nabij die plaatsen.

---

### Stelling 4.3.6 — Interpolerende eigenschap

> ⚠️ **Belangrijk voor examen:** Stelling en bewijs kennen.

**Stelling:** Zij $f$ een continue functie op $[a,b]$. Dan wordt de fout $r_n(x) = f(x) - y_n(x)$ van de kleinste-kwadratenbenadering van graad $n$ **nul in minstens $n+1$ punten** van het open interval $(a,b)$.

**Intuïtie:** Het residu staat loodrecht op de volledige benaderingsruimte $P_n[a,b]$. Als het residu maar $m \leq n$ tekenwisselingen zou hebben, dan zou het inproduct met een geschikte veelterm van graad $m \leq n$ niet-nul zijn — tegenspraak.

### Bewijs

Veronderstel dat $r_n(x)$ slechts $m \leq n$ tekenwisselingen heeft in $(a,b)$, in de punten $x_1 < \cdots < x_m$.

**Stap 1:** Definieer $V_m(x) = (x - x_1)(x - x_2) \cdots (x - x_m)$, een veelterm van graad $m \leq n$.

**Stap 2:** Het product $r_n(x) \cdot V_m(x)$ wisselt nergens van teken in $[a,b]$ (de tekenwisselingen van $r_n$ worden gecompenseerd). Dus:

$$
\int_a^b w(x)\, r_n(x)\, V_m(x)\, dx \neq 0. \tag{$*$}
$$

**Stap 3:** Maar het residu van een orthogonale projectie staat **loodrecht op de benaderingsruimte** (zie §3.4): $(r_n, p) = 0$ voor alle $p \in P_n[a,b]$.

Omdat $\deg(V_m) = m \leq n$, geldt $V_m \in P_n[a,b]$, dus:

$$
\int_a^b w(x)\, r_n(x)\, V_m(x)\, dx = (r_n, V_m) = 0.
$$

Dit is in tegenspraak met $(*)$. Dus $m \leq n$ is onmogelijk: $r_n$ heeft minstens $n+1$ tekenwisselingen en dus minstens $n+1$ nulpunten in $(a,b)$ (uit continuïteit van $f$). $\square$

**Dit zegt ons:** De kleinste-kwadratenbenadering interpoleert de functie $f$ automatisch in minstens $n+1$ punten — zonder dat we dat expliciet hebben opgelegd. De benaderingspunten liggen in de buurt van de nulpunten van $\phi_{n+1}(x)$.

---

### Evaluatie: het rekenschema van Smith (Clenshaw)

Gegeven de coëfficiënten $a_0, \ldots, a_n$ en de recursieparameters $\alpha_k, \beta_k, \lambda_k$, kan $y_n(x)$ efficiënt worden geëvalueerd zonder expliciet de veeltermen te berekenen:

$$
\begin{cases}
b_{n+2} = b_{n+1} = 0 \\
b_k = a_k + \lambda_{k+1}(x - \alpha_{k+1})\, b_{k+1} - \lambda_{k+2}\beta_{k+2}\, b_{k+2}, \quad k = n, n-1, \ldots, 0 \\
y_n(x) = b_0\, \phi_0(x) = b_0 \lambda_0.
\end{cases}
$$

Dit is het **rekenschema van Smith** (ook Clenshaw-algoritme genoemd).

**Correctheid:** De recursiebetrekking voor $b_k$ is algebraïsch equivalent aan de optelling $y_n(x) = \sum_k a_k \phi_k(x)$, dankzij de drietermsrecursiebetrekking voor $\phi_k$.

---

### Interpolatie in de nulpunten van $\phi_{n+1}$

Uit de foutschatting $r_n \approx -a_{n+1}\phi_{n+1}$ volgt dat de nulpunten van $\phi_{n+1}$ goede interpolatiepunten zijn. **Experiment (Voorbeeld 4.5):** De veelterm van graad 4 die $e^x$ interpoleert in de 5 nulpunten van de Legendre-veelterm van graad 5 is bijna identiek aan de echte kleinste-kwadratenbenadering van graad 4:

| Methode | Coëff. van $x^0$ | Coëff. van $x^4$ |
|---|---|---|
| Kleinste kwadraten (4.26) | 1.000031 | 0.043597 |
| Interpolatie in Legendre-punten (4.30) | 1.000000 | 0.043235 |

**Praktisch belang:** Bij equidistante interpolatiepunten kan de fout stijgen bij toenemende graad (het klassieke Runge-fenomeen, duidelijk zichtbaar voor $1/(1+6x^2)$). Interpolatie in de nulpunten van een orthogonale veelterm vermijdt dit probleem.

---

## 4.4 Klassieke orthogonale veeltermen

### Overzicht

De bekende families orthogonale veeltermen zijn:

| Naam | Interval | Gewichtsfunctie $w(x)$ |
|---|---|---|
| **Jacobi** $P_k^{(\alpha,\beta)}(x)$ | $[-1,1]$ | $(1-x)^\alpha(1+x)^\beta$, $\alpha,\beta > -1$ |
| **Legendre** $P_k(x)$ | $[-1,1]$ | $1$ |
| **Chebyshev 1e soort** $T_k(x)$ | $[-1,1]$ | $1/\sqrt{1-x^2}$ |
| **Chebyshev 2e soort** $U_k(x)$ | $[-1,1]$ | $\sqrt{1-x^2}$ |
| **Gegenbauer** $P_k^{(\lambda)}(x)$ | $[-1,1]$ | $(1-x^2)^{\lambda - 1/2}$ |
| **Laguerre** $L_k^{(\alpha)}(x)$ | $[0, \infty)$ | $x^\alpha e^{-x}$, $\alpha > -1$ |
| **Hermite** $H_k(x)$ | $(-\infty, \infty)$ | $e^{-x^2}$ |

---

### 4.4.1 Legendre-veeltermen

### Definitie

De **Legendre-veeltermen** $P_0(x), P_1(x), \ldots$ zijn de orthogonale veeltermen voor het scalair product $(p, q) = \int_{-1}^1 p(x) q(x)\, dx$, genormaliseerd zodat $P_k(1) = 1$.

**Expliciete waarden:**

$$
P_0(x) = 1, \quad P_1(x) = x, \quad P_2(x) = \frac{1}{2}(3x^2 - 1), \quad P_3(x) = \frac{1}{2}(5x^3 - 3x).
$$

**Eigenschappen:**
- $P_k$ heeft $k$ enkelvoudige nulpunten in $(-1, 1)$.
- $P_k$ is een even functie voor even $k$, oneven voor oneven $k$.
- $\|P_k\|^2 = \frac{2}{2k+1}$.

---

### 4.4.2 Chebyshev-veeltermen van de eerste soort

> ⚠️ **Belangrijk voor examen:** Definitie en equi-oscillatorische eigenschap.

### Definitie

De **Chebyshev-veeltermen** $T_0(x), T_1(x), \ldots$ zijn de orthogonale veeltermen voor

$$
(p, q) = \int_{-1}^1 \frac{p(x)\, q(x)}{\sqrt{1-x^2}}\, dx,
$$

genormaliseerd zodat $T_k(1) = 1$.

**Expliciete formule:**

$$
\boxed{T_k(x) = \cos(k \arccos x), \quad x \in [-1, 1].}
$$

**Expliciete waarden en recursiebetrekking:**

$$
T_0(x) = 1, \quad T_1(x) = x, \quad T_2(x) = 2x^2 - 1, \quad T_3(x) = 4x^3 - 3x.
$$

$$
T_{k+1}(x) = 2x\, T_k(x) - T_{k-1}(x).
$$

---

### 4.4.3 Equioscillatiestelling (minimax-benadering)

> ⚠️ **Belangrijk voor examen:** Karakterisering van de minimax-benadering en verband met Chebyshev.

**Probleem:** Soms wil men niet de kwadratische fout minimaliseren, maar de **maximale** fout:

$$
E_n = \max_{x \in [a,b]} |f(x) - y_n(x)| \quad \text{minimaliseren.}
$$

Dit noemt men het **minimaxcriterium**.

### Stelling 4.4.8 — Equioscillatiestelling

**Stelling:** Een veelterm $y_n$ van graad $n$ is een minimax-benadering voor $f$ **als en slechts als** er in $[a,b]$ minstens $n+2$ **extremaalpunten** zijn die afwisselen tussen $+$-punten en $-$-punten:

$$
\text{$+$-punt: } f(x_i) - y_n(x_i) = +E_n, \qquad \text{$-$-punt: } f(x_i) - y_n(x_i) = -E_n.
$$

**Intuïtie:** Als er minder dan $n+2$ alternerende extrema zijn, kan men de veelterm nog verbeteren. De fout moet precies $n+2$ keer zijn maximale waarde aannemen, afwisselend van teken.

**Verband met Chebyshev:** De Chebyshev-veelterm $T_{n+1}$ heeft precies $n+2$ extrema in $[-1,1]$ die afwisselen tussen $+1$ en $-1$. Omdat de fout van de Chebyshev-kleinste-kwadratenbenadering van graad $n$ geschat wordt door $r_n \approx -a_{n+1} T_{n+1}$, voldoet die fout bijna aan de equioscillatiestelling — de Chebyshev-benadering is in de praktijk bijna minimax. Dit verklaart waarom **interpolatie in Chebyshev-punten** zo populair is in numerieke software.

---

## Samenvatting

### De aanpak voor kleinste-kwadratenbenadering

**Gegeven:** functie $f$ op $[a,b]$, gewicht $w(x)$, gewenste graad $n$.

**Stap 1:** Bouw een orthogonale basis $\{\phi_0, \ldots, \phi_n\}$ via de drietermsrecursiebetrekking (§4.2.3), $O(n)$ scalaire producten.

**Stap 2:** Bereken de coëfficiënten $a_k = (\phi_k, f) / (\phi_k, \phi_k)$.

**Stap 3:** De beste benadering is $y_n(x) = \sum_{k=0}^n a_k \phi_k(x)$.

**Stap 4 (evaluatie):** Gebruik het rekenschema van Smith voor efficiënte evaluatie in punten.

### Sleutelresultaten

| Resultaat | Wat het zegt |
|---|---|
| Drietermsrecursiebetrekking (Stelling 4.2.3) | Orthogonale veeltermen efficiënt opbouwen in $O(n)$ |
| Nulpunten zijn reëel en enkelvoudig (Stelling 4.2.4) | Oscillerend gedrag, nulpunten zijn goede interpolatiepunten |
| Nulpunten = eigenwaarden tridiagonale matrix (Stelling 4.2.5) | Numerieke berekening van de nulpunten |
| Interpolerende eigenschap (Stelling 4.3.6) | Kleinste-kwadratenbenadering interpoleert automatisch in $\geq n+1$ punten |
| Equioscillatiestelling (Stelling 4.4.8) | Karakterisering van de minimax-benadering |
| Chebyshev ≈ minimax | Chebyshev-benadering is praktisch optimaal |

### Verband met Hoofdstuk 3

- Het **normaalstelsel** is de toepassing van Stelling 3.4.9 in de functieruimte $L^2[a,b]$.
- De **drietermsrecursiebetrekking** vervangt Gram-Schmidt met een veel efficiënter schema voor de speciale klasse van veeltermen.
- De **orthogonale projector** $P_D f = y_n$ is de beste benadering in de unitaire ruimte.
- De Hilbert-matrix is de grammatrix voor de monomiale basis — het hoge conditiegetal is precies de reden waarom we orthogonale basissen gebruiken (zie §3.4, Tabel 3.1).
