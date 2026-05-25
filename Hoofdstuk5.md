# Deel 5 – Benaderen door middel van splinefuncties

## Overzicht en rode draad

De centrale beperking van veeltermbenadering: één veelterm over het hele interval werkt slecht voor functies die lokaal sterk variëren. Hoge-graads veeltermen oscilleren sterk tussen de gegevens (Runge-fenomeen). De oplossing: splits het interval op en gebruik op elk stuk een **veelterm van lage graad**, aaneengekoppeld met de juiste continuïteit. Zo'n aaneenkoppeling heet een **splinefunctie**.

De rode draad:
1. Definieer splinefuncties en hun vrijheidsgraden.
2. Kies een goede **basis**: de afgeknotte-machtsbasis is eenvoudig maar numeriek onstabiel. De **B-splinebasis** is de correcte keuze.
3. Bestudeer de eigenschappen van B-splines (lokale drager, positiviteit, partitie van eenheid).
4. Ontwikkel efficiënte algoritmen: **de Boor** voor evaluatie, differentiatieschema voor afgeleiden.
5. Uitbreidingen: samenvallende knooppunten en periodieke splines.

---

## 5.1 Motivatie: wanneer falen veeltermen?

Drie situaties waarbij één globale veelterm geen goede benadering geeft:

1. **Functies met wisselend karakter:** sterk oscillerend op één deel, vlak elders. Een globale veelterm moet beide gelijktijdig modelleren, wat hoge graad vereist.
2. **Samengestelde functies:** de te benaderen functie bestaat zelf uit stukken met een scherpe knik (discontinuïteit in een afgeleide). Een globale veelterm heeft die discontinuïteit niet.
3. **Runge-fenomeen:** interpolerende veeltermen van hoge graad met equidistante punten oscilleren sterk — de fout kan zelfs toenemen bij toenemende graad (zie §4.3.3).

**Oplossing:** Verdeel $[a,b]$ in deelintervallen $[t_0, t_1], [t_1, t_2], \ldots, [t_{n-1}, t_n]$ en gebruik op elk stuk een veelterm van lage graad (typisch graad 1, 2 of 3). Koppel de stukken samen met de gewenste continuïteitsgraad.

**Visuele vloeiendheid:** Het menselijk oog kan discontinuïteiten in de **krommimg** (de tweede afgeleide) waarnemen. Om een visueel vloeiende curve te bekomen, moeten de eerste én tweede afgeleiden continu zijn in de knooppunten → **kubische splines** ($k = 3$, $C^2$-continuïteit).

---

## 5.2 Basisbegrippen

### Definitie van een splinefunctie

> ⚠️ **Belangrijk voor examen:** De definitie kennen, inclusief de twee voorwaarden.

### Definitie 5.2.1

Zij een **strikte knooppuntenrij** gegeven: $a = t_0 < t_1 < \cdots < t_{n-1} < t_n = b$.

Een **splinefunctie** $s(x)$ van graad $k > 0$ met knooppunten $t_0, t_1, \ldots, t_n$ is een functie op $[a,b]$ met:

1. Op elk deelinterval $[t_{i-1}, t_i]$ is $s(x)$ een **veelterm van graad $\leq k$**.
2. De functie $s(x)$ en haar afgeleiden $s', s'', \ldots, s^{(k-1)}$ zijn **continu op** $[a,b]$.

**Terminologie:**
- Graad $k = 0$: **stapfunctie** (discontinu in de knooppunten).
- Graad $k = 1$: **gebroken lijn** (continu, niet differentieerbaar in de knooppunten).
- Graad $k = 2$: **kwadratische spline** — parabolen aaneengekoppeld met $C^1$-continuïteit.
- Graad $k = 3$: **kubische spline** — veeltermen van graad 3 met $C^2$-continuïteit. Dit is de meest gebruikte klasse.

Een knooppunt $t_i$ is **actief** als $s^{(k)}$ er een discontinuïteit heeft (de spline is "maximaal ruw" in dat punt). Een veelterm zonder actieve knooppunten is een speciaal geval van een spline.

---

### Speciale klassen van randvoorwaarden

### Definitie 5.2.2 — Natuurlijke splinefunctie

Een spline van **oneven graad** $k = 2m+1$ ($m \geq 1$) is **natuurlijk** als:
$$
s^{(j)}(a) = s^{(j)}(b) = 0 \quad \text{voor } j = m+1, \ldots, 2m.
$$

**Intuïtie:** De hoogste afgeleiden zijn nul aan de randen — de spline wordt "glad" afgeknot. Voor kubische splines ($m=1$): $s''(a) = s''(b) = 0$.

### Definitie 5.2.3 — Periodieke splinefunctie

Een spline is **periodiek** als:
$$
s^{(j)}(a) = s^{(j)}(b) \quad \text{voor } j = 0, 1, \ldots, k-1.
$$

Gebruikt voor gesloten krommen en oppervlakken.

---

### De dimensie van de splineruimte

> ⚠️ **Belangrijk voor examen:** Vrijheidsgraden tellen.

**Telbewijs:**
- $n$ deelintervallen, elk met $k+1$ coëfficiënten → totaal $n(k+1)$ coëfficiënten.
- In elk van de $n-1$ **inwendige** knooppunten worden $k$ continuïteitsvoorwaarden opgelegd ($s, s', \ldots, s^{(k-1)}$ zijn continu).
- Vrijheidsgraden: $n(k+1) - (n-1)k = n + k$.

### Stelling 5.2.1

De vectorruimte van de splinefuncties van graad $k$ met knooppunten $t_0, \ldots, t_n$ heeft **dimensie $n + k$**.

---

### De afgeknotte-machtsbasis

Een eerste basis voor de splineruimte: start met de veeltermbasis op $[t_0, t_1)$, en voeg voor elk volgend interval een correctieterm toe. Dit leidt tot de **afgeknotte-machtsfuncties**:

$$
(x - t_i)_+^k = \begin{cases} 0 & \text{als } x \leq t_i \\ (x - t_i)^k & \text{als } x > t_i \end{cases}
$$

De volledige spline in afgeknotte-machtsvorm:

$$
s(x) = \sum_{i=0}^k a_i x^i + \sum_{i=1}^{n-1} c_i (x - t_i)_+^k.
$$

Dit geeft precies $k+1 + (n-1) = n+k$ coëfficiënten — klopt met de dimensie.

De functies $\{1, x, x^2, \ldots, x^k, (x-t_1)_+^k, \ldots, (x-t_{n-1})_+^k\}$ vormen een basis.

**Probleem:** Deze basis is verre van orthogonaal en geeft aanleiding tot numerieke instabiliteit bij berekeningen. Gebruik in de praktijk de **B-splinebasis**.

---

## 5.3 De B-splinevoorstelling

### 5.3.1 Gedeelde differenties

B-splines worden gedefinieerd via gedeelde differenties van afgeknotte-machtsfuncties.

> ⚠️ **Belangrijk voor examen:** Definitie en eigenschappen van gedeelde differenties.

### Definitie 5.3.4 — Gedeelde differentie

De **gedeelde differentie van orde 0** van $f$ in $x_i$:
$$
f[x_i] = f(x_i).
$$

De **gedeelde differentie van orde $j-i$** in de punten $x_i < x_{i+1} < \cdots < x_j$:
$$
f[x_i, x_{i+1}, \ldots, x_j] = \frac{f[x_{i+1}, \ldots, x_j] - f[x_i, \ldots, x_{j-1}]}{x_j - x_i}.
$$

**Rekenregel (differentietabel):** Elke waarde is het verschil van twee buren gedeeld door hun $x$-afstand. Voorbeeld voor $f[x_1, x_2, x_3, x_4, x_5]$:

$$
\begin{array}{ccccc}
f[x_1] & & & & \\
f[x_2] & f[x_1,x_2] & & & \\
f[x_3] & f[x_2,x_3] & f[x_1,x_2,x_3] & & \\
f[x_4] & f[x_3,x_4] & f[x_2,x_3,x_4] & f[x_1,x_2,x_3,x_4] & \\
f[x_5] & f[x_4,x_5] & f[x_3,x_4,x_5] & f[x_2,x_3,x_4,x_5] & f[x_1,x_2,x_3,x_4,x_5]
\end{array}
$$

**Sleuteleigenschappen:**

1. **Lineariteit:** $(af + bg)[x_i, \ldots, x_j] = a\, f[x_i, \ldots, x_j] + b\, g[x_i, \ldots, x_j]$.

2. **Newton-interpolatie:** De interpolerende veelterm van graad $j-i$ door $(x_i, f_i), \ldots, (x_j, f_j)$ is:
$$
p(x) = f[x_i] + f[x_i,x_{i+1}](x-x_i) + f[x_i,x_{i+1},x_{i+2}](x-x_i)(x-x_{i+1}) + \cdots
$$
de coëfficiënten zijn de gedeelde differenties.

3. **Continuïteit:** $f[x_i, \ldots, x_j]$ is continu in de argumenten als $f$ $(j-i)$ maal continu differentieerbaar is.

4. **Nul voor lage-graadsveeltermen:** $p_m[x_i, \ldots, x_j] = 0$ als $j - i > m$ (veelterm van graad $m$).

5. **Lineaire samenstelling:** $f[x_i, \ldots, x_j] = \sum_{k=i}^j \lambda_k f_k$ voor bepaalde gewichten $\lambda_k$.

6. **Leibniz-formule:** Als $f = g \cdot h$:
$$
f[x_i, \ldots, x_j] = \sum_{r=i}^j g[x_i, \ldots, x_r]\, h[x_r, \ldots, x_j].
$$

---

### 5.3.2 Definitie van B-splines

> ⚠️ **Belangrijk voor examen:** Definitie van de B-spline via gedeelde differentie.

### Definitie 5.3.5

De **gewone B-spline** van graad $k$ met knooppunten $t_i < t_{i+1} < \cdots < t_{i+k+1}$:

$$
M_{i,k+1}(x) = [t_i, t_{i+1}, \ldots, t_{i+k+1}]_t\, (t - x)_+^k.
$$

Dit is de gedeelde differentie van orde $k+1$ naar de variabele $t$ van de afgeknotte-machtsfunctie $(t-x)_+^k$, beschouwd als functie van $t$ (met $x$ als constante parameter).

De **genormaliseerde B-spline**:

$$
N_{i,k+1}(x) = (t_{i+k+1} - t_i)\, M_{i,k+1}(x).
$$

**Hoe lezen:** $M_{i,k+1}$ is geduid met twee indices: $i$ is het startpunt van de knooppuntenrij, $k+1$ is de orde (= graad + 1). De B-spline wordt bepaald door $k+2$ knooppunten $t_i, t_{i+1}, \ldots, t_{i+k+1}$.

**Waarom splinefuncties?** Via eigenschap 5 (lineaire samenstelling) schrijft men:
$$
M_{i,k+1}(x) = \sum_{s=i}^{i+k+1} \lambda_s (t_s - x)_+^k,
$$
wat een lineaire combinatie is van afgeknotte-machtsfuncties — dus een splinefunctie.

---

### Uitbreiding van de knooppuntenrij

Met $n+1$ knooppunten $t_0, \ldots, t_n$ zijn er slechts $n-k$ B-splines beschikbaar ($M_{0,k+1}, \ldots, M_{n-k-1,k+1}$), maar de splineruimte heeft dimensie $n+k$. We hebben $2k$ extra basisfuncties nodig.

**Oplossing:** Voeg $k$ knooppunten **links** van $t_0$ en $k$ knooppunten **rechts** van $t_n$ toe:

$$
t_{-k} \leq \cdots \leq t_{-1} \leq t_0 < t_1 < \cdots < t_{n-1} < t_n \leq t_{n+1} \leq \cdots \leq t_{n+k}.
$$

Dit geeft $n + 2k + 1$ knooppunten en de $n+k$ B-splines:

$$
N_{i,k+1}(x) \quad \text{voor } i = -k, -k+1, \ldots, n-1.
$$

**Veelgebruikte keuze:** de extra knooppunten samenlaten vallen met de randpunten: $t_{-k} = \cdots = t_0 = a$ en $t_n = \cdots = t_{n+k} = b$.

---

### 5.3.3 Eigenschappen van B-splines

> ⚠️ **Belangrijk voor examen:** De properties zijn de sleutel tot het begrip van B-splines. Kennen en begrijpen.

#### Eigenschap 5.1 — Basisgeval ($k = 0$)

$$
M_{i,1}(x) = \frac{1}{t_{i+1} - t_i} \cdot \mathbf{1}_{[t_i, t_{i+1})}(x), \qquad N_{i,1}(x) = \mathbf{1}_{[t_i, t_{i+1})}(x).
$$

De genormaliseerde B-spline van graad 0 is gewoon de karakteristieke functie van het interval $[t_i, t_{i+1})$.

---

#### Eigenschap 5.2 — Lokale drager (support)

> ⚠️ **Belangrijk voor examen:** Dit is de SLEUTELEIGENSCHAP van B-splines.

$$
M_{i,k+1}(x) = 0 \quad \text{voor } x \leq t_i \text{ en voor } x \geq t_{i+k+1}.
$$

**Bewijs (schets):**
- Voor $x \geq t_{i+k+1}$: alle termen $(t_s - x)_+^k = 0$ voor $s \leq i+k+1$, dus $M_{i,k+1}(x) = 0$.
- Voor $x \leq t_i$: het "+"-teken kan worden weggelaten en $(t-x)^k$ is een veelterm van graad $k$ in $t$. De gedeelde differentie van orde $k+1$ van een veelterm van graad $k$ is nul. $\square$

**Belang:** Elke B-spline $N_{i,k+1}$ is nul buiten het interval $[t_i, t_{i+k+1}]$, dat slechts $k+1$ knooppuntenintervallen omspant. In elk punt $x \in [t_j, t_{j+1})$ zijn **hoogstens $k+1$ B-splines** tegelijk niet-nul. Dit geeft de splineberekeningen een **bandsgewijze structuur** — een grote numerieke voordeel.

---

#### Eigenschap 5.3 — Recursiebetrekking

> ⚠️ **Belangrijk voor examen:** Kern van alle B-splineberekeningen.

$$
M_{i,k+1}(x) = \frac{x - t_i}{t_{i+k+1} - t_i}\, M_{i,k}(x) + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_i}\, M_{i+1,k}(x),
$$

$$
\boxed{N_{i,k+1}(x) = \frac{x - t_i}{t_{i+k} - t_i}\, N_{i,k}(x) + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_{i+1}}\, N_{i+1,k}(x).}
$$

**Bewijs (idee):** Schrijf $(t-x)_+^k = (t-x)_+^{k-1} \cdot (t-x)$ en pas de Leibniz-formule toe op dit product. De factor $(t-x)$ is een graad-1 veelterm, zodat de meeste termen in de som verdwijnen. Na vereenvoudiging en gebruik van de recursieve definitie van gedeelde differenties volgt de recursiebetrekking. $\square$

**Interpretatie:** De recursiebetrekking drukt een B-spline van graad $k$ uit als een **convexe combinatie** van twee B-splines van graad $k-1$. De coëfficiënten $\frac{x-t_i}{t_{i+k}-t_i}$ en $\frac{t_{i+k+1}-x}{t_{i+k+1}-t_{i+1}}$ tellen op tot 1 en zijn beide niet-negatief voor $x \in [t_i, t_{i+k+1}]$.

---

#### Eigenschap 5.4 — Positiviteit

$$
M_{i,k+1}(x) > 0 \quad \text{voor } t_i < x < t_{i+k+1} \quad (k \geq 1).
$$

**Belang:** Samen met de partitie-van-eenheid eigenschap (hieronder) garandeert dit dat een spline in B-splinevoorstelling de **convexe romp** van zijn controlepunten niet verlaat.

---

#### Eigenschap 5.5 — Randwaarden

$$
M_{i,k+1}^{(j)}(t_i) = M_{i,k+1}^{(j)}(t_{i+k+1}) = 0 \quad \text{voor } j = 0, \ldots, k-1 \quad (k \geq 1).
$$

De B-spline en haar eerste $k-1$ afgeleiden zijn nul aan de randen van zijn drager. Dit volgt rechtstreeks uit de continuïteitseigenschap van de spline en uit de lokale-drager-eigenschap.

---

#### Eigenschap 5.6 — Partitie van eenheid

> ⚠️ **Belangrijk voor examen:** Centrale eigenschap voor stabiliteit en geometrische interpretatie.

$$
\sum_{i=-k}^{n-1} N_{i,k+1}(x) = 1 \quad \text{voor alle } x \in [t_0, t_n].
$$

**Bewijs (schets):** Voor $x \in [t_j, t_{j+1})$ zijn door de lokale-drager-eigenschap enkel de B-splines $N_{j-k,k+1}, \ldots, N_{j,k+1}$ niet-nul. Via de recursiebetrekking reduceert de som stapsgewijs tot $N_{j,1}(x) = 1$. $\square$

**Belang:** De som van alle B-splines is altijd gelijk aan 1. Dit betekent dat een spline $s(x) = \sum c_i N_{i,k+1}(x)$ een **gewogen gemiddelde** is van de coëfficiënten $c_i$ — het gemiddelde wordt gewogen door de B-splinewaarden. Dit geeft de **convexe-romeigenschap** en maakt de voorstelling geometrisch intuïtief (de spline ligt altijd "tussen" de controlepunten).

---

#### Eigenschap 5.7 — Afgeleide van een genormaliseerde B-spline

$$
N_{i,k+1}'(x) = k \left(\frac{N_{i,k}(x)}{t_{i+k} - t_i} - \frac{N_{i+1,k}(x)}{t_{i+k+1} - t_{i+1}}\right) \quad (k \geq 1).
$$

---

### 5.3.4 Algoritmen voor evaluatie en differentiatie

Een splinefunctie in B-splinevoorstelling:

$$
s(x) = \sum_{i=-k}^{n-1} c_i\, N_{i,k+1}(x).
$$

#### Evaluatie: het algoritme van de Boor

> ⚠️ **Belangrijk voor examen:** Structuur en werking van het algoritme kennen.

### Stelling 5.3.3 (de Boor)

Zij $x \in [t_j, t_{j+1})$. Dan is $s(x) = c_j^{[k]}$, waarbij:

$$
c_i^{[0]} = c_i, \qquad c_i^{[r]} = \alpha_{i,r}\, c_i^{[r-1]} + (1 - \alpha_{i,r})\, c_{i-1}^{[r-1]},
$$

met

$$
\alpha_{i,r} = \frac{x - t_i}{t_{i+k+1-r} - t_i}.
$$

**Bewijs (idee):** Pas de recursiebetrekking voor $N_{i,k+1}$ toe op de som $s(x) = \sum c_i N_{i,k+1}(x)$. De som reduceert stapsgewijs:

$$
s(x) = \sum c_i^{[r]} N_{i,k+1-r}(x) = \cdots = \sum c_i^{[k]} N_{i,1}(x).
$$

Voor $x \in [t_j, t_{j+1})$ is $N_{j,1}(x) = 1$ en alle andere $N_{i,1}(x) = 0$. Dus $s(x) = c_j^{[k]}$. $\square$

**Tabelschema** (voor $x \in [t_j, t_{j+1})$, $k$ ronden):

$$
\begin{array}{ccccc}
c_{j-k}^{[0]} & & & & \\
c_{j-k+1}^{[0]} & c_{j-k+1}^{[1]} & & & \\
\vdots & \vdots & \ddots & & \\
c_{j-1}^{[0]} & c_{j-1}^{[1]} & \cdots & c_{j-1}^{[k-1]} & \\
c_j^{[0]} & c_j^{[1]} & \cdots & c_j^{[k-1]} & c_j^{[k]} = s(x)
\end{array}
$$

**Interpretatie:** Het algoritme van de Boor is de analoog van het Horner-schema voor veeltermen, maar dan voor B-splines. Het is numeriek stabiel doordat alle gewichten $\alpha_{i,r} \in [0,1]$ zijn.

**Lokaal:** Dankzij de lokale-drager-eigenschap zijn enkel de $k+1$ coëfficiënten $c_{j-k}, \ldots, c_j$ betrokken bij de evaluatie in $x \in [t_j, t_{j+1})$.

---

#### Differentiatie van een splinefunctie

### Stelling 5.3.4

De $r$-de afgeleide van $s(x) = \sum_{i=-k}^{n-1} c_i N_{i,k+1}(x)$ is zelf een spline van graad $k-r$:

$$
s^{(r)}(x) = \sum_{i=-(k-r)}^{n-1} c_i^{(r)}\, N_{i,k+1-r}(x),
$$

waarbij de coëfficiënten recursief worden bepaald door:

$$
c_i^{(0)} = c_i, \qquad c_i^{(r)} = (k+1-r)\, \frac{c_i^{(r-1)} - c_{i-1}^{(r-1)}}{t_{i+k+1-r} - t_i}.
$$

**Evaluatie van de afgeleide in $x$:** Bereken eerst de coëfficiënten $c_i^{(r)}$ via de formule hierboven (een trapeziumvormige tabel), en pas daarna het de Boor-algoritme toe op de resulterende coëfficiënten.

---

## 5.4 Uitbreidingen

### 5.4.1 Samenvallende knooppunten

**Motivatie:** Om functies met knikken of sprongen te modelleren, wil men soms dat een spline minder glad is in bepaalde punten. Dit bereikt men door knooppunten samen te laten vallen.

**Uitbreiding van gedeelde differenties:** Als $x_i = x_{i+1} = \cdots = x_{i+l}$ (meervoudigheid $l+1$):

$$
f[x_i, x_{i+1}, \ldots, x_{i+l}] = \frac{f^{(l)}(x_i)}{l!}.
$$

**Continuïteitsregels voor splines met meervoudige knooppunten:**

| Meervoudigheid $l$ van knooppunt $t$ | Continuïteit van $s$ in $t$ |
|---|---|
| $l = 1$ (enkelvoudig) | $C^{k-1}$: $s, s', \ldots, s^{(k-1)}$ continu |
| $l = 2$ | $C^{k-2}$ |
| $\vdots$ | $\vdots$ |
| $l = k$ | $C^0$: alleen $s$ zelf continu |
| $l = k+1$ | mogelijk discontinu in $s$ zelf |
| $l \geq k+2$ | B-spline wordt nul (uitsluiten) |

**Vuistregel:** Elke extra meervoudigheid kost één graad continuïteit.

**B-splines met samenvallende knooppunten** worden gedefinieerd met de veralgemeende definitie van gedeelde differenties. Het de Boor-algoritme blijft geldig.

**Toepassing:** Met een $(k+1)$-voudig knooppunt aan de rand ($t_0 = t_1 = \cdots = t_k = a$ en $t_n = \cdots = t_{n+k} = b$) bekomt men een B-splinebasis waarbij de eerste en laatste B-spline interpoleren in $a$ respectievelijk $b$. Dit is de standaardkeuze voor niet-periodieke splines.

---

## Samenvatting en sleutelresultaten

### Waarom B-splines?

| Eigenschap | Belang |
|---|---|
| **Lokale drager** (E5.2) | Elke $N_{i,k+1} = 0$ buiten $[t_i, t_{i+k+1}]$; slechts $k+1$ B-splines actief per punt |
| **Positiviteit** (E5.4) | $N_{i,k+1}(x) > 0$ op $(t_i, t_{i+k+1})$; stabiele numerieke eigenschappen |
| **Partitie van eenheid** (E5.6) | $\sum N_i = 1$; convexe romp van controlepunten bevat de spline |
| **Recursiebetrekking** (E5.3) | Efficiënte berekening via de Boor |
| **Stabiele basis** | Geen slechte conditionering, in tegenstelling tot afgeknotte machtsbasis |

### Vergelijking van de twee basissen

| | Afgeknotte-machtsbasis | B-splinebasis |
|---|---|---|
| Dimensie | $n+k$ ✓ | $n+k$ ✓ |
| Numerieke stabiliteit | slecht | goed |
| Lokale drager | neen | ja |
| Evalutie-algoritme | direct | de Boor |

### De Boor vs. de Boor-differentiatie

**Evaluatie:** $k$ ronden van convexe combinaties $\to$ $O(k^2)$ bewerkingen.

**Differentiatie:** Eerst $r$ ronden van verschil-operaties op de coëfficiënten, dan de Boor op de resulterende coëfficiënten.

### Sleutelformules

| Concept | Formule |
|---|---|
| Dimensie splineruimte | $n + k$ |
| Aantal B-splines | $n + k$ (met uitgebreide knooppuntenrij) |
| Uitgebreide knooppuntenrij | $n + 2k + 1$ punten |
| B-splinedefinitie | $M_{i,k+1}(x) = [t_i,\ldots,t_{i+k+1}]_t (t-x)_+^k$ |
| Recursiebetrekking | $N_{i,k+1} = \frac{x-t_i}{t_{i+k}-t_i} N_{i,k} + \frac{t_{i+k+1}-x}{t_{i+k+1}-t_{i+1}} N_{i+1,k}$ |
| Partitie van eenheid | $\sum_{i=-k}^{n-1} N_{i,k+1}(x) = 1$ op $[t_0,t_n]$ |
| Spline in B-splinevorm | $s(x) = \sum_{i=-k}^{n-1} c_i N_{i,k+1}(x)$ |
| Afgeleide coëfficiënten | $c_i^{(1)} = k \cdot \frac{c_i - c_{i-1}}{t_{i+k} - t_i}$ |
