# Modelvragen Deel Samaey — Uitgewerkte Antwoorden

> KU Leuven, Numerieke Benadering 2025–2026.
> Aanpak: eerst begrijpen wat en waarom, dan de formules.

---

## V1 — Kubische splinefunctie

### (a) Definitie

**Wat is een spline?**

Stel je wil een vloeiende curve tekenen door een reeks punten. Als je één grote veelterm over het hele interval gebruikt, gaat die sterk oscilleren tussen de punten — dat heet het Runge-fenomeen. De oplossing is simpel: verdeel het interval in kleine stukjes en gebruik op elk stukje een lage-graads veelterm. Om ervoor te zorgen dat alles mooi aan elkaar kleeft, eisen we dat de stukjes vloeiend op elkaar aansluiten.

Dat is precies wat een spline doet.

**Formeel (Definitie 5.2.1):**

Zij $a = t_0 < t_1 < \cdots < t_n = b$ een rij van knooppunten die het interval opdelen.

Een **splinefunctie van graad $k$** is een functie $s(x)$ op $[a,b]$ die voldoet aan:
1. Op elk stukje $[t_{i-1}, t_i]$ is $s(x)$ een **veelterm van graad $\leq k$**.
2. $s$ en haar afgeleiden $s', s'', \ldots, s^{(k-1)}$ zijn **continu** op het hele interval $[a,b]$.

**Waarom kubisch ($k=3$)?**

Het menselijk oog merkt een sprong in de kromming op. Kromming hangt af van $s''$. Door $s''$ continu te eisen, moeten we tot de tweede afgeleide aansluiten — dat vraagt om graad 3 (kubisch). Zo is een kubische spline visueel altijd vloeiend.

> ⚠️ **Examentip:** Ken de twee eisen: (1) veelterm per stuk, (2) continuïteit t/m de $(k-1)$-de afgeleide.

---

### (b) Dimensie van de splineruimte

**Wat verwacht je intuïtief?**

Je hebt $n$ stukjes, elk met een kubische veelterm (4 coëfficiënten) → normaal zou je $4n$ vrijheidsgraden hebben. Maar de continuïteitseisen koppelen de stukjes aan elkaar en verminderen dat getal. Het resultaat is $n + k$.

**Bewijs (telbewijs):**

**Stap 1 — Totaal zonder continuïteit:**
$n$ stukken, elk met $k+1$ coëfficiënten: $\rightarrow n(k+1)$ vrijheidsgraden.

**Stap 2 — Aftrekken van continuïteitseisen:**
In elk van de $n-1$ inwendige knooppunten $t_1, \ldots, t_{n-1}$ moeten $k$ afgeleiden overeenkomen:
$$\text{beperkingen} = (n-1) \cdot k.$$

**Stap 3 — Resultaat:**
$$\dim = n(k+1) - (n-1)k = nk + n - nk + k = n + k. \quad \square$$

Voor kubische splines ($k=3$): dimensie $= n + 3$.

> ⚠️ **Examentip:** Dit telbewijs staat bijna zeker op het examen. Leer de drie stappen.

---

### (c) Twee basissen: afgeknotte machten en B-splines

#### Afgeknotte-machtsbasis

**Wat is het idee?**

Je wil een functie die op een bepaald knooppunt $t_i$ "inschakelt" en daarvoor nul is. Dat lukt met een afgeknotte machtsfunctie: rechts van $t_i$ gedraagt ze zich als een gewone macht, links is ze nul.

**Formule:**
$$
(x - t_i)_+^k = \begin{cases} 0 & \text{als } x \leq t_i \\ (x - t_i)^k & \text{als } x > t_i \end{cases}
$$

**De volledige basis:**
$$\{1,\, x,\, \ldots,\, x^k,\quad (x-t_1)_+^k,\, \ldots,\, (x-t_{n-1})_+^k\}.$$
Dat zijn $k+1 + (n-1) = n+k$ functies — klopt met de dimensie.

**Nadeel:** Deze basis is numeriek instabiel (slecht geconditioneerd). In de praktijk gebruik je de B-splinebasis.

#### B-splinebasis

**Wat is het idee?**

Een B-spline is een "bobbeltje" dat alleen niet-nul is over een klein stukje van het interval — de zogenaamde **lokale drager**. Dat maakt berekeningen heel efficient: de meeste basisfuncties zijn nul op een gegeven punt.

**Definitie (via gedeelde differentie):**
$$N_{i,k+1}(x) = (t_{i+k+1} - t_i)\, [t_i, \ldots, t_{i+k+1}]_t\, (t - x)_+^k.$$

**Hoe bereken je ze? — Recursiebetrekking:**

Begin met de eenvoudigste B-splines (graad 0: gewoon 1 op een stukje, 0 elders) en bouw daar hogere graden van op:
$$N_{i,1}(x) = \begin{cases} 1 & \text{als } x \in [t_i, t_{i+1}) \\ 0 & \text{anders} \end{cases}$$

$$N_{i,k+1}(x) = \frac{x - t_i}{t_{i+k} - t_i}\, N_{i,k}(x) + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_{i+1}}\, N_{i+1,k}(x).$$

Je combineert telkens twee lagere B-splines met gewichten die optellen tot 1.

**Drie sleuteleigenschappen:**

| Eigenschap | Wat het betekent |
|---|---|
| Lokale drager | $N_{i,k+1}(x) = 0$ buiten $[t_i, t_{i+k+1}]$ — elke B-spline is slechts op een klein stuk actief |
| Positiviteit | $N_{i,k+1}(x) > 0$ binnenin die drager |
| Partitie van eenheid | $\sum_i N_{i,k+1}(x) = 1$ — de B-splines verdelen de "verantwoordelijkheid" keurig |

---

### (d) B-splinebasis voor benadering van data

**Wat is de situatie?**

Je hebt meetpunten $(x_r, f_r)$ en wil een vloeiende splinecurve die zo goed mogelijk door die punten past (kleinste kwadraten).

**Hoe pak je dat aan?**

Schrijf de spline als een lineaire combinatie van B-splines:
$$s(x) = \sum_i c_i\, N_{i,k+1}(x).$$

Dan zoek je de coëfficiënten $c_i$ die de gewogen fout minimaliseren:
$$\min_{c} \sum_{r=1}^N w_r\, \bigl(f_r - s(x_r)\bigr)^2.$$

Differentieer naar $c_j$ en stel gelijk aan nul — dat geeft het **normaalstelsel** $Gc = b$:
$$G_{j,i} = \sum_{r} w_r N_{j}(x_r) N_{i}(x_r), \qquad b_j = \sum_{r} w_r N_{j}(x_r) f_r.$$

**Keuzes die je moet maken:**
- Graad $k$: bijna altijd $k=3$ (kubisch)
- Aantal knooppunten: te weinig → slechte fit, te veel → overfitting
- Ligging knooppunten: idealiter meer knooppunten waar de data sterk varieert

> ⚠️ **Examentip:** **Schönberg-Whitney-voorwaarde** — er moet minstens één datapunt liggen in elk knooppunteninterval. Zonder dit wordt het stelsel singulier.

---

### (e) Spaarsheidspatroon van het lineaire stelsel

**Waarom is de matrix niet vol?**

Elke B-spline is alleen niet-nul over $k+1$ deelintervallen. Dat betekent dat twee B-splines $N_i$ en $N_j$ alleen tegelijk niet-nul zijn als hun dragers overlappen — en dat gebeurt alleen als $|i-j| \leq k$.

**Gevolg:** De Grammatrix $G$ is een **bandmatrix** met bandbreedte $2k+1$. Buiten die band zijn alle elementen nul:
$$G_{j,i} = 0 \quad \text{als } |j - i| > k.$$

**Spaarsheidspatroon voor $k=3$ (bandbreedte 7):**
```
   [* * * *        ]
   [* * * * *      ]
   [* * * * * *    ]
   [* * * * * * *  ]
   [  * * * * * * *]
   [    * * * * * *]
   [      * * * * *]
```

**Praktisch voordeel:** In plaats van $O(N^3)$ bewerkingen voor een vol stelsel, kost Cholesky op een bandmatrix slechts $O(Nk^2)$ — veel sneller!

> ⚠️ **Examentip:** Op het examen: teken gewoon het patroon, je hoeft de waarden niet te kennen.

---

## V2 — Tweedimensionale veeltermbenaderering

### (a) Orthogonale basis voor het 2D-benaderingsprobleem

**Wat is het probleem?**

Je hebt meetdata $f_{i,j}$ op een rechthoekig rooster en wil er een 2D-veelterm $z(x,y) = \sum_{k,l} c_{k,l} x^k y^l$ doorheen passen. Als je gewoon monomiale basisfuncties gebruikt, krijg je een enorm slecht geconditioneerd stelsel (bijna-Hilbert). Oplossing: gebruik orthogonale basisfuncties.

**Het slimme inzicht: tensorproduct**

Als het gewicht scheidbaar is — $w_{i,j} = w_i^{(1)} w_j^{(2)}$ — dan kan je 2D-orthogonaliteit terugbrengen tot twee keer 1D-orthogonaliteit. Je bouwt de 2D-basisfunctie gewoon als een product:
$$\phi_{k,l}(x,y) = p_k(x) \cdot q_l(y),$$
waarbij $\{p_k\}$ orthogonaal zijn in $x$ en $\{q_l\}$ orthogonaal zijn in $y$ (elk apart, via Forsythe).

**Bewijs van 2D-orthogonaliteit:**
$$\sum_{i,j} w_{i,j}\, \phi_{k,l}(x_i, y_j)\, \phi_{u,v}(x_i, y_j) = \underbrace{\left(\sum_i w_i^{(1)} p_k(x_i) p_u(x_i)\right)}_{=0 \text{ als } k \neq u} \cdot \left(\sum_j w_j^{(2)} q_l(y_j) q_v(y_j)\right).$$

Zodra $k \neq u$ (of $l \neq v$) is één van de twee factoren nul → het product is nul → de 2D-basisfuncties zijn orthogonaal. $\square$

> ⚠️ **Examentip:** De sleutel is de scheidbaarheid van het gewicht. Zonder die eigenschap werkt dit niet.

---

### (b) Coëfficiënten van de beste benadering

**Waarom worden de coëfficiënten eenvoudig?**

Normaal moet je een groot normaalstelsel oplossen voor alle coëfficiënten tegelijk. Maar met een orthogonale basis diagonaliseert het normaalstelsel volledig: elke coëfficiënt $a_{k,l}$ kan onafhankelijk van de anderen berekend worden.

**Formule:**
$$a_{k,l} = \frac{\displaystyle\sum_{i,j} w_i^{(1)} w_j^{(2)}\, f_{i,j}\, p_k(x_i)\, q_l(y_j)}{\|p_k\|^2 \|q_l\|^2}.$$

**Herkomst:** Het normaalstelsel is $\sum_{k',l'} a_{k',l'} \langle \phi_{k',l'}, \phi_{k,l} \rangle = \langle f, \phi_{k,l} \rangle$. Door orthogonaliteit zijn alle termen met $(k',l') \neq (k,l)$ nul. Wat overblijft is $a_{k,l} \cdot \|p_k\|^2 \|q_l\|^2 = \langle f, \phi_{k,l} \rangle$, en dus de formule hierboven. $\square$

> ⚠️ **Examentip:** Kernboodschap: orthogonale basis → diagonaal stelsel → elke coëfficiënt apart te berekenen.

---

### (c) Efficiënte berekening van de coëfficiënten

**Het slimme algoritme: herhaalde 1D-benadering**

In plaats van een duur 2D-stelsel op te lossen, kan je de berekening opsplitsen in twee opeenvolgende 1D-stappen:

1. **Stap 1 — langs de $x$-richting:** voor elke kolom $j$ van de data, doe een 1D Forsythe-benadering in $x$:
$$b_{k,j} = \frac{\displaystyle\sum_i w_i^{(1)} f_{i,j}\, p_k(x_i)}{\|p_k\|^2}.$$

2. **Stap 2 — langs de $y$-richting:** voor elke rij $k$ van de tussenresultaten, doe een 1D Forsythe-benadering in $y$:
$$a_{k,l} = \frac{\displaystyle\sum_j w_j^{(2)} b_{k,j}\, q_l(y_j)}{\|q_l\|^2}.$$

```
Data f[i,j] (M×N)  →  b[k,j] ((m+1)×N)  →  a[k,l] ((m+1)×(n+1))
                       Forsythe in x           Forsythe in y
```

---

### (d) Vergelijking rekencomplexiteit

**Directe methode:** lost een normaalstelsel op met $(m+1)(n+1)$ onbekenden. Kost $O((mn)^3)$ — extreem duur bij grote graad.

**Herhaalde 1D Forsythe:** lineaire complexiteit $O(MN(m+n))$ — veel beter!

| Methode | Rekenkost | Stabiliteit |
|---|---|---|
| Directe methode | $O((mn)^3)$ | Slecht (Hilbert-achtige matrix) |
| Herhaalde 1D Forsythe | $O(MN(m+n))$ | Uitstekend (diagonaal stelsel) |

> ⚠️ **Examentip:** Sleutelwoord: "herhaalde eendimensionale benadering". Dankzij scheidbaarheid reduceert 2D naar 1D+1D.

---

## V3 — Discreet benaderingsprobleem: veeltermbenadering

### (a) Scalair product voor discrete kleinste kwadraten

**Wat is het probleem?**

Je hebt $N$ meetpunten $(x_i, f_i)$ en wil een veelterm van lage graad $n \ll N$ die zo goed mogelijk past — niet door de punten heen (dat geeft interpolatie), maar in de kleinste-kwadratenzin.

**Het discrete scalair product** is het wiskundig gereedschap dat "dicht bij de data" formaliseert:
$$(f, g) = \sum_{i=1}^N w_i\, f(x_i)\, g(x_i), \qquad w_i > 0.$$

De gewichten $w_i$ bepalen hoe zwaar elk datapunt meetelt. Het benaderingsprobleem wordt dan: minimaliseer de gewogen fout:
$$E(c) = \|f - y_n\|^2 = \sum_{i=1}^N w_i (f_i - y_n(x_i))^2.$$

---

### (b) Normaalstelsel = overgedetermineerd stelsel

**Twee kanten van dezelfde medaille**

Er zijn twee manieren om hetzelfde probleem te bekijken, en ze geven hetzelfde stelsel.

**Aanpak 1 — Normaalstelsel:**
Schrijf $y_n(x) = \sum_{k=0}^n c_k \phi_k(x)$. Eis dat de fout loodrecht staat op alle basisfuncties:
$$Gc = b, \qquad G_{jk} = \sum_i w_i \phi_j(x_i) \phi_k(x_i), \quad b_j = \sum_i w_i \phi_j(x_i) f_i.$$

**Aanpak 2 — Overgedetermineerd stelsel:**
Stapel alle vergelijkingen $f_i \approx \sum_k c_k \phi_k(x_i)$ op in matrixvorm: $Ac \approx f$ met $A_{ik} = \phi_k(x_i)$. Weeg ze met $D = \text{diag}(\sqrt{w_i})$: $DAc \approx Df$. De normaalvergelijkingen van dit stelsel zijn $A^T W A c = A^T W f$.

**Bewijs van de equivalentie:**
$(A^TWA)_{jk} = \sum_i w_i A_{ij} A_{ik} = \sum_i w_i \phi_j(x_i) \phi_k(x_i) = G_{jk}$.

Beide formuleringen geven letterlijk hetzelfde stelsel. $\square$

> ⚠️ **Examentip:** Dit verband is de brug tussen de abstracte normaalstelsel-theorie en de concrete matrixformulering.

---

### (c) Detecteren van overfitting

**Wat is overfitting?**

Als je een model te complex maakt (te hoge graad), past het perfect op de trainingsdata — inclusief de meetruis. Het model leert de ruis in plaats van de echte structuur. Gevolg: het werkt goed op bekende punten, maar slecht op nieuwe punten.

**Hoe herken je het?**

Splits de data in twee groepen: **trainingsdata** (om het model te bouwen) en **testdata** (om het te evalueren). Kijk hoe beide fouten evolueren als je de graad verhoogt:

```
Fout
 |  E_test
 |    \   /
 |     \ /   ← optimale graad
 |      *
 |        \_____ E_train
 |_____________________ graad n
```

- **Trainingsfout $E_{\text{train}}$:** daalt altijd als je de graad verhoogt (meer vrijheidsgraden).
- **Testfout $E_{\text{test}}$:** daalt eerst, maar stijgt opnieuw bij te hoge graad.

De **optimale graad** is daar waar $E_{\text{test}}$ minimaal is.

---

### (d) Tikhonov regularisatie

**Wat is het idee?**

In plaats van de complexiteit te beperken door de graad vast te kiezen, voeg je een **strafterm** toe aan de kostfunctie die grote coëfficiënten ontmoedigt. Grote coëfficiënten geven sterk oscillerende veeltermen — precies wat we willen vermijden.

**De geregulariseerde kostfunctie:**
$$E_{\text{reg}}(c) = \underbrace{\sum_{i} w_i (f_i - y_n(x_i))^2}_{\text{trouw aan de data}} + \lambda \underbrace{\|c\|_2^2}_{\text{straf op grote coëfficiënten}}, \qquad \lambda > 0.$$

De parameter $\lambda$ regelt de balans:

| $\lambda$ | Effect |
|---|---|
| $\lambda = 0$ | Geen regularisatie — gewone kleinste kwadraten |
| $\lambda$ klein | Lichte straf op grote coëfficiënten |
| $\lambda \to \infty$ | Alle $c_k \to 0$ (het model vervalt volledig) |

---

### (e) Wijziging aan het normaalstelsel

**Wat verandert er?**

De toevoeging van de strafterm $\lambda\|c\|^2$ geeft slechts één kleine aanpassing aan het normaalstelsel: er wordt $\lambda I$ opgeteld bij de normaalmatrix.

**Afleiding:**
Schrijf $E_{\text{reg}} = \|D(Ac - f)\|^2 + \lambda\|c\|^2$. Uitwerken en differentiëren naar $c$:
$$\nabla_c E_{\text{reg}} = 2(A^TWA + \lambda I)c - 2A^TWf = 0.$$

**Het geregulariseerde normaalstelsel:**
$$\boxed{(A^TWA + \lambda I)\hat{c} = A^TWf.}$$

**Waarom is dit altijd inverteerbaar?**
De eigenwaarden van $A^TWA + \lambda I$ zijn $\sigma_i^2 + \lambda \geq \lambda > 0$. De matrix is dus altijd positief definiet, zelfs als $A$ niet vol rang heeft.

| | Gewoon normaalstelsel | Geregulariseerd |
|---|---|---|
| Matrix | $A^TWA$ | $A^TWA + \lambda I$ |
| Inverteerbaar? | Niet altijd | Altijd (voor $\lambda > 0$) |
| Conditiegetal | Hoog | Lager |

> ⚠️ **Examentip:** De aanpassing is enkel $+\lambda I$ op de diagonaal. Klein verschil, grote gevolgen.

---

## V4 — Overfitting en regularisatie

### (a) Overfitting: definitie en detectie

**Simpel gezegd:**

Overfitting treedt op als een model te nauw aansluit bij de ruizige trainingsdata in plaats van bij de onderliggende functie. Het leert de ruis in plaats van de structuur.

Symptomen:
- Lage fout op de trainingsdata
- Hoge fout op nieuwe (test)data

**Hoe vaststelllen:**
- Splits data in trainings- en testset
- Evalueer $E_{\text{test}}$ als functie van de modelcomplexiteit
- Overfitting = $E_{\text{test}}$ stijgt terwijl $E_{\text{train}}$ blijft dalen

**Alternatief:** Leave-one-out kruisvalidatie — bouw $N$ modellen, elk met één punt weggelaten, en evalueer op dat weggelaten punt.

> ⚠️ **Examentip:** Onderscheid van onderfitting: bij onderfitting zijn *beide* fouten hoog.

---

### (b) Tikhonov vs. LASSO

**Beide regulariseren, maar op een andere manier**

Beide methodes voegen een strafterm $\lambda F(c)$ toe, maar de keuze van $F$ bepaalt het gedrag volledig.

**Tikhonov (L2-regularisatie):** straft de *grootte* van elke coëfficiënt.
- Alle coëfficiënten worden kleiner, maar nooit precies nul.
- Leidt tot een gewoon lineair stelsel: $(A^TWA + \lambda I)\hat{c} = A^TWf$.

**LASSO (L1-regularisatie):** straft de *som van de absolute waarden*.
- Duwt sommige coëfficiënten naar precies nul — de irrelevante features vallen weg.
- Leidt tot een niet-lineair stelsel (moet iteratief opgelost worden).

| Eigenschap | Tikhonov (L2) | LASSO (L1) |
|---|---|---|
| Coëfficiënten | Klein maar nooit nul | Vaak exact nul |
| Stelsel | Lineair | Niet-lineair |
| Variabelenselectie | Nee | Ja (automatisch) |

**Waarom is LASSO schaars? — Geometrisch argument:**

Stel je het optimalisatieprobleem voor als: minimaliseer $E(c)$ met als beperking dat $F(c) \leq K$.

```
   L2-bol (Tikhonov)          L1-bol (LASSO)

         ___                       /\
        /   \                     /  \
       | c*  |   raakpunt         / c* \   raakpunt in hoek
        \___/    willekeurig     \    /   (c₂ = 0!)
                                  \  /
                                   \/
```

De L1-bol heeft scherpe hoeken op de assen. De elliptische niveaukrommen van $E(c)$ raken de L1-bol bijna altijd in zo'n hoek, waardoor een coëfficiënt precies nul wordt. Bij de gladde L2-bol is dat niet het geval.

> ⚠️ **Examentip:** Ken dit geometrisch argument — het staat bijna zeker op het examen.

---

### (c) Regularisatie bij diepe neurale netwerken

**Waarom regulariseren in neurale netwerken?**

Grote gewichten $W^{[l]}$ maken neuronen gevoelig voor kleine veranderingen in de invoer — een teken van overfitting. Door grote gewichten te bestraffen, wordt het netwerk robuuster.

**De geregulariseerde kostfunctie:**
$$\text{Kost}_{\text{reg}} = \frac{1}{N}\sum_{i} \|\mathbf{y}^{(i)} - \mathbf{a}^{[L](i)}\|^2 + \frac{\lambda}{N}\sum_{l=2}^L \|W^{[l]}\|_2^2.$$

**Merk op:** de biases $\mathbf{b}^{[l]}$ worden *niet* geregulariseerd — die bepalen enkel een verschuiving, niet de gevoeligheid.

**Andere technieken:**
- **Dropout:** willekeurig neuronen tijdelijk uitschakelen tijdens het trainen
- **Early stopping:** stop zodra $E_{\text{test}}$ begint te stijgen
- **Data augmentatie:** kunstmatig meer trainingsdata aanmaken

---

## V5 — Diep neuraal netwerk

### (a) Architectuur: 3 lagen, 2 inputs, 3 verborgen, 2 outputs

**Visueel:**
```
Invoer [1]     Verborgen [2]     Uitvoer [3]

  x₁ ○——————→ ○ n₂₁ ——→ ○ n₃₁ → y₁
     |        ↗ ○ n₂₂        ↗
     |       ↗  ○ n₂₃ ——→ ○ n₃₂ → y₂
  x₂ ○——————
```
(Elke laag is volledig verbonden met de volgende.)

**Aantal parameters:**
- Laag 2: $W^{[2]} \in \mathbb{R}^{3 \times 2}$ (6 gewichten) + $\mathbf{b}^{[2]} \in \mathbb{R}^3$ (3 biases) = **9 parameters**
- Laag 3: $W^{[3]} \in \mathbb{R}^{2 \times 3}$ (6 gewichten) + $\mathbf{b}^{[3]} \in \mathbb{R}^2$ (2 biases) = **8 parameters**
- **Totaal: 17 parameters**

---

### (b) Voorwaartse propagatie

**Wat gebeurt er?**

Een invoer $\mathbf{x}$ wordt laag per laag doorheen het netwerk "gepropageerd". In elke laag doe je een lineaire transformatie (gewichten + bias), gevolgd door een niet-lineaire activatiefunctie $\sigma$.

**Berekening:**

Laag 2:
$$\mathbf{z}^{[2]} = W^{[2]} \mathbf{x} + \mathbf{b}^{[2]}, \qquad \mathbf{a}^{[2]} = \sigma(\mathbf{z}^{[2]}).$$

Laag 3:
$$\mathbf{z}^{[3]} = W^{[3]} \mathbf{a}^{[2]} + \mathbf{b}^{[3]}, \qquad \mathbf{a}^{[3]} = \sigma(\mathbf{z}^{[3]}) = F(\mathbf{x}).$$

**Dimensies expliciet — schrijf deze altijd op bij het examen:**
$$\underbrace{W^{[2]}}_{3 \times 2} \cdot \underbrace{\mathbf{x}}_{2 \times 1} + \underbrace{\mathbf{b}^{[2]}}_{3 \times 1} = \underbrace{\mathbf{z}^{[2]}}_{3 \times 1} \xrightarrow{\sigma} \underbrace{\mathbf{a}^{[2]}}_{3 \times 1}$$
$$\underbrace{W^{[3]}}_{2 \times 3} \cdot \underbrace{\mathbf{a}^{[2]}}_{3 \times 1} + \underbrace{\mathbf{b}^{[3]}}_{2 \times 1} = \underbrace{\mathbf{z}^{[3]}}_{2 \times 1} \xrightarrow{\sigma} \underbrace{\mathbf{a}^{[3]}}_{2 \times 1}$$

**Samengestelde formule:**
$$F(\mathbf{x}) = \sigma\!\Bigl(W^{[3]} \sigma\!\bigl(W^{[2]} \mathbf{x} + \mathbf{b}^{[2]}\bigr) + \mathbf{b}^{[3]}\Bigr).$$

> ⚠️ **Examentip:** Schrijf bij elke stap de dimensies expliciet.

---

### (c) Trainingsalgoritme bij veel datapunten

**Het probleem:**

De volledige gradiënt over $N$ punten berekenen kost $O(N)$ per stap — veel te duur als $N$ groot is.

**De oplossing: mini-batch SGD**

Gebruik op elke stap slechts een kleine, willekeurige steekproef (mini-batch) van $m \ll N$ punten om de gradiënt te schatten. Dat is veel goedkoper en werkt verrassend goed.

**Algoritme:**
```
voor elke epoch:
    permuteer de N trainingspunten willekeurig
    voor k = 1 tot N/m:
        kies mini-batch van m punten
        schat gradiënt via backpropagation: g ≈ (1/m) Σ ∇C(c)
        update: c ← c - α · g
```

**Te kiezen hyperparameters:**
- **Learning rate $\alpha$:** te groot → divergentie; te klein → trage convergentie
- **Batchgrootte $m$:** grote batch = nauwkeurigere gradiënt maar duurder per stap
- **Aantal epochs:** hoe vaak je de volledige dataset doorloopt
- **Regularisatieparameter $\lambda$**

**Voordeel van mini-batch SGD:** de ruis in de gradiëntschatting helpt paradoxaal genoeg om aan lokale minima te ontsnappen.

> ⚠️ **Examentip:** Ken de terminologie: epoch, mini-batch, learning rate, backpropagation.

---

### (d) Backpropagation: wat en waarom?

**Wat is het?**

Backpropagation is de slimme manier om de gradiënt van de kostfunctie t.o.v. alle gewichten te berekenen. Het past de kettingregel toe, maar dan van achter naar voor doorheen de lagen — zodat tussenresultaten hergebruikt worden.

**Hoe?**

Sla tijdens de voorwaartse pass $\mathbf{z}^{[l]}$ en $\mathbf{a}^{[l]}$ op. Bereken dan achterwaarts de "fout per laag" $\boldsymbol{\delta}^{[l]}$:
$$\boldsymbol{\delta}^{[L]} = \nabla_{\mathbf{a}^{[L]}} C \odot \sigma'(\mathbf{z}^{[L]}),$$
$$\boldsymbol{\delta}^{[l]} = \bigl(W^{[l+1]}\bigr)^T \boldsymbol{\delta}^{[l+1]} \odot \sigma'(\mathbf{z}^{[l]}).$$

De gradiënten volgen direct:
$$\frac{\partial C}{\partial W^{[l]}} = \boldsymbol{\delta}^{[l]} \bigl(\mathbf{a}^{[l-1]}\bigr)^T, \qquad \frac{\partial C}{\partial \mathbf{b}^{[l]}} = \boldsymbol{\delta}^{[l]}.$$

**Waarom zo snel?**

Naïef zou je één voorwaartse pass per parameter nodig hebben: $O(P)$ passes voor $P$ parameters. Backpropagation hergebruikt de $\boldsymbol{\delta}^{[l]}$-waarden en heeft dezelfde kost als **slechts één voorwaartse pass** — ongeacht het aantal parameters.

---

### (e) Waarom stochastische steilste helling?

**Kort antwoord:** de exacte gradiënt over $N$ punten berekenen is $O(N)$ per stap. SGD vervangt dit door een mini-batch van $m \ll N$ punten: $O(m)$ per stap. Daardoor kan je veel meer stappen zetten in dezelfde tijd. Bovendien helpt de ruis in de gradiëntschatting om lokale minima te ontwijken.

---

## V6 — Methode van de Toegevoegde Gradiënten (Conjugate Gradients)

**Probleemstelling:** Los $Ax = b$ op met $A$ **symmetrisch positief definiet** (SPD).

### (a) Algoritme in pseudocode + uitleg

**Wat is CG?**

Stel je probeert een punt te vinden dat een kwadraat functie $f(x) = \frac{1}{2}x^TAx - x^Tb$ minimaliseert. Dat is wiskundig equivalent met $Ax = b$ oplossen (want $\nabla f = Ax - b = 0$). CG doet dit door slim gekozen zoekrichtingen te gebruiken die "toegevoegd" (conjugaat) zijn — ze storen elkaars voortgang niet.

**Pseudocode:**
```
Start: x⁽⁰⁾ = 0,   r⁽⁰⁾ = p⁽⁰⁾ = b
voor k = 1, 2, ...:
    α⁽ᵏ⁾ = ‖r⁽ᵏ⁻¹⁾‖² / (p⁽ᵏ⁻¹⁾ᵀ A p⁽ᵏ⁻¹⁾)   ← optimale stapgrootte
    x⁽ᵏ⁾ = x⁽ᵏ⁻¹⁾ + α⁽ᵏ⁾ · p⁽ᵏ⁻¹⁾             ← stap in zoekrichting
    r⁽ᵏ⁾ = r⁽ᵏ⁻¹⁾ - α⁽ᵏ⁾ · A p⁽ᵏ⁻¹⁾            ← update residu
    als ‖r⁽ᵏ⁾‖ < ε: stop
    β⁽ᵏ⁾ = ‖r⁽ᵏ⁾‖² / ‖r⁽ᵏ⁻¹⁾‖²                ← conjugeringsfactor
    p⁽ᵏ⁾ = r⁽ᵏ⁾ + β⁽ᵏ⁾ · p⁽ᵏ⁻¹⁾               ← nieuwe zoekrichting
```

**Uitleg van de rollen:**
- $r^{(k)} = b - Ax^{(k)}$: het **residu** — hoe ver we nog van de oplossing zijn
- $p^{(k)}$: de **zoekrichting** — een mix van het residu en de vorige richting om "zigzag" te vermijden
- $\alpha^{(k)}$: de **optimale stapgrootte** langs $p^{(k)}$ (exacte lijnzoeking)
- $\beta^{(k)}$: zorgt dat opeenvolgende richtingen A-toegevoegd blijven

---

### (b) CG bouwt een Krylov-ruimte op

**Wat is een Krylov-ruimte?**

Na $k$ stappen heeft CG alleen vectoren gebruikt die te schrijven zijn als veeltermen in $A$ toegepast op $b$: $b, Ab, A^2b, \ldots$ De ruimte opgespannen door deze vectoren heet de Krylov-ruimte:
$$\mathcal{K}_k(A, b) = \text{span}\{b, Ab, A^2b, \ldots, A^{k-1}b\}.$$

**Stelling:** Na $k-1$ stappen geldt:
$$\mathcal{K}_k(A, b) = \text{span}\{p^{(0)}, \ldots, p^{(k-1)}\} = \text{span}\{r^{(0)}, \ldots, r^{(k-1)}\}.$$

**Bewijs (inductie):**

Basisgeval: $p^{(0)} = r^{(0)} = b$ → $\mathcal{K}_1 = \langle b \rangle$. ✓

Inductiestap: stel geldig voor stap $k-1$. Dan:
$$r^{(k-1)} = r^{(k-2)} - \alpha^{(k-1)} A p^{(k-2)} \in \mathcal{K}_{k-1} + A \mathcal{K}_{k-1} = \mathcal{K}_k.$$
$$p^{(k-1)} = r^{(k-1)} + \beta^{(k-1)} p^{(k-2)} \in \mathcal{K}_k. \quad \square$$

> ⚠️ **Examentip:** Sleutel: vermenigvuldiging met $A$ verhoogt de Krylov-graad met 1.

---

### (c) CG vindt de beste benadering in de Krylov-ruimte

**Wat betekent "beste benadering"?**

$x^{(k)}$ is het element van $\mathcal{K}_k(A,b)$ dat de **A-norm van de fout** minimaliseert:
$$\|e^{(k)}\|_A = \min_{x \in \mathcal{K}_k(A,b)} \|x - x^*\|_A, \qquad \|v\|_A = \sqrt{v^T A v}.$$

**Bewijs:**

Schrijf een willekeurig punt als $x = x^{(k)} + \delta x$ met $\delta x \in \mathcal{K}_k$:
$$\|x - x^*\|_A^2 = \|e^{(k)}\|_A^2 + 2 e^{(k)T} A\, \delta x + \|\delta x\|_A^2.$$

Nu is $Ae^{(k)} = Ax^{(k)} - b = -r^{(k)}$. De residu's staan loodrecht op de Krylov-ruimte (Eigenschap 11.3): $r^{(k)T} \delta x = 0$.

Dus de gekruiste term is nul, en $\|\delta x\|_A^2 \geq 0$ → minimum voor $\delta x = 0$, i.e. $x = x^{(k)}$. $\square$

---

### (d) Verband met veeltermbenadering + versnelling

**Het veeltermverband:**

De A-normfout bij stap $k$ is equivalent met een veeltermbenaderingsprobleem:
$$\|e^{(k)}\|_A^2 = \min_{p \in \mathcal{N}_k} \|p(A) e^{(0)}\|_A^2, \qquad \mathcal{N}_k = \{p : \deg p \leq k,\ p(0) = 1\}.$$

**Wat betekent dit praktisch?**

De convergentie hangt ervan af hoe goed je de eigenwaarden van $A$ kunt "onderdrukken" met een veelterm die in $0$ gelijk is aan $1$. Dat gaat beter als:
- $A$ weinig **verschillende eigenwaarden** heeft → CG convergeert exact in $k$ stappen (met $k$ = aantal eigenwaarden)
- De eigenwaarden **geclusterd** zijn → de veelterm is makkelijker te vinden

**Preconditioning: hoe convergentie versnellen?**

Kies een inverteerbare matrix $D$ zodat $D^{-1}AD^{-T}$ beter geclusterde eigenwaarden heeft. Los dan het equivalent geconditioneerde stelsel op:
$$(D^{-1}AD^{-T}) y = D^{-1}b, \qquad y = D^T x.$$

Een goede preconditioner $D$ maakt het geconditioneerde systeem bijna als de eenheidsmatrix → CG convergeert in weinig iteraties.

> ⚠️ **Examentip:** Kernboodschap: cluster de eigenwaarden → snellere convergentie.

---

## V7 — Niet-lineair optimalisatieprobleem

### (a) Sterkste daling = richting van de negatieve gradiënt

**Intuïtie:**

De gradiënt $\nabla f(x)$ staat loodrecht op de niveauverzameling door $x$ en wijst de richting van sterkste *stijging*. De richting van sterkste *daling* is dus $-\nabla f(x)$.

**Formeel bewijs:**

De richtingsafgeleide in richting $s$ (met $\|s\|_2 = 1$) is:
$$\frac{d}{d\alpha} f(x + \alpha s)\bigg|_{\alpha=0} = s^T \nabla f(x).$$

We willen dit zo **negatief mogelijk** maken. Via de Cauchy-Schwarz-ongelijkheid:
$$s^T \nabla f(x) \geq -\|s\|_2 \cdot \|\nabla f(x)\|_2 = -\|\nabla f(x)\|_2,$$
met gelijkheid wanneer $s = -\nabla f(x) / \|\nabla f(x)\|_2$.

**Conclusie:** de sterkste daling treedt op in de richting $s = -\nabla f(x)$. $\square$

> ⚠️ **Examentip:** Dit is Eigenschap 11.1. Ken het bewijs: richtingsafgeleide + Cauchy-Schwarz.

---

### (b) Lijnzoekmethode vs. Armijo-backtracking

**Wat is het doel?**

Bij gradiëntafdaling moet je de stapgrootte $\alpha^{(k)}$ kiezen in de iteratie $x^{(k)} = x^{(k-1)} - \alpha^{(k)} \nabla f(x^{(k-1)})$.

#### Exacte lijnzoekmethode

**Idee:** kies $\alpha^{(k)}$ zo dat $f$ zo laag mogelijk is langs de gekozen richting:
$$\alpha^{(k)} = \arg\min_{\alpha > 0} f(x^{(k-1)} + \alpha p^{(k)}).$$

**Nadeel:** vereist een extra optimalisatieprobleem per stap, en geeft bij steilste afdaling zigzag-convergentie (elke stap staat loodrecht op de vorige).

**Bewijs van loodrechtheid:** na exacte lijnzoeking is $\nabla f(x^{(k)})^T p^{(k)} = 0$. Maar de volgende richting is $p^{(k+1)} = -\nabla f(x^{(k)})$, dus $p^{(k+1)} \perp p^{(k)}$.

#### Backtracking met de Armijo-voorwaarde

**Idee:** geen exacte minimalisatie — slechts een *voldoende* daling eisen. Dit is goedkoper en werkt in de praktijk even goed.

**De Armijo-voorwaarde:**
$$f\!\left(x + \alpha p\right) \leq f(x) + c \cdot \alpha \cdot \nabla f(x)^T p, \qquad c \in (0, 1).$$

Intuïtief: de nieuwe functiewaarde moet **onder de gecorrigeerde raaklijn** liggen (met helling $c$ keer de richtingsafgeleide).

**Algoritme:** start met $\alpha = 1$, halveer tot de voorwaarde voldaan is.

**Drie eigenschappen die Samaey test:**

1. **Stopt altijd:** de raaklijn $f(x) + \alpha \nabla f \cdot p$ heeft negatieve helling (want $\nabla f \cdot p < 0$ bij afdaling). De gecorrigeerde raaklijn is minder steil. Voor kleine genoeg $\alpha$ valt de functiewaarde daar zeker onder → algoritme stopt in eindig aantal stappen.

2. **Geen te grote stappen:** het rechterlid daalt lineair in $\alpha$ (want $c \cdot \nabla f \cdot p < 0$). Voor grote $\alpha$ eist de voorwaarde een grote daling in $f$ → vermijdt stappen die over het minimum springen.

3. **Groot aanvaardingsgebied ver van optimum:** ver van het minimum is $|\nabla f \cdot p|$ groot → de gecorrigeerde raaklijn daalt snel → er is een groot bereik van $\alpha$-waarden die de voorwaarde voldoen → weinig halveerstappen nodig.

**Vergelijking:**

| | Exacte lijnzoeking | Armijo-backtracking |
|---|---|---|
| Rekenkost | Hoog (extra minimaliseringsprobleem) | Laag (enkel functie-evaluaties) |
| Garantie | Optimale stap langs richting | Voldoende daling |
| Zigzag | Ja (bij steilste afdaling) | Minder uitgesproken |
| Gebruik | Theoretische analyse | Praktische implementaties |

> ⚠️ **Examentip:** Ken de Armijo-voorwaarde letterlijk en weet de drie geometrische eigenschappen te benoemen.

---

## Samenvatting Sleutelformules

| Concept | Formule |
|---|---|
| Dimensie splineruimte | $n + k$ |
| B-splinerecursie | $N_{i,k+1} = \frac{x-t_i}{t_{i+k}-t_i}N_{i,k} + \frac{t_{i+k+1}-x}{t_{i+k+1}-t_{i+1}}N_{i+1,k}$ |
| Partitie van eenheid | $\sum_i N_{i,k+1}(x) = 1$ |
| 2D-orthogonale coëfficiënten | $a_{k,l} = \frac{\langle f, p_k q_l \rangle}{\|p_k\|^2 \|q_l\|^2}$ |
| Tikhonov normaalstelsel | $(A^TWA + \lambda I)\hat{c} = A^TWf$ |
| Voorwaartse propagatie | $\mathbf{a}^{[l]} = \sigma(W^{[l]}\mathbf{a}^{[l-1]} + \mathbf{b}^{[l]})$ |
| CG-stapgrootte | $\alpha^{(k)} = \frac{\|r^{(k-1)}\|^2}{p^{(k-1)T}Ap^{(k-1)}}$ |
| CG-conjugering | $\beta^{(k)} = \frac{\|r^{(k)}\|^2}{\|r^{(k-1)}\|^2}$ |
| Richtingsafgeleide | $\nabla_s f(x) = s^T \nabla f(x)$ |
| Armijo-voorwaarde | $f(x + \alpha p) \leq f(x) + c\alpha\nabla f^T p$, $c \in (0,1)$ |
