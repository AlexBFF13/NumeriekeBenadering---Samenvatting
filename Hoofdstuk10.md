# Deel 10 – Niet-lineaire benaderingsproblemen en neurale netwerken

## Overzicht

In de eerste delen van de cursus was elk benaderingsprobleem **lineair**: de te benaderen functie werd voorgesteld als een lineaire combinatie van *vaste* basisfuncties, en de enige onbekenden waren de coëfficiënten. Dat leidde steeds tot een lineair stelsel.

Dit hoofdstuk verbreekt die aanname. In **niet-lineaire** benaderingsproblemen hangen de basisfuncties zelf ook af van onbekende parameters. De onbekenden zijn dan niet langer lineair gekoppeld, en het resulterende optimalisatieprobleem is niet langer een lineair stelsel — het is een **niet-lineair, niet-convex optimalisatieprobleem**.

Het meest concrete en actuele voorbeeld van zo'n aanpak is het **diep neuraal netwerk** (DNN): een samengestelde, niet-lineaire functie waarvan de parameters (gewichten en biases) worden bepaald door een optimalisatieprocedure die men **trainen** noemt.

**Rode draad:**
1. Wat maakt een benaderingsprobleem niet-lineair? (§10.1)
2. De bouwblokken van een neuraal netwerk: activatiefuncties en lagen (§10.2.1–10.2.2)
3. Hoe werkt een volledig netwerk? Uitgewerkt voorbeeld (§10.2.3)
4. De kostfunctie en het optimalisatieprobleem (§10.2.4)
5. Uitbreidingen: convolutionele netwerken en regularisatie (§10.3)

---

## 10.1 Benaderingsproblemen met niet-lineaire parameterafhankelijkheid

### Wanneer is een benaderingsprobleem niet-lineair?

Een lineair benaderingsprobleem heeft de vorm:

$$y_n(x) = \sum_{k=0}^n c_k \phi_k(x),$$

waarbij de basisfuncties $\phi_k(x)$ **volledig gespecificeerd** zijn — enkel de coëfficiënten $c_k$ zijn onbekend.

Een niet-lineair benaderingsprobleem ontstaat wanneer de basisfuncties zelf ook **onbekende parameters** bevatten. De onbekenden zitten dan niet langer lineair in de uitdrukking.

### Voorbeeld 1: periodiek signaal met onbekende frequenties

Gegeven een tijdreeks $(t_i, y_i)_{i=1}^N$. Als het signaal periodiek is, kan men benaderen door:

$$y_n(t) = \sum_{k=0}^n a_k \sin(\omega_k t + \varphi_k). \tag{10.1}$$

- Als $\omega_k$ en $\varphi_k$ **bekend** zijn → lineair probleem in $a_k$.
- Als $\omega_k$ en/of $\varphi_k$ ook **onbekend** zijn → niet-lineair probleem: de functie $\sin(\omega_k t + \varphi_k)$ hangt niet-lineair af van de parameters $\omega_k$ en $\varphi_k$.

Het kleinste-kwadratenprobleem wordt dan:

$$(\hat{a}_k, \hat{\omega}_k, \hat{\varphi}_k)_{k=0}^n = \underset{a_k,\omega_k,\varphi_k}{\arg\min} \sum_{i=1}^N \left|y_i - y_n(t_i)\right|^2. \tag{10.2}$$

Dit is een **niet-lineair optimalisatieprobleem**.

### Voorbeeld 2: radioactief verval met onbekende vervalkonstanten

Gegeven stralingsdata $(t_i, y_i)_{i=1}^N$ van een mengsel van radioactieve stoffen:

$$y_n(t) = \sum_{k=0}^n a_k e^{-\lambda_k t}. \tag{10.3}$$

- Als de stoffen (en dus $\lambda_k$) **bekend** zijn → lineair probleem in $a_k$.
- Als $\lambda_k$ ook **onbekend** zijn → niet-lineair probleem.

De coëfficiënten $a_k$ geven de beginhoeveelheden van elke stof weer; de $\lambda_k$ zijn gerelateerd aan de halfwaardetijden.

### Waarschuwing: niet-lineaire transformaties zijn gevaarlijk

Het kan verleidelijk zijn om het niet-lineaire probleem te lineariseren via een transformatie. Voorbeeld: $y(t) = a e^{-\lambda t}$ wordt via logaritme $\log y = \log a - \lambda t$ — een lineair probleem in $\log a$ en $\lambda$.

**Dit is een slecht idee.** De kleinste-kwadratenaanpak veronderstelt dat de meetfouten een bepaalde structuur hebben (typisch: normale verdeling). Een niet-lineaire transformatie (zoals de logaritme) **vervormt de structuur van de meetfouten**, waardoor de kleinste-kwadratenoplossing op de getransformeerde data iets anders minimaliseert dan wat je eigenlijk wil.

> ⚠️ **Belangrijk voor examen:** Weet het verschil tussen lineaire en niet-lineaire benaderingsproblemen. Weet waarom niet-lineaire transformaties om te lineariseren gevaarlijk zijn.

---

## 10.2 Diepe neurale netwerken

### Motiverende taak: classificatie

Beschouw een verzameling gelabelde datapunten in $\mathbb{R}^2$: sommige behoren tot categorie A (cirkels), andere tot categorie B (kruisjes). Het doel: construeer een functie $F: \mathbb{R}^2 \to \mathbb{R}^2$ die elk nieuw punt correct classificeert.

Dit is **gesuperviseerd leren**: we leren van data waarbij de labels (categorieën) door mensen zijn toegekend.

### 10.2.1 Activatiefuncties

**Intuïtie:** Het idee van een neuraal netwerk is om via herhaalde toepassing van een eenvoudige niet-lineaire functie (de activatiefunctie) zeer complexe, niet-lineaire afbeeldingen te kunnen construeren.

Enkele gangbare activatiefuncties:

| Naam | Voorschrift | Eigenschappen |
|---|---|---|
| Lineair | $f(x) = x$ | Geen niet-lineariteit — nutteloos voor diepte |
| Binaire stap | $f(x) = \mathbf{1}[x \geq 0]$ | Niet differentieerbaar in $x=0$ |
| **Sigmoïde** | $f(x) = \dfrac{1}{1+e^{-x}}$ | Glad, output $\in (0,1)$, goedkope afgeleide |
| Tanh | $f(x) = \tanh(x)$ | Output $\in (-1,1)$, symmetrisch |
| Arctangens | $f(x) = \arctan(x)$ | Vergelijkbaar met tanh |
| **ReLU** | $f(x) = \max(0, x)$ | Stuksgewijs lineair, meest gebruikt in de praktijk |

**De sigmoïde in detail:**

$$f(x) = \frac{1}{1+e^{-x}}. \tag{10.8}$$

Een bijzonder handig kenmerk: de afgeleide kan uitgedrukt worden in termen van de functie zelf:

$$f'(x) = f(x)\bigl(1 - f(x)\bigr). \tag{10.9}$$

Dit maakt de berekening van afgeleiden goedkoop — want $f(x)$ is toch al berekend in de voorwaartse doorgang.

**Schaling en verschuiving:** Via een lineaire transformatie op de invoer kan men de locatie en helling van de activatiefunctie aanpassen:

$$f(ax + b), \quad a \text{ = schaalfactor}, \quad b \text{ = bias (verschuiving)}.$$

> ⚠️ **Belangrijk voor examen:** Ken de sigmoïde en haar afgeleide. Begrijp waarom niet-lineariteit van de activatiefunctie essentieel is (zonder het zijn alle lagen samen één lineaire afbeelding).

### 10.2.2 Lagen van neuronen

Een **laag** van neuronen neemt een invoervector $\mathbf{x}$ en produceert een uitvoervector via:

$$\text{output} = \sigma(W\mathbf{x} + \mathbf{b}), \tag{10.10}$$

waarbij:
- $W \in \mathbb{R}^{n_{\text{out}} \times n_{\text{in}}}$: de **gewichtsmatrix** — bepaalt hoe sterk elk invoerneuron bijdraagt aan elk uitvoerneuron
- $\mathbf{b} \in \mathbb{R}^{n_{\text{out}}}$: de **biasvector** — voegt een verschuiving toe per uitvoerneuron
- $\sigma$: de activatiefunctie, componentsgewijs toegepast

**Detailniveau van één neuron $i$:**

$$\text{output}_i = \sigma\!\left(\sum_j W_{ij} x_j + b_i\right).$$

Elk neuron vormt een gewogen som van alle inkomende signalen, voegt een bias toe, en past de activatiefunctie toe. Dit is de directe analogie met een biologisch neuron: de gewogen som is de "activatiedrempel", de activatiefunctie bepaalt of het neuron "vuurt".

### 10.2.3 Uitgewerkt voorbeeld: 4-laags netwerk

Beschouw een netwerk met vier lagen voor het classificatieprobleem in $\mathbb{R}^2$:

| Laag | Neuronen | Invoerdim. | Uitvoerdim. | Parameters |
|---|---|---|---|---|
| 1 (input) | 2 | — | $\mathbb{R}^2$ | — |
| 2 | 2 | $\mathbb{R}^2$ | $\mathbb{R}^2$ | $W^{[2]} \in \mathbb{R}^{2\times2}$, $\mathbf{b}^{[2]} \in \mathbb{R}^2$ → 6 params |
| 3 | 3 | $\mathbb{R}^2$ | $\mathbb{R}^3$ | $W^{[3]} \in \mathbb{R}^{3\times2}$, $\mathbf{b}^{[3]} \in \mathbb{R}^3$ → 9 params |
| 4 (output) | 2 | $\mathbb{R}^3$ | $\mathbb{R}^2$ | $W^{[4]} \in \mathbb{R}^{2\times3}$, $\mathbf{b}^{[4]} \in \mathbb{R}^2$ → 8 params |

**Totaal: 23 parameters.**

De volledige afbeelding van het netwerk is:

$$F(\mathbf{x}) = \sigma\!\Bigl(W^{[4]}\,\sigma\!\Bigl(W^{[3]}\,\sigma(W^{[2]}\mathbf{x} + \mathbf{b}^{[2]}) + \mathbf{b}^{[3]}\Bigr) + \mathbf{b}^{[4]}\Bigr) \in \mathbb{R}^2. \tag{10.11}$$

**Interpretatie van de uitvoer:** De twee componenten van $F(\mathbf{x})$ stellen de "scores" voor categorie A en B voor. Classificeer $\mathbf{x}$ als categorie A als $F_1(\mathbf{x}) > F_2(\mathbf{x})$, en als B anders.

**Kernpunt:** De samengesteldheid van lagen (diepte) stelt het netwerk in staat zeer complexe, niet-lineaire beslissingsgrenzen te leren — iets wat geen enkele lineaire methode kan.

### 10.2.4 De kostfunctie en het optimalisatieprobleem

**Doelwaarden:** Voor elk trainingsdata-punt $\mathbf{x}^{(i)}$ definiëren we de gewenste uitvoer:

$$\mathbf{y}(\mathbf{x}^{(i)}) = \begin{cases} [1, 0]^T & \text{als } \mathbf{x}^{(i)} \in \text{categorie A} \\ [0, 1]^T & \text{als } \mathbf{x}^{(i)} \in \text{categorie B.} \end{cases} \tag{10.12}$$

**Kostfunctie (gemiddelde kwadratische afwijking):**

$$\text{Kost}(\mathbf{c}) = \frac{1}{N}\sum_{i=1}^N \frac{1}{2}\left\|\mathbf{y}(\mathbf{x}^{(i)}) - F(\mathbf{x}^{(i)};\mathbf{c})\right\|_2^2, \tag{10.13–10.14}$$

waarbij $\mathbf{c}$ alle gewichten en biases omvat. De factor $\frac{1}{2}$ is er voor het gemak bij differentiatie (valt weg), $\frac{1}{N}$ normaliseert zodat de orde van de kost onafhankelijk is van het aantal trainingsdata.

**Het optimalisatieprobleem:**

$$\hat{\mathbf{c}} = \underset{\mathbf{c}}{\arg\min}\ \text{Kost}(\mathbf{c}). \tag{10.6}$$

Dit is een **niet-lineair, niet-convex kleinste-kwadratenprobleem**. De niet-lineariteit komt van de samengestelde activatiefuncties; de niet-convexiteit betekent dat er meerdere lokale minima kunnen zijn.

> ⚠️ **Belangrijk voor examen:** Weet de kostfunctie op te schrijven en te begrijpen. Begrijp waarom dit een niet-lineair optimalisatieprobleem is (de parameters $\mathbf{c}$ zitten niet-lineair in $F$). Weet dat er geen garantie is op het vinden van het globaal minimum.

**Twee fundamentele moeilijkheden:**
1. **Niet-convex:** men kan niet garanderen dat het gevonden minimum globaal is. Het algoritme kan vast zitten in een lokaal minimum.
2. **Hoge dimensionaliteit:** moderne netwerken hebben miljoenen tot miljarden parameters — uitputtend zoeken is onmogelijk.

Dit motiveert Hoofdstuk 11 (optimalisatie-algoritmen), waar methodes zoals gradient descent en toegevoegde gradiënten worden behandeld.

---

## 10.3 Uitbreidingen bij neurale netwerken

### 10.3.1 Convolutionele neurale netwerken (CNNs)

**Het schaalprobleem van volledig verbonden netwerken:**

Een kleurenafbeelding van $200 \times 200$ pixels heeft $200 \times 200 \times 3 = 120\,000$ invoercomponenten. De gewichtsmatrix van de eerste verborgen laag zou dan $n_{\text{out}} \times 120\,000$ elementen hebben — voor $n_{\text{out}} = 1000$ zijn dat al **120 miljoen parameters** in één laag. Dit is onuitvoerbaar.

**Oplossing: gestructureerde (ijle) gewichtsmatrices.**

Een CNN legt een specifiek **sparsiteitspatroon** op aan de gewichtsmatrices. In plaats van elke invoercomponent te verbinden met elk uitvoerneuron, past elk neuron een kleine **lokale filter (kernel)** toe op een *patch* van de invoer.

**Voorbeeld in 1D:** Vermenigvuldig invoervector $\mathbf{x} \in \mathbb{R}^6$ met:

$$W = \begin{pmatrix} 1 & -1 & & & & \\ & 1 & -1 & & & \\ & & 1 & -1 & & \\ & & & 1 & -1 & \\ & & & & 1 & -1 \end{pmatrix} \in \mathbb{R}^{5 \times 6}. \tag{10.15}$$

De $k$-de component van de uitvoer is $x_k - x_{k+1}$ — een **discrete afgeleide**. Dit detecteert abrupte overgangen (randen) in het signaal. In 2D kan men op analoge wijze randen, curven of kleurvlekken detecteren.

**Convolutie:** De bewerking die hierboven beschreven is, heeft de vorm van een **convolutie**. De $k$-de component van de convolutie van $\mathbf{x}$ met filter $g$ is:

$$y_k = \sum_n x_n\, g_{k-n}.$$

Het voorbeeld hierboven: filter $g_0 = 1$, $g_{-1} = -1$, alle andere $g_k = 0$.

**Voordelen van CNNs:**
- **Parameterdeling:** dezelfde filter wordt op elk deel van het beeld toegepast → enorm minder parameters
- **Translatie-invariantie:** een rand op een andere positie in het beeld geeft hetzelfde filterantwoord
- **Efficiënte berekening:** convoluties kunnen ultrasnell berekend worden via FFT-gebaseerde algoritmen

**In de praktijk:** beelddata is een 3D-tensor (breedte × hoogte × kleurkanalen). CNN-filters zijn kleine tensors die op *patches* van die invoertensor worden geschoven (2D-convolutie per kleurkanaal).

> ⚠️ **Belangrijk voor examen:** Begrijp waarom volledig verbonden netwerken niet schalen voor beelden, en hoe CNNs dit oplossen via ijle gewichtsmatrices met een convolutiestructuur. Ken het 1D-voorbeeld.

### 10.3.2 Overfitting vermijden bij neurale netwerken

Neurale netwerken — zeker diepe netwerken met veel parameters — zijn bijzonder gevoelig voor overfitting. Dezelfde principes uit Hoofdstuk 7 zijn van toepassing.

**Geregulariseerde kostfunctie voor een neuraal netwerk:**

$$\text{Kost}_{\text{reg}}(\mathbf{c}) = \frac{1}{N}\sum_{i=1}^N \left\|\mathbf{y}(\mathbf{x}^{(i)}) - \mathbf{a}^{[L]}(\mathbf{x}^{(i)})\right\|_2^2 + \frac{\lambda}{N}\sum_{l=2}^L \left\|W^{[l]}\right\|_2^2. \tag{10.16}$$

- De eerste term: de gewone kostfunctie (datafidelity).
- De tweede term: L2-regularisatie op de gewichtsmatrices (niet op de biases).
- $\lambda > 0$: de regularisatieparameter.

**Motivatie:** Grote gewichten maken neuronen gevoelig — kleine perturbaties in de invoer geven grote schommelingen in de uitvoer. Dat is een teken van overfitting. Door grote gewichten te penaliseren dwing je het netwerk tot robuustere, generaliserendere representaties.

**Praktisch effect op de optimalisatie:** Het toevoegen van de regularisatieterm brengt slechts een **kleine, goedkope wijziging** aan het optimalisatie-algoritme — het is eenvoudig te implementeren.

**Andere technieken tegen overfitting bij neurale netwerken** (ter algemene kennis):
- **Dropout:** tijdens training worden willekeurig neuronen tijdelijk uitgeschakeld, wat het netwerk dwingt robuustere representaties te leren.
- **Early stopping:** het trainingsproces wordt gestopt zodra de generalisatiefout (op testdata) stijgt, ook al daalt de trainingsfout nog.
- **Data augmentatie:** het kunstmatig uitbreiden van de trainingsset (bijv. spiegelen of roteren van afbeeldingen).

---

## Samenvatting

| Concept | Kern |
|---|---|
| Niet-lineair benaderingsprobleem | Basisfuncties hangen zelf af van onbekende parameters → geen lineair stelsel |
| Gevaar van lineariseren | Niet-lineaire transformaties vervormen de foutstructuur |
| Activatiefunctie | De niet-lineaire bouwsteen van een netwerk; sigmoïde: $\sigma(x) = \frac{1}{1+e^{-x}}$, afgeleide: $\sigma'(x) = \sigma(x)(1-\sigma(x))$ |
| Laag van neuronen | $\sigma(W\mathbf{x} + \mathbf{b})$: lineaire transformatie gevolgd door activatiefunctie |
| Diep netwerk | Samengestelde laagsgewijze afbeelding $F(\mathbf{x}; \mathbf{c})$ — niet-lineair in parameters $\mathbf{c}$ |
| Kostfunctie | $\frac{1}{N}\sum_i \frac{1}{2}\|\mathbf{y}^{(i)} - F(\mathbf{x}^{(i)};\mathbf{c})\|_2^2$ |
| Trainen | $\min_\mathbf{c}\ \text{Kost}(\mathbf{c})$ — niet-convex, geen garantie globaal minimum |
| CNN | IJle gewichtsmatrices met convolutiestructuur → schaalbaar voor beelden |
| Regularisatie bij NNs | $+ \frac{\lambda}{N}\sum_l \|W^{[l]}\|_2^2$ toevoegen aan de kostfunctie |

### Verbinding met de rest van de cursus

- **Hoofdstuk 1:** Niet-lineaire benadering als derde categorie naast vectorbenadering en functiebenadering — nu uitgewerkt.
- **Hoofdstuk 7:** Regularisatietechnieken (L2/Tikhonov) worden rechtstreeks toegepast op neurale netwerken in §10.3.2.
- **Hoofdstuk 11:** De optimalisatie-algoritmen (gradient descent, conjugate gradients, Gauss-Newton) die nodig zijn om het netwerk te trainen.
- **Hoofdstuk 12:** De singuliere-waardenontbinding geeft inzicht in lage-rangbenaderingen, wat verwant is aan de structuur die neurale netwerken proberen te leren.
