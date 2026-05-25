# Hoofdstuk1 – Inleiding

## Overzicht

Dit hoofdstuk legt het fundament voor de hele cursus. Het antwoordt op de vraag: *wat is een benaderingsprobleem eigenlijk, en waarom doen we dat?* De vier bouwstenen die hier geïntroduceerd worden — de te benaderen functie, de klasse van benaderingsfuncties, het benaderingscriterium, en het benaderingsalgoritme — keren in elk volgend hoofdstuk terug. Dit hoofdstuk bevat weinig wiskundige bewijzen; het gaat eerder over het opbouwen van intuïtie en het kaderen van het hele vakgebied.

---

## 1.1 Context: numerieke analyse in de ingenieurspraktijk

De rol van numerieke analyse begrijp je het best als je de stappen beschouwt die een ingenieur doorloopt van idee tot product:

1. **Beschrijf het fysisch probleem** (bijv. een windturbine ontwerpen)
2. **Stel een wiskundig model op** (differentiaalvergelijkingen, fysische wetten)
3. **Stel een numeriek model op** — één dat op een computer uitgerekend kan worden
4. **Ontwikkel numerieke methodes** om dat model op te lossen
5. **Implementeer de methodes** in bruikbare software

*Deze cursus focust op stappen 3 en 4.* Meer concreet behandelt ze twee samenhangende onderwerpen:

- **Benaderingstheorie** — hoe stel je een goede benadering op voor een functie (gekend of ongekend, afgeleid uit data)?
- **Numerieke lineaire algebra** — hoe bereken je die benadering efficiënt en nauwkeurig?

Het sleutelinzicht van de cursus is dat heel verschillende benaderingsproblemen — van tijdreeksen over beeldverwerking tot differentiaalvergelijkingen — op **dezelfde onderliggende wiskundige principes** berusten. Door die principes te begrijpen, begrijp je een heel breed spectrum aan methodes in één keer.

---

## 1.2 De componenten van een benaderingsprobleem

Elk benaderingsprobleem heeft vijf vragen waarop je een antwoord moet geven:

| Vraag | Inhoud |
|---|---|
| **Wat** benaderen we? | De aard van de data: continu of discreet, 1D of hoger, exact of met ruis |
| **Waarom** benaderen we? | Het doel: compressie, ruisverwijdering, interpolatie, classificatie, ... |
| **Waarmee** benaderen we? | De klasse van benaderingsfuncties: veeltermen, splines, trigonometrische functies, ... |
| **Hoe beoordelen** we de kwaliteit? | Het benaderingscriterium: minimax, kleinste kwadraten, interpolatie, ... |
| **Hoe berekenen** we de beste benadering? | Het algoritme en de implementatie |

De eerste twee vragen komen in dit hoofdstuk aan bod. De rest vormt het onderwerp van de rest van de cursus.

### 1.2.1 Wat wensen we te benaderen?

Data kan **continu** zijn (een functie van tijd of ruimte) of **discreet** (opgemeten op een eindig aantal punten). In de praktijk is data vrijwel altijd discreet — denk aan foto's, sensormetingen, tijdreeksen. Toch is het nuttig om die data te modelleren als een continue onderliggende functie, omdat dat toelaat na te denken over het te benaderen object los van de meetresolutie.

Data kan ook geclassificeerd worden naar **dimensie**: een tijdreeks is 1D, een beeld is 2D, een video is 3D. Die dimensie beïnvloedt welke benaderingstechniek gepast is.

Wiskundig stelt men data altijd voor als een **vector** in een vectorruimte. Een essentiële vraag is dan: bevat die vector exact de ware waarden, of zit er meetruis op?

### 1.2.2 Waarvoor dient de benadering?

Er zijn meerdere typische doelstellingen:

- **Datacompressie**: de data voorstellen met minder coëfficiënten (en eventueel reconstrueren vanuit die compressie).
- **Ruisverwijdering**: de onderliggende gladde functie terugvinden achter ruizige metingen.
- **Data-interpretatie**: betekenis extraheren, bijv. een aankoopaanbeveling als lineaire combinatie van gebruikersprofielen, of beeldclassificatie (kat/geen kat).
- **Data-manipulatie**: de afgeleide of integraal berekenen, interpoleren, extrapoleren (weersvoorspelling, risicoschatting in verzekeringen).

> ⚠️ **Belangrijk voor examen:** Ken de vier doelstellingen en weet ze te illustreren met een voorbeeld. De keuze van benaderingscriterium en -klasse hangt direct af van het doel.

---

## 1.3 Van data naar functiebenadering

### 1.3.1 Benadering van vectoren in eindigdimensionale ruimtes

Beschouw een signaalvector $\mathbf{b} = (b_1, \ldots, b_m) \in \mathbb{R}^m$ en een set referentiesignalen $\mathbf{a}_1, \ldots, \mathbf{a}_n \in \mathbb{R}^m$ (met $n < m$). De vraag is: hoe goed kan $\mathbf{b}$ voorgesteld worden als lineaire combinatie van de $\mathbf{a}_i$?

Wiskundig: zoek de **beste benadering** van $\mathbf{b}$ in de deelruimte

$$\mathcal{D} = \langle \mathbf{a}_1, \ldots, \mathbf{a}_n \rangle \subset \mathbb{R}^m,$$

d.w.z. vind coëfficiënten $x_1, \ldots, x_n$ zodat

$$\mathbf{b} \approx \sum_{i=1}^n x_i \, \mathbf{a}_i.$$

De vector $\mathbf{x} = (x_1, \ldots, x_n)$ is dan een **laagdimensionale representatie** van $\mathbf{b}$. Dit is het onderwerp van **Hoofdstuk 2**.

### 1.3.2 Lineaire benadering van functies

Het equivalent voor continue functies: benader een gegeven functie $f(x)$ door een eindige lineaire combinatie van basisfuncties $\phi_k(x)$:

$$y_n(x) = \sum_{k=0}^n a_k \, \phi_k(x). \tag{1.1}$$

De keuze van de basisfuncties is een aparte wetenschap:

| Situatie | Geschikte basis |
|---|---|
| Periodieke functies | $1, \cos(kx), \sin(kx)$ — Fourier-reeks |
| Functie met verval-/hyperbolisch karakter | Negatieve machten van $x$: $1, x^{-1}, x^{-2}, \ldots$ |
| Functie met exponentieel verval (radioactiviteit) | $e^{-\lambda_k t}$, waarbij $\lambda_k$ de vervalkonstanten zijn |
| Algemeen gladde functie | Veeltermen (Hoofdstuk 4), splines (Hoofdstuk 5) |

Een derde reden om andere basisfuncties te kiezen is **numeriek van aard**: sommige basissen leiden tot slecht geconditioneerde stelsels, waarbij afrondingsfouten zich sterk voortplanten. *Orthogonale* basissen vermijden dit probleem.

### 1.3.3 Niet-lineaire benaderingen

Als de basisparams $\lambda_k$ in (1.5) *niet* op voorhand gekend zijn, heb je een **niet-lineair** benaderingsprobleem. Typisch voorbeeld: radioactief verval waarbij je noch de vervalkonstanten, noch de coëfficiënten kent. Deze problemen zijn dikwijls slecht geconditioneerd.

**Rationale benaderingen** (teller/noemer beide veeltermen) zijn niet-lineair en nuttig voor functies met singuliere punten:

$$y_{m,n}(x) = \frac{a_0 + a_1 x + \cdots + a_m x^m}{b_0 + b_1 x + \cdots + b_n x^n}. \tag{1.6}$$

**Diepe neurale netwerken** zijn ook niet-lineaire benaderingstechnieken — ze worden behandeld in Deel III van de cursus.

---

## 1.4 Veeltermbenadering als illustratief voorbeeld

Dit onderdeel illustreert de vier componenten van een benaderingsprobleem aan de hand van een concreet geval: benader $e^x$ op $[-1, 1]$ door een veelterm van graad 4.

### 1.4.1 Benadering van een continue functie

#### Maclaurin-reeks (afgebroken Taylor-reeks)

De Maclaurin-reeks van $e^x$ afbreken na de vijfde term geeft:

$$y_4(x) = 1 + x + \frac{1}{2}x^2 + \frac{1}{6}x^3 + \frac{1}{24}x^4.$$

**Gedrag van de fout $r(x) = e^x - y_4(x)$:**

- De fout is nul in $x = 0$ (de ontwikkelingspunt) en de eerste vier afgeleiden van de benadering zijn er gelijk aan die van $e^x$.
- De fout **stijgt monotoon** naarmate $|x|$ toeneemt.
- De maximale fout treedt op aan de intervalgrens $x = 1$: fout $\approx 10^{-2}$.

**Intuïtie:** Een Taylor-reeks is lokaal optimaal. Ze is extreem nauwkeurig dicht bij het ontwikkelingspunt, maar doet het slechter aan de randen. Dit betekent dat een groot gedeelte van de benaderingsruimte "verspild" wordt aan extra precisie in het midden, terwijl de randen te weinig aandacht krijgen. *Een meer gelijkmatige foutsverdeling is wenselijk.*

#### Veelterminterpolatie (Lagrange)

In plaats van de reeks, vertrek je vanuit vijf tabelwaarden $(x_i, e^{x_i})$ voor equidistante punten $x_i \in [-1,1]$. De **interpolerende veelterm** $y_{N-1}(x)$ is de unieke veelterm van graad $N-1$ die de functiewaarden exact reproduceert:

$$y_{N-1}(x_i) = f_i \quad \text{voor } i = 1, \ldots, N. \tag{1.7}$$

**Lagrange-interpolatieformule:**

$$y_{N-1}(x) = \sum_{k=1}^N l_k(x) \, f_k, \tag{1.8}$$

met de **Lagrange-basisfuncties**:

$$l_k(x) = \prod_{\substack{j=1 \\ j \neq k}}^N \frac{x - x_j}{x_k - x_j} = \frac{\omega(x)}{(x - x_k)\,\omega'(x_k)}, \tag{1.9}$$

waarbij $\omega(x) = \prod_{j=1}^N (x - x_j)$.

De functie $l_k(x)$ is de unieke veelterm van graad $N-1$ met $l_k(x_k) = 1$ en $l_k(x_j) = 0$ voor $j \neq k$.

**Resultaat:** De maximale fout daalt van $\approx 10^{-2}$ (Maclaurin) naar $\approx 10^{-3}$ (interpolatie met equidistante punten), bij dezelfde graad 4. Maar equidistante punten zijn *niet* optimaal — er bestaan betere keuzes (Chebyshev-punten, zie Hoofdstuk 4).

> ⚠️ **Belangrijk voor examen:** Ken de Lagrange-formule en weet waarom interpolatie beter is dan Taylor voor globale benadering. Weet ook dat equidistante knooppunten niet optimaal zijn.

#### Andere benaderingscriteria

Er zijn twee fundamentele criteria die *alle* functiewaarden over het interval in rekening brengen, in plaats van slechts $N$ punten:

**Minimaxcriterium (Chebyshev-benadering):**

$$\min_{a_0, \ldots, a_4} \max_{x \in [-1,1]} \left| e^x - (a_0 + a_1 x + \cdots + a_4 x^4) \right|. \tag{1.12}$$

Dit minimaliseert de grootste afwijking. Het probleem is moeilijk op te lossen: het is een minimalisatieprobleem in een vijfdimensionale ruimte.

**Kleinste-kwadratencriterium:**

$$\min_{a_0, \ldots, a_4} \int_{-1}^1 w(x) \left[ e^x - (a_0 + a_1 x + \cdots + a_4 x^4) \right]^2 dx. \tag{1.13}$$

Hier is $w(x) > 0$ een **gewichtsfunctie** die bepaalt hoeveel nadruk er op elk punt van het interval wordt gelegd. Dit probleem is *veel* eenvoudiger op te lossen dan het minimax-probleem.

**Effect van de gewichtsfunctie:** Met $w(x) \equiv 1$ is de maximale fout $\approx 1{,}2 \times 10^{-3}$; met de **Chebyshev-gewichtsfunctie**

$$w(x) = \frac{1}{\sqrt{1-x^2}} \tag{1.15}$$

(die extra gewicht geeft aan de randen) daalt de maximale fout naar $\approx 5{,}9 \times 10^{-4}$ en is de fout **gelijkmatig verspreid** over het interval. Dit is de wenselijke situatie.

> ⚠️ **Belangrijk voor examen:** Weet de twee criteria te formuleren (minimax vs. kleinste kwadraten), weet dat de gewichtsfunctie de kwaliteit van de benadering beïnvloedt, en begrijp waarom kleinste kwadraten in de praktijk populairder is (eenvoudiger op te lossen).

### 1.4.2 Benadering van een discrete functie

Wanneer de data niet een continue functie is maar een tabel $(x_i, f_i)_{i=1}^N$ — bijv. opgemeten signalen met ruis — dan is interpolatie van hoge graad een slechte keuze. Er zijn vier redenen:

1. **Ruis volgen:** Een interpolant gaat exact door elk punt, inclusief de meetfout. Dat is precies wat je *niet* wilt bij ruizige data.
2. **Onvoldoende nauwkeurigheid in de data:** Als data slechts op 5 cijfers nauwkeurig is, heeft het geen zin een benadering te vragen die nauwkeuriger is.
3. **Datacompressie:** Je wil typisch een compacte functie met veel minder parameters dan datapunten.
4. **Numerieke problemen:** Veeltermen van hoge graad zijn slecht geconditioneerd.

Het alternatief is de **discrete kleinste-kwadratenbenadering**:

$$\min_{a_0, \ldots, a_n} \sum_{i=1}^N w_i \left[ f_i - \sum_{k=0}^n a_k x_i^k \right]^2. \tag{1.17}$$

Dit geeft een zacht verlopende benadering die de ruis "wegfiltert" (*curve fitting*). De gewichten $w_i$ laten toe sommige datapunten meer te laten doorwegen.

---

## 1.5 De beste benadering als optimalisatieprobleem

### Algemene formulering

> ⚠️ **Belangrijk voor examen:** Zorg dat je de algemene optimalisatieformulering begrijpt — ze is de rode draad van de hele cursus.

Gegeven data $(x_i, f_i)_{i=1}^N$ en basisfuncties $\{\phi_k\}_{k=1}^n$, zoek de coëfficiënten die

$$\underset{\{a_k\}}{\arg\min} \left\| \mathbf{f} - \sum_{k=1}^n a_k \boldsymbol{\phi}_k \right\| \tag{1.18}$$

minimaliseren, waarbij $\mathbf{f} = (f_i)_i$ en $\boldsymbol{\phi}_k = (\phi_k(x_i))_i \in \mathbb{R}^N$.

De keuze van de **norm** $\|\cdot\|$ in (1.18) bepaalt het soort benadering:
- $\ell^\infty$-norm (maximumafwijking) → minimaxcriterium
- $\ell^2$-norm (som van kwadraten) → kleinste-kwadraatcriterium

In Hoofdstuk 3 behandelen we dit in een abstracte, algemene setting.

### 1.5.1 Overfitting en oversmoothing

De kwaliteit van de beste benadering hangt sterk af van de **verhouding** tussen:
- het aantal te bepalen parameters $n$ (de dimensie van het model)
- het aantal datapunten $N$
- de ruis op de data

**Voorbeeld:** Zoek een benadering voor $f(x) = \sin(2\pi x)$ op $[0,1]$ gegeven 11 ruizige meetpunten met standaardafwijking $0{,}2$.

- **Graad 10** (= interpolerend, 11 parameters voor 11 punten): de veelterm reproduceert de data exact, inclusief de ruis. Ze oscilleert wild. Dit is **overfitting** — het model past zich aan aan de ruis in plaats van aan de onderliggende functie.
- **Graad 1** (lineair, 2 parameters): de benadering is veel te glad en mist de oscillerende structuur van $\sin(2\pi x)$. Dit is **oversmoothing** — het model is te simpel.
- **Graad 5**: een goed evenwicht — de benadering volgt de onderliggende functie goed zonder de ruis te volgen.

**Intuïtie:** Meer parameters geven meer flexibiliteit, maar ook meer kans dat je ruis gaat fitten. Minder parameters is robuuster voor ruis, maar kan de werkelijke structuur missen. Het vinden van de juiste balans is een centraal thema in de cursus.

> ⚠️ **Belangrijk voor examen:** Weet de begrippen *overfitting* en *oversmoothing* te definiëren en illustreren. Begrijp hoe de dimensie van het model ten opzichte van de data de kwaliteit van de benadering bepaalt.

### 1.5.2 Regularisatie

In de praktijk is de ware onderliggende functie ongekend, waardoor je niet meteen weet of je met overfitting of oversmoothing te maken hebt. Je wilt het model **zo eenvoudig mogelijk houden** maar **niet eenvoudiger dan nodig**.

Bij neurale netwerken is het probleem zelfs extremer: het aantal parameters kan (veel) groter zijn dan het aantal datapunten. In dat geval is er niet eens een unieke beste benadering.

Al deze gevallen vallen onder de noemer **slecht gestelde problemen** (*ill-posed problems*). De oplossing is **regularisatie**: een straf op de complexiteit van het model toevoegen aan het optimalisatiecriterium. Een typische geregulariseerde objectieffunctie heeft de vorm:

$$\underset{\{a_k\}}{\arg\min} \left\| \mathbf{f} - \sum_{k=1}^n a_k \boldsymbol{\phi}_k \right\|_\alpha + \lambda \|\mathbf{a}\|_\beta. \tag{1.19}$$

- De eerste term: **datafidelity** — de benadering moet de data goed reproduceren.
- De tweede term: **regularisatieterm** — de benadering mag niet te complex worden.
- $\lambda > 0$: de **regularisatieparameter** die de balans bepaalt tussen de twee termen.

Regularisatie wordt uitgewerkt in **Hoofdstukken 7 en 10**.

---

## 1.6 Besluit: de vier componenten van een benaderingsprobleem

> ⚠️ **Belangrijk voor examen:** Elk benaderingsprobleem heeft precies deze vier componenten. Zorg dat je ze kan benoemen, beschrijven en illustreren.

Elk benaderingsprobleem bevat steeds vier onlosmakelijk verbonden componenten:

1. **De te benaderen functie (of data)**  
   Kan continu of discreet zijn, exact of ruizig, één- of meerdimensionaal. Soms is de "ware" te benaderen functie zelfs ongekend (bijv. bij beeldrestauratie).

2. **De klasse van benaderingsfuncties**  
   Wordt gekozen op basis van:
   - Wiskundige eenvoud (veeltermen, splines)
   - Fysische motivatie (trigonometrisch, exponentieel)
   - Numerieke eigenschappen (orthogonale basisfuncties)
   - Compactheid van representatie (wavelets)

3. **Het benaderingscriterium**  
   Bepaalt wat "beste benadering" betekent. Keuze hangt af van het doel:
   - Interpolatie (exact door datapunten)
   - Kleinste kwadraten (gemiddelde fout minimaliseren)
   - Minimax (maximale fout minimaliseren)

4. **Het benaderingsalgoritme**  
   Hoe berekenen we de optimale coëfficiënten efficiënt en numeriek stabiel? Dit leidt tot computerprogramma's die de benadering genereren.

---

## Verbinding met de rest van de cursus

| Hoofdstuk | Onderwerp | Verband met Deel 1 |
|---|---|---|
| 2 | Vectorbenadering | Beste benadering van $\mathbf{b}$ in een deelruimte $\mathcal{D}$ |
| 3 | Functiebenadering (abstract) | Abstracte normen en criteria; unitaire ruimten |
| 4 | Veeltermbenadering | Concrete uitwerking van (1.1) met orthogonale polynomen |
| 5 | Splinebenadering | Stuksgewijze veeltermbenaderingen |
| 6 | Discrete benadering uit data | Concrete uitwerking van (1.17) |
| 7 | Regularisatie | Concrete uitwerking van (1.19) |
| 10 | Neurale netwerken | Niet-lineaire benaderingen, overfitting bij groot parameteraantal |
