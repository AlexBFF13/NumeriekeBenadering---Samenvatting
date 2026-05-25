# Deel 7 – Regularisatietechnieken

## Overzicht

In de vorige hoofdstukken werd de **beste benadering** gezocht in een *vooraf gekozen* deelruimte — bijv. veeltermen van graad $n$, of splines met vaste knopen. De keuze van die deelruimte (welke graad? hoeveel knopen?) werd als gegeven beschouwd. Hoofdstuk 7 confronteert ons met de vraag: hoe kies je die deelruimte eigenlijk?

Het antwoord dat dit hoofdstuk geeft is **regularisatie**: in plaats van de complexiteit van het model op voorhand vast te leggen, voeg je een straf toe aan de kostfunctie die te complexe modellen ontmoedigt. Zo kan je zoeken in een grote ruimte (bv. hoge-graadsveeltermen), maar toch een eenvoudig, stabiel model bekomen.

**Rode draad:**
1. Overfitting begrijpen — wanneer treedt het op en hoe herken je het?
2. Detecteren via de generalisatiefout (training- vs. testdata)
3. Vermijden via regularisatie: de gewijzigde kostfunctie
4. Twee concrete smaken: L2-regularisatie (Tikhonov/ridge) en L1-regularisatie (LASSO)
5. Geometrische interpretatie via constrained optimalisatie

---

## 7.1 Inleiding en motivatie

Tot nu toe zochten we de beste benadering

$$y_n(x) = \sum_{k=0}^n c_k \phi_k(x)$$

door de **kostfunctie** te minimaliseren:

$$E(\mathbf{c}) = \sum_{i=1}^N w_i\bigl(f_i - y_n(x_i)\bigr)^2 = \|\mathbf{D}(A\mathbf{c} - \mathbf{f})\|_2^2, \tag{7.3}$$

waarbij $A_{ik} = \phi_k(x_i)$, $\mathbf{f} = (f_i)_i$, en $D$ een diagonaalmatrix met gewichten $\sqrt{w_i}$.

Dit geeft het **normaalstelsel** $(A^TWA)\mathbf{c} = A^TW\mathbf{f}$, met $W = D^2$.

**Het probleem:** de graad $n$ (of meer algemeen: de keuze van de deelruimte) moeten we zelf opgeven. Te laag → slecht model. Te hoog → **overfitting**. Regularisatie is een techniek om dit probleem te omzeilen door de complexiteit van het model te *penaliseren* in de kostfunctie zelf, in plaats van de complexiteit op voorhand hard te begrenzen.

---

## 7.2 Het probleem van overfitting

### Wat is overfitting?

We vertrekken van data $(x_i, f_i)_{i=1}^N$ waarbij de waargenomen waarden **ruizig** zijn:

$$f_i = f(x_i) + \varepsilon_i, \quad \varepsilon_i \sim \mathcal{N}(0, \sigma^2).$$

De ware onderliggende functie $f(x)$ is glad; de metingen $f_i$ bevatten echter toevallige schommelingen door $\varepsilon_i$.

**Overfitting** treedt op wanneer het model te nauw aansluit bij de ruizige data in plaats van bij de onderliggende functie:
- het model leert de **ruis** in de trainingsdata in plaats van de **structuur** van $f$.
- het model presteert goed op de trainingspunten, maar slecht op nieuwe punten.

### Illustratie: $f(x) = \sin(2\pi x)$ met $N = 11$ ruizige metingen, $\sigma = 0.1$

| Graad $n$ | Gedrag |
|---|---|
| 1 | **Onderfitting (oversmoothing):** het model is te simpel om de oscillatie te vatten |
| 3 | Goede balans — volgt de onderliggende structuur goed |
| 6 | Ogenschijnlijk ook goed — maar begint wat te reageren op de ruis |
| 9 | **Overfitting:** de veelterm oscilleert hevig en past perfect op de ruizige data, maar beschrijft $f$ slecht |

> ⚠️ **Belangrijk voor examen:** Begrijp wat overfitting is, wanneer het optreedt, en hoe je het herkent. Ken het onderscheid tussen overfitting en onderfitting (oversmoothing).

**Kernobservatie:** Wanneer $n = N - 1$ is er een unieke interpolerende veelterm. Die reproduceert de data exact — inclusief alle meetfouten. Hogere-graadsveeltermen bevatten genoeg vrijheidsgraden om de ruis te leren in plaats van de functie.

---

## 7.3 Detecteren van overfitting: de generalisatiefout

### Waarom visuele inspectie niet altijd werkt

In hoge dimensies (bijv. foto's als vectoren in $\mathbb{R}^D$ met $D$ groot) kan je de benadering niet visualiseren. En in de praktijk ken je de ware functie $f$ niet. Je kunt overfitting dus niet altijd met het blote oog vaststellen.

### Trainings- vs. testdata

De oplossing: gebruik een **tweede dataset** om de kwaliteit van de benadering objectief te meten.

| Naam | Definitie | Gebruik |
|---|---|---|
| **Trainingsdata** | $(x_i, f_i)_{i=1}^N$ | Opstellen van de benadering (minimaliseren van $E_{\text{train}}$) |
| **Testdata** | $(\tilde{x}_i, \tilde{f}_i)_{i=1}^{\tilde{N}}$ | Evalueren van de benadering op ongeziene data |

**Trainingsfout:**
$$E_{\text{train}} = \sum_{i=1}^N w_i\bigl(f_i - y_n(x_i)\bigr)^2. \tag{7.5}$$

**Generalisatiefout (testfout):**
$$E_{\text{test}} = \sum_{i=1}^{\tilde{N}} \tilde{w}_i\bigl(\tilde{f}_i - y_n(\tilde{x}_i)\bigr)^2. \tag{7.6}$$

### Gedrag als functie van de graad

> ⚠️ **Belangrijk voor examen:** Ken dit gedrag en kun je het verklaren.

- De **trainingsfout** daalt **monotoon** als de graad stijgt. Dit volgt wiskundig uit het feit dat de beste benadering in een grotere deelruimte altijd minstens zo goed is als in een kleinere.
- De **generalisatiefout** daalt aanvankelijk mee, maar **begint opnieuw te stijgen** bij te hoge graden — het signaal van overfitting.

De optimale graad is die waarbij de **generalisatiefout minimaal** is.

### Opmerking: bias-variance trade-off

De trainingsfout is een deterministische grootheid (de **bias**: hoe goed het model de trainingsdata beschrijft). De generalisatiefout is een toevalsvariabele (de **variantie**: hoeveel de benadering schommelt afhankelijk van de testpunten). Bij overfitting is de variantie groot. De keuze voor model-complexiteit is dan ook een afweging: hoge bias (te eenvoudig model) vs. hoge variantie (te complex model).

### Alternatief: leave-one-out kruisvalidatie

Een manier om testdata te vermijden: stel $N$ benaderingen op, elke keer met één datapunt weggelaten, en kijk hoe goed elk model presteert op het weggelaten punt. Dit is rekenkundig duur maar statistisch robuust.

---

## 7.4 Overfitting vermijden: principe van regularisatie

### Het basisidee

Een hogere-graadsveelterm *kan* elke lagere-graadsveelterm nabootsen door de hogere coëfficiënten op nul te zetten. Het probleem is: we weten niet welke coëfficiënten dat moeten zijn.

**Idee van regularisatie:** in plaats van coëfficiënten hard op nul te zetten, **penaliseren** we grote coëfficiënten in de kostfunctie. Dit geeft een balans tussen:
- de data goed beschrijven (lage trainingsfout)
- een eenvoudig model houden (kleine coëfficiënten)

### De geregulariseerde kostfunctie

$$E_{\text{reg}}(\mathbf{c}) = \underbrace{\sum_{i=1}^N w_i\bigl(f_i - y_n(x_i)\bigr)^2}_{\text{datafidelity}} + \lambda \underbrace{F(\mathbf{c})}_{\text{regularisatieterm}}, \tag{7.8}$$

waarbij:
- $\lambda > 0$ de **regularisatieparameter** (hyperparameter) is — bepaalt de verhouding tussen datafidelity en regularisatie
- $F(\mathbf{c})$ de **penalisatiefunctie** is — meet de "complexiteit" van de coëfficiënten

### De rol van $\lambda$

| $\lambda$ | Effect |
|---|---|
| $\lambda = 0$ | Geen regularisatie — gewone kleinste-kwadratenoplossing |
| $\lambda$ klein | Kleine invloed van regularisatie — coëfficiënten krimpen lichtjes |
| $\lambda$ groot | Sterke regularisatie — coëfficiënten worden sterk naar nul gedrukt |
| $\lambda \to \infty$ | Alle coëfficiënten gaan naar nul: $\mathbf{c} \to \mathbf{0}$ |

De "goede" waarde van $\lambda$ wordt gekozen door de generalisatiefout te vergelijken voor verschillende waarden van $\lambda$.

### Twee keuzes voor $F(\mathbf{c})$

> ⚠️ **Belangrijk voor examen:** Ken de twee vormen en hun eigenschappen.

**L2-penalisatie (2-norm):**
$$F(\mathbf{c}) = \|\mathbf{c}\|_2^2 = \sum_{k=0}^n c_k^2 \quad \Rightarrow \quad \text{Tikhonovregularisatie (ridge regression)}. \tag{7.9}$$

**L1-penalisatie (1-norm):**
$$F(\mathbf{c}) = \|\mathbf{c}\|_1 = \sum_{k=0}^n |c_k| \quad \Rightarrow \quad \text{LASSO-regularisatie}. \tag{7.10}$$

### Belangrijk: basis-afhankelijkheid

Anders dan bij de gewone beste benadering (die uniek is ongeacht de keuze van basis), **hangt de geregulariseerde oplossing af van de gekozen basis**. De reden: de waarde van $F(\mathbf{c})$ hangt af van de grootte van de coëfficiënten, en die verschilt als je een andere basis gebruikt om dezelfde functie voor te stellen. Dit is een fundamenteel verschil t.o.v. het ongeregulariseerde benaderingsprobleem.

---

## 7.5 Tikhonovregularisatie (L2-regularisatie / ridge regression)

### Definitie

$$E_{\text{reg}}(\mathbf{c}) = \|\mathbf{D}(A\mathbf{c} - \mathbf{f})\|_2^2 + \lambda\|\mathbf{c}\|_2^2 = (A\mathbf{c} - \mathbf{f})^T W (A\mathbf{c} - \mathbf{f}) + \lambda \mathbf{c}^T\mathbf{c}. \tag{7.11}$$

### Afleiding van het geregulariseerde normaalstelsel

We minimaliseren door de gradiënt naar $\mathbf{c}$ op nul te stellen:

$$\nabla_\mathbf{c} E_{\text{reg}}(\mathbf{c}) = 0.$$

**Stap 1:** Schrijf $E_{\text{reg}}$ uit:

$$E_{\text{reg}}(\mathbf{c}) = \mathbf{c}^T(A^TWA + \lambda I)\mathbf{c} - 2\mathbf{f}^TWA\mathbf{c} + \mathbf{f}^TW\mathbf{f}. \tag{7.12}$$

→ We gebruiken dat $\mathbf{f}^TWA\mathbf{c} = \mathbf{c}^TA^TW\mathbf{f}$ (reëel geval).

**Stap 2:** Differentieer naar $\mathbf{c}$ en stel gelijk aan nul:

$$\nabla_\mathbf{c} E_{\text{reg}} = 2(A^TWA + \lambda I)\mathbf{c} - 2A^TW\mathbf{f} = \mathbf{0}.$$

**Conclusie:** Het **geregulariseerde normaalstelsel** is:

$$\boxed{(A^TWA + \lambda I)\,\hat{\mathbf{c}} = A^TW\mathbf{f}.} \tag{7.14}$$

> ⚠️ **Belangrijk voor examen:** Weet dit stelsel af te leiden en het te vergelijken met het gewone normaalstelsel $A^TWA\mathbf{c} = A^TW\mathbf{f}$.

### Vergelijking met het gewone normaalstelsel

| | Normaalstelsel | Geregulariseerd normaalstelsel |
|---|---|---|
| Matrix | $A^TWA$ | $A^TWA + \lambda I$ |
| Altijd inverteerbaar? | Nee (als $A$ niet volle rang heeft) | Ja (voor $\lambda > 0$) |
| Conditie | $\kappa(A)^2$ | Beter: $\lambda$ verhoogt de kleinste eigenwaarden |

**Cruciale eigenschap:** Door $\lambda I$ op te tellen, worden alle eigenwaarden van de systeemmatrix minstens $\lambda$. Dit maakt het stelsel **altijd goed gesteld** (inverteerbaar), zelfs als $A^TWA$ singulier of bijna singulier is. Regularisatie lost dus ook het probleem van slecht geconditioneerde normaalsystemen op.

**Praktisch voordeel:** Hetzelfde type lineair stelsel als het normaalstelsel → zelfde algoritmen (Cholesky-factorisatie, QR, ...) zijn toepasbaar.

### Veralgemeende Tikhonovregularisatie

In de meest algemene vorm wordt de penalisatie per coëfficiënt gewogen:

$$E_{\text{reg,Tikh}} = \|\mathbf{D}(A\mathbf{c} - \mathbf{f})\|_2^2 + \|\Gamma\mathbf{c}\|_2^2, \tag{7.15}$$

waarbij $\Gamma$ een diagonaalmatrix is met gewichten $(\gamma_k)_k$. Dit laat toe bepaalde coëfficiënten harder te penaliseren dan andere. Het nadeel is dat er meer hyperparameters gekozen moeten worden (de hele matrix $\Gamma$).

**Terminologie:**
- **Ridge regression** = Tikhonov met scalaire $\lambda$ (vereenvoudigd)
- **Tikhonovregularisatie** = de veralgemeende vorm met matrix $\Gamma$

In de praktijk worden deze termen vaak door elkaar gebruikt.

---

## 7.6 LASSO-regularisatie (L1-regularisatie)

### Motivatie

Tikhonovregularisatie maakt coëfficiënten **klein**, maar zelden **exact nul**. De reden is wiskundig: de $\ell^2$-penalisatie $\|\mathbf{c}\|_2^2$ is differentieerbaar en duwt coëfficiënten monotoon naar nul, maar raakt nooit precies nul (tenzij $\lambda \to \infty$).

**LASSO** (Least Absolute Shrinkage and Selection Operator) gebruikt de $\ell^1$-norm als penalisatie, die wel in staat is om coëfficiënten exact op nul te zetten:

$$E_{\text{reg,L1}}(\mathbf{c}) = \sum_{i=1}^N w_i\bigl(f_i - y_n(x_i)\bigr)^2 + \lambda\sum_{k=0}^n |c_k|. \tag{7.16}$$

### Gevolg voor de oplossing

De $\ell^1$-norm is **niet differentieerbaar** in de punten $c_k = 0$. Dit heeft twee gevolgen:

1. De optimaliteitsvoorwaarde $\nabla_\mathbf{c} E_{\text{reg,L1}} = \mathbf{0}$ is **geen lineair stelsel meer** → andere optimalisatiemethodes zijn nodig (bijv. coördinaat-afdaling, subgradiëntmethodes).

2. Precies die niet-differentieerbaarheid in $c_k = 0$ zorgt ervoor dat LASSO de neiging heeft coëfficiënten **exact op nul** te stellen — wat leidt tot **ijle (sparse) oplossingen**.

> ⚠️ **Belangrijk voor examen:** Weet waarom LASSO schaarse oplossingen geeft, en waarom Tikhonov dat niet doet. Het verschil zit in de geometrie van de respectieve eenheidsbollen (zie §7.7).

---

## 7.7 Regularisatie als optimalisatie met beperkingen

### Equivalente formulering

Het regularisatieprobleem

$$\min_\mathbf{c}\ E(\mathbf{c}) + \lambda F(\mathbf{c})$$

is equivalent met het **geconstrained optimalisatieprobleem**:

$$\min_\mathbf{c}\ E(\mathbf{c}) \quad \text{subject to}\quad F(\mathbf{c}) \leq K, \tag{7.17–7.18}$$

voor een passende waarde van $K$ die overeenkomt met de gekozen $\lambda$ (via de Lagrange-dualiteit). De waarde van $K$ is kleiner naarmate $\lambda$ groter is.

**Intuïtie:** $F(\mathbf{c}) \leq K$ beperkt de zoekruimte tot de $K$-bol in de respectieve norm. We zoeken het punt in die bol dat de ongeregulariseerde kostfunctie $E(\mathbf{c})$ minimaliseert.

### Geometrische interpretatie: waarom LASSO schaars is

Beschouw het geval $n = 2$ (twee coëfficiënten $c_1, c_2$):

**Ridge regression ($\ell^2$-bol):**
- De toegestane set $\{\mathbf{c} : c_1^2 + c_2^2 \leq K\}$ is een **schijf** (gladde, ronde rand).
- De contourlijnen van $E(\mathbf{c})$ (ellipsen) raken de bol typisch op een punt **niet** op een as.
- Gevolg: de optimale $\mathbf{c}^*$ heeft typisch **beide coëfficiënten niet-nul**.

**LASSO ($\ell^1$-bol):**
- De toegestane set $\{\mathbf{c} : |c_1| + |c_2| \leq K\}$ is een **ruit** (met scherpe hoeken op de assen).
- De contourlijnen van $E(\mathbf{c})$ raken de ruit het vaakst in een **hoekpunt**, d.w.z. precies op een as.
- Gevolg: de optimale $\mathbf{c}^*$ heeft typisch **één coëfficiënt exact nul** → **schaarse oplossing**.

**Intuïtie in woorden:** De hoeken van de $\ell^1$-bol liggen op de coördinaatassen. Elliptische niveaukrommen die "van buiten naar binnen" bewegen, raken de ruitvorm het vaakst in een hoek — op een as. Op een as is één coëfficiënt nul. De gladde $\ell^2$-bol heeft geen hoeken, dus de raakpunten liggen willekeurig op de rand.

```
       L2-bol (ridge)           L1-bol (LASSO)

           *                        /\
         *   *                     /  \
        *  c* *     c* toevallig  * c* \   c* raakt hoek
         *   *      op rand         \  /   → c_2 = 0
           *                        \/
```

> ⚠️ **Belangrijk voor examen:** Dit geometrische argument is het kernargument waarom LASSO schaarse oplossingen geeft. Ken het en kun je het uitleggen aan de hand van een schets.

---

## Samenvatting

| Concept | Formule / resultaat |
|---|---|
| Gewone kostfunctie | $E(\mathbf{c}) = \|\mathbf{D}(A\mathbf{c} - \mathbf{f})\|_2^2$ |
| Trainingsfout | $E_{\text{train}} = \sum_i w_i(f_i - y_n(x_i))^2$ |
| Generalisatiefout | $E_{\text{test}} = \sum_i \tilde{w}_i(\tilde{f}_i - y_n(\tilde{x}_i))^2$ |
| Geregulariseerde kostfunctie | $E_{\text{reg}}(\mathbf{c}) = E(\mathbf{c}) + \lambda F(\mathbf{c})$ |
| **Tikhonov/ridge** ($F = \|\cdot\|_2^2$) | $(A^TWA + \lambda I)\hat{\mathbf{c}} = A^TW\mathbf{f}$ |
| **LASSO** ($F = \|\cdot\|_1$) | Niet-lineair stelsel, vereist iteratieve methodes |
| Regularisatie als constraint | $\min E(\mathbf{c})$ s.t. $F(\mathbf{c}) \leq K$ |
| Waarom LASSO schaars? | $\ell^1$-bol heeft hoeken op de assen → raakpunten liggen op assen |
| Waarom Tikhonov altijd inverteerbaar? | $\lambda > 0$ verhoogt alle eigenwaarden: $A^TWA + \lambda I$ is s.p.d. |

### Keuze van $\lambda$: praktische aanpak

1. Kies een reeks van waarden $\lambda_1 < \lambda_2 < \cdots < \lambda_K$.
2. Stel voor elke $\lambda$ de geregulariseerde oplossing $\hat{\mathbf{c}}(\lambda)$ op.
3. Evalueer de **generalisatiefout** op testdata voor elk $\lambda$.
4. Kies de $\lambda$ die de generalisatiefout minimaliseert.

### Verbinding met de rest van de cursus

- De geregulariseerde kostfunctie $(7.8)$ is een speciaal geval van de algemene formulering uit §1.5.2 (zie formule 1.19).
- Tikhonovregularisatie leidt tot een lineair stelsel dat nauw verwant is aan het normaalstelsel uit §2.4.1 — hetzelfde type algoritmen zijn toepasbaar.
- Regularisatie bij neurale netwerken (Hoofdstuk 10) volgt dezelfde principes, maar in een niet-lineaire context.
