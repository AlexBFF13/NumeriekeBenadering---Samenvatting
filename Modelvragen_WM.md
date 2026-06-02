# Modelexamen Wim Michiels — Uitgewerkte antwoorden

> KU Leuven, Numerieke Benadering 2025–2026.
> Aanpak: eerst begrijpen wat en waarom, dan de formules.

---

## V1 – Projector in $\mathbb{C}^m$

### (a) Definitie

**Wat is een projector?**

Stel je doet een schaduwprojectie: als je de schaduw al op de muur hebt, en je projecteert ze nogmaals, verandert er niets. Dat is precies het idee: een projector is een lineaire afbeelding die "stopt na één keer". Dubbel projecteren geeft hetzelfde als één keer.

**Formeel:**
Een matrix $P \in \mathbb{C}^{m \times m}$ is een **projector** als zij **idempotent** is:
$$P^2 = P.$$

**Waarop wordt geprojecteerd?** Op de **beeldruimte** $\mathcal{R}(P)$ — de ruimte van alle vectoren die je via $P$ kan bereiken.

**In welke richting?** Langs de **nulruimte** $\mathcal{N}(P)$ — de richting die volledig "weggeprojecteerd" wordt.

**De unieke ontbinding:** Elke vector $\mathbf{v}$ valt uiteen in een stuk op de beeldruimte en een stuk in de nulruimte:
$$\mathbf{v} = \underbrace{P\mathbf{v}}_{\in\,\mathcal{R}(P)} + \underbrace{(I-P)\mathbf{v}}_{\in\,\mathcal{N}(P)}.$$

Merk op: als $\mathbf{v}$ al in $\mathcal{R}(P)$ zit, dan $P\mathbf{v} = \mathbf{v}$ — projectie verandert niets aan een al-geprojecteerde vector.

> ⚠️ **Examentip:** Ken de definitie en de unieke ontbinding. Het idempotentie-bewijs gaat altijd via $P^2 = P$.

---

### (b) Expliciete uitdrukking voor $P_{S_1, S_2}$

**Wat is de situatie?**

We willen projecteren op een deelruimte $S_1$, maar niet loodrecht — schuin, in een richting zodat het residu loodrecht staat op een andere ruimte $S_2$.

Stel $A = [a_1 \cdots a_n]$ (basis van $S_1$) en $B = [b_1 \cdots b_n]$ (basis van $S_2$).

**Aanpak — drie stappen:**

We schrijven $P\mathbf{v} = A\mathbf{x}$ voor een onbekend $\mathbf{x}$ (het resultaat moet in $S_1$ zitten). Het residu $\mathbf{v} - A\mathbf{x}$ moet loodrecht staan op alle vectoren van $S_2$:
$$B^*(\mathbf{v} - A\mathbf{x}) = \mathbf{0}.$$

Los op voor $\mathbf{x}$:
$$B^*A\mathbf{x} = B^*\mathbf{v} \quad\Rightarrow\quad \mathbf{x} = (B^*A)^{-1}B^*\mathbf{v}.$$

De projector wordt:
$$\boxed{P_{S_1,S_2} = A(B^*A)^{-1}B^*.}$$

**Verificatie (idempotentie):**
$$P^2 = A(B^*A)^{-1}\underbrace{B^*A}_{=(B^*A)}(B^*A)^{-1}B^* = A(B^*A)^{-1}B^* = P. \checkmark$$

**Speciaal geval — orthogonale projector:** Als $B = A$ (projecteer loodrecht):
$$P_{S_1} = A(A^*A)^{-1}A^*.$$

> ⚠️ **Examentip:** Ken de algemene formule en haar afleiding in de drie stappen. Het speciaal geval met $B=A$ levert de orthogonale projector.

---

### (c) Bewijs: $P = P^*$ als en slechts als orthogonale projector

**Wat betekent "orthogonale projector"?**

Een projector heet orthogonaal als de beeldruimte loodrecht staat op de nulruimte: $\mathcal{R}(P) \perp \mathcal{N}(P)$. Dat is precies de situatie waarbij de projectie de "kortste" weg kiest.

**De wiskundige vertaling** van die loodrechtheid is simpelweg $P = P^*$.

**Bewijs $(\Leftarrow)$: als $P = P^*$, dan is $P$ orthogonaal**

Neem $\mathbf{x} = P\mathbf{u} \in \mathcal{R}(P)$ en $\mathbf{y} \in \mathcal{N}(P)$ willekeurig. Dan:
$$\langle \mathbf{x}, \mathbf{y} \rangle = (P\mathbf{u})^*\mathbf{y} = \mathbf{u}^*P^*\mathbf{y} \overset{P^*=P}{=} \mathbf{u}^*(P\mathbf{y}) = \mathbf{u}^*\mathbf{0} = 0.$$

Dus $\mathcal{R}(P) \perp \mathcal{N}(P)$. $\square$

**Bewijs $(\Rightarrow)$: als $P$ orthogonaal is, dan $P = P^*$**

Kies een orthonormale basis $\{\mathbf{q}_1,\ldots,\mathbf{q}_n\}$ van $\mathcal{R}(P)$ en $\{\mathbf{q}_{n+1},\ldots,\mathbf{q}_m\}$ van $\mathcal{N}(P)$. Omdat $\mathcal{R}(P) \perp \mathcal{N}(P)$, is de samengevoegde matrix $Q$ unitair. Dan:
$$Q^*PQ = \text{diag}(1,\ldots,1,0,\ldots,0).$$
Deze diagonaalmatrix is reëel symmetrisch, dus $(Q^*PQ)^* = Q^*PQ$, wat geeft $P = P^*$. $\square$

> ⚠️ **Examentip:** De kernstap is $\mathbf{u}^*P^*\mathbf{y} = \mathbf{u}^*P\mathbf{y} = 0$. Dit bewijs wordt bijna altijd gevraagd.

---

### (d) De uitmiddelingsmatrix $T$

**Wat doet $T$?**

$T$ vervangt elke component door het gemiddelde van alle componenten:
$$T\mathbf{x} = \frac{x_1 + x_2 + \cdots + x_m}{m}\begin{pmatrix}1\\\vdots\\1\end{pmatrix}.$$

In matrixvorm: $T = \tfrac{1}{m}\mathbf{1}\mathbf{1}^T$ met $\mathbf{1} = (1,\ldots,1)^T$.

**Is $T$ een projector?** Bereken $T^2$:
$$T^2 = \frac{1}{m^2}\mathbf{1}(\underbrace{\mathbf{1}^T\mathbf{1}}_{=m})\mathbf{1}^T = \frac{1}{m}\mathbf{1}\mathbf{1}^T = T. \checkmark$$

**Is $T$ een orthogonale projector?** Bereken $T^T$:
$$T^T = \left(\frac{1}{m}\mathbf{1}\mathbf{1}^T\right)^T = \frac{1}{m}\mathbf{1}\mathbf{1}^T = T. \checkmark$$

Omdat $T = T^T$ en $T^2 = T$: $T$ is een **orthogonale projector** — hij projecteert loodrecht op de ruimte van constante vectoren $\langle\mathbf{1}\rangle$.

**Intuïtie:** Het residu $(I-T)\mathbf{x}$ is de vector met het gemiddelde weggehaald — die heeft som nul en staat loodrecht op $\mathbf{1}$.

> ⚠️ **Examentip:** Twee korte berekeningen: $T^2 = T$ (projector) en $T = T^T$ (orthogonaal).

---

## V2 – Householder-transformaties en QR-factorisatie

### (a) Wat is een Householder-transformatie?

**Het idee in gewone taal:**

Stel je hebt een vector $\mathbf{x}$ en je wil hem met één operatie omzetten naar een veelvoud van $\mathbf{e}_1 = (1,0,\ldots,0)^T$ — alle componenten behalve de eerste op nul zetten. Dit doe je door $\mathbf{x}$ te **spiegelen** t.o.v. een slim gekozen hypervlak.

Een reflectie t.o.v. een hypervlak loodrecht op $\mathbf{v}$: de component langs $\mathbf{v}$ wordt omgekeerd, de rest blijft. Dit is equivalent met twee keer de projectie op $\langle\mathbf{v}\rangle$ aftrekken:
$$\boxed{F_m = I - 2\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}.}$$

**Verificatie — unitair:** $(F_m)^* F_m = (I - 2\tfrac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}})^2 = I. \checkmark$

**Keuze van $\mathbf{v}$:** Om $F_m\mathbf{x} = -\text{sign}(x_1)\|\mathbf{x}\|_2\,\mathbf{e}_1$ te bekomen, kies:
$$\mathbf{v} = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2\,\mathbf{e}_1.$$

**Geometrisch:** $F_m$ spiegelt $\mathbf{x}$ t.o.v. het hypervlak loodrecht op $\mathbf{v}$, zodat $\mathbf{x}$ precies op de $\mathbf{e}_1$-as terechtkomt.

> ⚠️ **Examentip:** Householder is een *reflector* — unitair en zelfinvers ($F_m^2 = I$). Zij brengt een volledige vector op de eerste standaardbasisvector.

---

### (b) Waarom rekening houden met het teken van $x_1$?

**Het probleem:**

Er zijn twee keuzen: $\mathbf{v} = \mathbf{x} \pm \|\mathbf{x}\|_2\mathbf{e}_1$. Zonder tekenconventie kies je willekeurig.

**De numerieke ramp bij de verkeerde keuze:**

Als $x_1 > 0$ en je kiest $\mathbf{v} = \mathbf{x} - \|\mathbf{x}\|_2\mathbf{e}_1$, dan is:
$$v_1 = x_1 - \|\mathbf{x}\|_2 \approx 0 \quad \text{als } \mathbf{x} \approx \|\mathbf{x}\|_2\,\mathbf{e}_1.$$

Dit is **gevaarlijke aftrekking** van bijna gelijke getallen — je verliest significante cijfers.

**De correcte keuze:**

Kies het teken zodat je altijd *optelt*:
$$\mathbf{v} = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2\,\mathbf{e}_1 \quad\Rightarrow\quad v_1 = x_1 + |x_1| \geq |x_1|.$$

Dan is $v_1$ nooit klein in vergelijking met $\|\mathbf{x}\|_2$.

> ⚠️ **Examentip:** De reden is altijd "gevaarlijke aftrekking vermijden" — niet meer en niet minder.

---

### (c) Algoritme voor QR-factorisatie: stap voor stap

**Wat is het doel?**

Breng een matrix $A \in \mathbb{C}^{m \times n}$ in bovendriehoekse vorm door opeenvolgend Householder-reflecties toe te passen. Elke stap brengt één kolom in orde.

#### Berekening van $R$

**Stap $j$:** beschouw de subkolom $\mathbf{x} = A_{j:m,j}$ (de elementen onder de diagonaal in kolom $j$).

1. Bereken $\mathbf{v}_j = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2\,\mathbf{e}_1$.
2. Pas de Householder-reflector toe op de relevante submatrix (efficiënt, zonder de $m\times m$ matrix expliciet te vormen):
$$A_{j:m,j:n} \leftarrow A_{j:m,j:n} - 2\frac{\mathbf{v}_j(\mathbf{v}_j^* A_{j:m,j:n})}{\mathbf{v}_j^*\mathbf{v}_j}.$$

Na $n$ stappen is $A$ bovendriehoeks: dat is $R$.

**Pseudocode:**
```
for j = 1 to n:
    x = A[j:m, j]
    v = x + sign(x[1]) * ||x||₂ * e₁
    A[j:m, j:n] -= 2*(v*v'*A[j:m,j:n]) / (v'*v)
```

#### Berekening van $Q$

Omdat $Q_{n,m}\cdots Q_{1,m} A = R$ en elke Householder-transformatie zelfinvers is:
$$A = \underbrace{Q_{1,m}\cdots Q_{n,m}}_{=Q}\, R.$$

In de praktijk bouwt men $Q$ op door de opgeslagen vectoren $\mathbf{v}_j$ achteraf op $I_m$ toe te passen (in omgekeerde volgorde).

**Rekenkost:**
- Enkel $R$: $\sim 2mn^2 - \tfrac{2}{3}n^3$ bewerkingen
- Inclusief $Q$: $\sim 4m^2n - \tfrac{2}{3}n^3$ bewerkingen

> ⚠️ **Examentip:** Ken de pseudocode en de efficiënte toepassing via $F_m\mathbf{q} = \mathbf{q} - 2\tfrac{\mathbf{v}^*\mathbf{q}}{\mathbf{v}^*\mathbf{v}}\mathbf{v}$ — dit vermijdt het expliciet opbouwen van de $m\times m$ matrix.

---

## V3 – Overgedetermineerd stelsel $Ax \approx b$ via QR

### (a) Kleinste-kwadratenoplossing via QR

**Wat is het probleem?**

Je hebt meer vergelijkingen dan onbekenden ($m > n$): geen exacte oplossing bestaat. Je zoekt de $\hat{x}$ waarvoor de residunorm $\|Ax - b\|_2$ minimaal is — de "best mogelijke" benadering.

**Geometrisch:** zoek het punt $A\hat{x}$ in de kolomruimte van $A$ dat het dichtst bij $b$ ligt.

**Hoe pak je het op met QR?**

Bereken $A = QR$ (volledige QR) en gebruik dat orthogonale transformaties de 2-norm bewaren:
$$\|Ax - b\|_2 = \|Q^T(Ax - b)\|_2 = \left\|\begin{pmatrix}\hat{R}x - b_1 \\ -b_2\end{pmatrix}\right\|_2^2 = \|\hat{R}x - b_1\|_2^2 + \|b_2\|_2^2,$$

waarbij $Q^Tb = \begin{pmatrix}b_1 \\ b_2\end{pmatrix}$ met $b_1 \in \mathbb{R}^n$.

Enkel de eerste term hangt af van $x$. Het minimum wordt bereikt door:
$$\hat{R}\hat{x} = b_1 \quad \text{(opgelost via achterwaartse substitutie).}$$

De **minimale restfout** is altijd $\|b_2\|_2$ — dit deel kan geen enkele $x$ elimineren.

**Algoritme:**
1. Bereken $A = \hat{Q}\hat{R}$ (onvolledige QR)
2. Bereken $b_1 = \hat{Q}^T b$
3. Los $\hat{R}\hat{x} = b_1$ op via achterwaartse substitutie

> ⚠️ **Examentip:** De sleutelstap is de splitsing in $b_1$ en $b_2$ na de orthogonale transformatie.

---

### (b) Onvolledige of volledige QR?

**De onvolledige (dunne) QR** $A = \hat{Q}\hat{R}$ volstaat.

**Waarom?** De oplossing $\hat{x}$ hangt alleen af van $b_1 = \hat{Q}^T b$ (de eerste $n$ componenten). De rest $b_2$ geeft enkel de minimale fout, maar je hebt die niet nodig om $\hat{x}$ te berekenen.

**Geometrisch:** $b_2 = (I - \hat{Q}\hat{Q}^T)b$ is de projectie van $b$ op het orthogonale complement van $\mathcal{R}(A)$. Geen keuze van $x$ kan dit deel elimineren — het is de intrinsieke, onvermijdelijke restfout.

**Praktisch voordeel:** Rekenkost $O(mn^2)$ i.p.v. $O(m^2n)$ voor de volledige factorisatie.

---

### (c) $A$ is $10^7 \times 10$ — welk algoritme?

**Antwoord: Gram-Schmidt** (gewijzigd Gram-Schmidt).

**Waarom niet Householder?**

Bij Householder heb je de volledige $Q \in \mathbb{R}^{10^7 \times 10^7}$ nodig — dat zijn $10^{14}$ getallen. Totaal onhaalbaar qua geheugen.

Bij Gram-Schmidt bouw je alleen $\hat{Q} \in \mathbb{R}^{10^7 \times 10}$ op — slechts $10^8$ getallen. Dat is haalbaar.

| Eigenschap | Gram-Schmidt | Householder | Givens |
|---|---|---|---|
| Rekenkost $\hat{Q}$ + $\hat{R}$ | $O(2mn^2)$ | $O(4m^2n)$ | $O(6m^2n)$ |
| Geheugen | $O(mn)$ | $O(m^2)$ voor volle $Q$ | $O(m^2)$ |

Voor $m = 10^7$, $n = 10$: Householder kost een factor $\sim 10^6$ meer dan Gram-Schmidt.

> ⚠️ **Examentip:** Bij $m \gg n$: altijd Gram-Schmidt. Reden: geheugen en kost van volledige $Q$ bij Householder.

---

### (d) Waarom beter dan het normaalstelsel oplossen?

**Het probleem met het normaalstelsel** $A^TA\hat{x} = A^Tb$:

Het kwadreren van $A$ kwadreert ook het conditiegetal:
$$\kappa(A^T A) = \kappa(A)^2.$$

**Gevolg:** een fout op $b$ of $A$ wordt $\kappa(A)^2$ keer vergroot bij het normaalstelsel, maar slechts $\kappa(A)$ keer bij de QR-aanpak.

| Methode | Relatieve fout |
|---|---|
| Normaalstelsel | $O(\kappa(A)^2 \cdot \epsilon_{\text{mach}})$ |
| QR-aanpak | $O(\kappa(A) \cdot \epsilon_{\text{mach}})$ |

Voor $\kappa(A) = 1000$: normaalstelsel geeft fouten $10^6 \times$ groter dan QR. Het verschil tussen zinvol en volledig foutief.

> ⚠️ **Examentip:** Onthoud $\kappa(A^TA) = \kappa(A)^2$ en de vergelijking van de fouten.

---

## V4 – De maximumnorm op $C[a,b]$

### (a) Bewijs dat $\|f\|_\infty = \max_{x\in[a,b]}|f(x)|$ een norm is

**Wat betekent "norm"?** Een norm meet de "grootte" van een vector (hier: functie). Ze moet aan vier eigenschappen voldoen.

**Eigenschap 1 — Niet-negativiteit:** $|f(x)| \geq 0$, dus ook het maximum is $\geq 0$. ✓

**Eigenschap 2 — Nulvector:** $\|f\|_\infty = 0$ als en slechts als $|f(x)| = 0$ voor alle $x$, dus $f \equiv 0$. (Cruciaal: continue functie op gesloten begrensd interval — het maximum is bereikbaar.) ✓

**Eigenschap 3 — Homogeniteit:** $\|\alpha f\|_\infty = \max_x|\alpha f(x)| = |\alpha|\max_x|f(x)| = |\alpha|\,\|f\|_\infty$. ✓

**Eigenschap 4 — Driehoeksongelijkheid:** Voor alle $x$:
$$|f(x) + g(x)| \leq |f(x)| + |g(x)| \leq \|f\|_\infty + \|g\|_\infty.$$
Neem het maximum over $x$: $\|f + g\|_\infty \leq \|f\|_\infty + \|g\|_\infty$. ✓

$I(f) = \|f\|_\infty$ is dus een norm op $C[a,b]$. $\square$

> ⚠️ **Examentip:** De driehoeksongelijkheid (eigenschap 4) is de meest gevraagde stap — schrijf ze expliciet uit.

---

### (b) Is het een strikte norm?

**Wat betekent "strikt"?** Een norm is strikt als gelijkheid in de driehoeksongelijkheid alleen geldt als één functie een niet-negatief veelvoud is van de andere.

**Antwoord: Nee**, de maximumnorm is **niet** strikt.

**Tegenvoorbeeld:** $[a,b] = [0,1]$, $f(x) = 1$, $g(x) = x$. Dan:
$$\|f\|_\infty = 1, \quad \|g\|_\infty = 1, \quad \|f+g\|_\infty = \max_{[0,1]}(1+x) = 2 = \|f\|_\infty + \|g\|_\infty,$$
maar $f$ is geen veelvoud van $g$.

**Gevolg:** De beste benadering in de maximumnorm is niet altijd uniek.

**Intuïtie:** De eenheidsbol van de maximumnorm heeft "rechte stukken" — twee functies die tegelijk hun maximum bereiken zijn even "groot". Bij de $L^2$-norm (bol = echt rond) kan dat niet.

---

### (c) Voldoet deze norm aan de parallellogramgelijkheid?

**Wat is de parallellogramgelijkheid?** $\|f+g\|^2 + \|f-g\|^2 = 2(\|f\|^2 + \|g\|^2)$.

Als een norm hieraan voldoet, komt ze voort uit een scalair product.

**Antwoord: Nee.**

**Tegenvoorbeeld:** $f(x) = 1$, $g(x) = x$ op $[0,1]$:
$$\|f+g\|_\infty^2 + \|f-g\|_\infty^2 = 4 + 1 = 5, \qquad 2(\|f\|_\infty^2 + \|g\|_\infty^2) = 2(1+1) = 4.$$

$5 \neq 4$ — de gelijkheid geldt niet.

**Gevolg:** De maximumnorm komt **niet** voort uit een scalair product. Je kan geen $(\cdot,\cdot)$ vinden zodat $\sqrt{(f,f)} = \|f\|_\infty$.

**Vergelijk:** De $L^2$-norm $\|f\|_2 = \sqrt{\int_a^b f^2\,dx}$ voldoet wél — zij is de geïnduceerde norm van $(f,g) = \int_a^b fg\,dx$.

> ⚠️ **Examentip:** De parallellogramgelijkheid is het criterium (stelling van Jordan-von Neumann). Hetzelfde tegenvoorbeeld $f=1$, $g=x$ werkt altijd voor de maximumnorm.

---

## V5 – Gram-Schmidt in unitaire ruimte

### (a) Gram-Schmidt algoritme + verband met orthogonale projecties

**Wat is het doel?**

Je hebt een basis $\{a_1,\ldots,a_n\}$ en wil er een **orthonormale** basis van maken: vectoren die loodrecht op elkaar staan én lengte 1 hebben. Gram-Schmidt doet dit één vector per keer.

**Het idee per stap:**

Om $q_j$ te bepalen, begin je met $a_j$ en trek je er alle componenten van af die al in de richting van de eerder bepaalde $q_1,\ldots,q_{j-1}$ liggen. Wat overblijft staat automatisch loodrecht op alle vorige vectoren.

**Algoritme:**
```
voor j = 1 tot n:
    v_j = a_j
    voor i = 1 tot j-1:
        r_{ij} = (q_i, a_j)         ← scalair product
        v_j -= r_{ij} · q_i        ← verwijder component langs q_i
    r_{jj} = ||v_j||
    q_j = v_j / r_{jj}             ← normaliseer
```

**Verband met projecties:**

Stap $j$ berekent precies:
$$v_j = a_j - \underbrace{\sum_{i=1}^{j-1}(q_i, a_j)\,q_i}_{=P_{\langle q_1,\ldots,q_{j-1}\rangle}\,a_j}.$$

We trekken de **loodrechte projectie** van $a_j$ op de reeds gebouwde ruimte af. Na normalisatie geeft dit $q_j$.

**Relatie met QR:** De coëfficiënten $r_{ij}$ vormen een bovendriehoekse matrix $\hat{R}$, en we krijgen de abstracte analogon van $A = \hat{Q}\hat{R}$:
$$a_j = \sum_{i=1}^j r_{ij}\,q_i.$$

> ⚠️ **Examentip:** De geometrische interpretatie "projecteer en trek af" is essentieel. $v_j$ staat loodrecht op de vorige ruimte precies omdat we de projectie aftrekken.

---

### (b) Aantal scalaire producten

**Per stap $j$:** $j-1$ scalaire producten voor de coëfficiënten $r_{ij}$, plus 1 voor de norm.

**Totaal:**
$$\sum_{j=1}^n j = \frac{n(n+1)}{2} = O(n^2) \text{ scalaire producten.}$$

**Vergelijk:** ook het normaalstelsel kost $O(n^2)$ scalaire producten, maar is numeriek gevaarlijker ($\kappa(G) = \kappa(A)^2$).

> ⚠️ **Examentip:** De drietermsrecursiebetrekking voor orthogonale veeltermen (V6) reduceert dit tot slechts $O(n)$ scalaire producten — enorme besparing voor die specifieke klasse.

---

### (c) Herorthogonalisatie

**Wat is het?**

Na Gram-Schmidt pas je het algoritme nogmaals toe op de uitvoer $\{q_1,\ldots,q_n\}$.

**Wanneer nodig?**

Gram-Schmidt is wiskundig correct maar **numeriek niet stabiel** bij slecht geconditioneerde bases (vectoren liggen bijna in dezelfde richting). Afrondingsfouten tasten de orthogonaliteit aan:
$$(q_i, q_j) \approx \epsilon_{\text{mach}} \cdot \kappa(A).$$

Na herorthogonalisatie dalen de fouten naar $\epsilon_{\text{mach}}^2 \cdot \kappa(A)^2$ — in de meeste praktische gevallen verwaarloosbaar.

**Gewijzigd Gram-Schmidt (MGS) — numeriek stabielere variant:**

Het verschil: bij klassiek GS worden alle componenten berekend t.o.v. de **originele** $a_j$. Bij MGS wordt na elke aftrekking de **bijgewerkte** $v_j$ gebruikt voor de volgende stap:

```
voor j = 1 tot n:
    v_j = a_j
    voor i = 1 tot j-1:
        r_{ij} = (q_i, v_j)         ← inproduct met bijgewerkte v_j!
        v_j -= r_{ij} · q_i
    r_{jj} = ||v_j||
    q_j = v_j / r_{jj}
```

Wiskundig identiek aan klassiek GS, maar numeriek significant robuuster: fouten worden na elke deelstap opgeschoond.

> ⚠️ **Examentip:** Sleutelwoord: **Modified Gram-Schmidt (MGS)**. Het verschil: "stapsgewijs verwijderen t.o.v. bijgewerkte vector" i.p.v. "alles t.o.v. originele vector".

---

### (d) Orthogonale veeltermen opstellen met Gram-Schmidt: efficiënt?

**Antwoord: Nee**, niet de meest efficiënte methode.

Gram-Schmidt kost $O(n^2)$ scalaire producten — elk een integraal. Dat is duur.

**Betere aanpak: drietermsrecursiebetrekking (V6).**

Voor orthogonale veeltermen geldt dat elke nieuwe veelterm $\phi_k$ alleen afhangt van $\phi_{k-1}$ en $\phi_{k-2}$. Dat vraagt slechts **2 scalaire producten per stap**: totaal $O(n)$ i.p.v. $O(n^2)$.

**Waarom kan dit?** Veeltermen hebben een speciale structuur: het scalair product van $x\phi_{k-1}$ met $\phi_j$ is nul voor alle $j$ behalve $j = k-2, k-1, k$. Die structuur benut Gram-Schmidt niet, maar de drietermsrecursie wel.

---

## V6 – Orthogonale rij veeltermen $\{\phi_0, \phi_1, \phi_2, \ldots\}$

### (a) Nulpunten van $\phi_n$

**Stelling:** $\phi_n(x)$ heeft precies **$n$ enkelvoudige reële nulpunten**, allemaal in het **open interval** $(a,b)$.

**Vier argumenten:**

1. **Reëel:** alle coëfficiënten van $\phi_n$ zijn reëel.
2. **Precies $n$:** $\phi_n$ is van graad $n$ (hoogstens $n$ nulpunten). De orthogonaliteit t.o.v. constanten $(\phi_n, 1) = 0$ dwingt minstens één nulpunt in $(a,b)$. Een volledig argument toont dat er minstens $n$ tekenwisselingen zijn.
3. **Enkelvoudig:** meervoudige nulpunten leiden tot een contradictie met de orthogonaliteitseis.
4. **Open interval:** de randpunten zijn nooit nulpunten (volgt uit normalisatie en positief gewicht).

> ⚠️ **Examentip:** Weet dit te formuleren. Het volledige bewijs hoeft niet, maar ken het tegenspraakargument via het teken van $\phi_n \cdot q$.

---

### (b) Nulpunten zijn eigenwaarden van een tridiagonale matrix

**Wat is het idee?**

De drietermsrecursiebetrekking koppelt opeenvolgende veeltermen. Als je die recursie als een matrixvergelijking herschrijft, zie je dat de nulpunten van $\phi_n$ precies de situaties zijn waarbij die matrixvergelijking een eigenwaardevergelijking wordt.

**De drietermsrecursie:**
$$\nu_{k-1}\phi_{k-2}(x) + \alpha_k\phi_{k-1}(x) + \mu_k\phi_k(x) = x\phi_{k-1}(x), \quad k = 1,\ldots,n.$$

**Stap 1:** Definieer $\Phi(x) = [\phi_0(x), \phi_1(x), \ldots, \phi_{n-1}(x)]^T$ en schrijf de recursie voor $k=1,\ldots,n$ als:
$$A\Phi(x) = x\Phi(x) - \mu_n\phi_n(x)\,\mathbf{e}_n,$$
met $A$ de tridiagonale matrix met $\alpha_k$ op de diagonaal en $\mu_k, \nu_k$ op de zij-diagonalen.

**Stap 2:** Evalueer in een nulpunt $x_j$ van $\phi_n$ (dus $\phi_n(x_j) = 0$). De stoorterm valt weg:
$$A\Phi(x_j) = x_j\Phi(x_j).$$

**Conclusie:** $x_j$ is een **eigenwaarde** van $A$ met eigenvector $\Phi(x_j)$. Omdat er $n$ nulpunten zijn en $A$ een $n\times n$ matrix is, zijn de nulpunten en eigenwaarden precies dezelfde $n$ getallen. $\square$

**Praktisch belang:** De nulpunten worden berekend door het eigenwaardenprobleem van $A$ op te lossen via het QR-algoritme — efficiënt en nauwkeurig in $O(n^2)$ bewerkingen.

> ⚠️ **Examentip:** Sleutelstap: evalueer in nulpunt $\phi_n(x_j) = 0$ → stoorterm valt weg → eigenwaardevergelijking.

---

### (c) Interpolatie in nulpunten van Chebyshev: goede aanpak?

**Antwoord: Ja, uitstekend.**

**Waarom zijn Chebyshev-nulpunten goed?**

De Chebyshev-veelterm $T_{n+1}$ oscilleert gelijkmatig tussen $-1$ en $+1$ over $[-1,1]$ — precies $n+2$ extrema met gelijke amplitude. Die eigenschap geeft drie voordelen:

1. **Bijna gelijk aan kleinste-kwadratenbenadering:** de interpolatiefout in Chebyshev-punten is bijna identiek aan de optimale $L^2$-benadering.
2. **Bijna minimax:** de fout oscilleert gelijkmatig — bijna optimaal in maximumnorm.
3. **Vermijdt Runge-fenomeen:** bij equidistante punten explodeert de fout bij hoge graden. Chebyshev-punten concentreren zich aan de randen en vermijden dit.

**De Chebyshev-punten:**
$$x_k = \cos\!\left(\frac{(2k-1)\pi}{2(n+1)}\right), \quad k = 1,\ldots,n+1.$$

**Concreet voorbeeld ($e^x$ op $[-1,1]$, graad 4):**

| Methode | Max. fout |
|---|---|
| Taylor | $\approx 10^{-2}$ |
| Equidistante interpolatie | $\approx 10^{-3}$ |
| Chebyshev-interpolatie | $\approx 3 \times 10^{-4}$ |
| Kleinste kwadraten | $\approx 3 \times 10^{-4}$ |

> ⚠️ **Examentip:** Drie punten: (1) bijna KK-benadering, (2) bijna minimax, (3) vermijdt Runge. Alles volgt uit de gelijkmatige oscillatie van $T_{n+1}$.

---

## V7 – Grafen: directe verbindingen en meest centrale knoop

### Directe verbindingen via de verbindingsmatrix

**Wat is een verbindingsmatrix?**

Voor een netwerk van $N$ knopen maak je een matrix $A$ met $A_{ij} = 1$ als knoop $i$ en $j$ verbonden zijn, en $0$ anders. De matrix is symmetrisch (ongerichte graaf) en heeft nullen op de diagonaal (geen zelfverbindingen).

**Aantal directe buren van knoop $i$:**

Tel gewoon de enen in rij $i$:
$$\text{graad}(i) = \sum_j A_{ij} = (A^2)_{ii}.$$

(Want $(A^2)_{ii} = \sum_j A_{ij}A_{ji} = \sum_j A_{ij}^2 = \sum_j A_{ij}$ omdat $A_{ij} \in \{0,1\}$.)

**Lussen van lengte $n$:**

$(A^n)_{ij}$ telt het aantal paden van lengte $n$ van knoop $i$ naar knoop $j$. Voor $i=j$: het aantal lussen van lengte $n$ die starten en eindigen in $i$.

---

### Definitie meest centrale knoop

**Het probleem met enkel directe verbindingen:**

Een knoop kan veel buren hebben die onderling slecht verbonden zijn. Een andere knoop met minder directe buren maar goed verbonden met invloedrijke knopen is in de praktijk centraler.

**De oplossing: combineer lussen van alle lengtes** met afnemend gewicht:
$$\text{centraliteit}(i) = \sum_{n=0}^\infty \frac{1}{n!}(A^n)_{ii} = (e^A)_{ii}.$$

De factor $1/n!$ zorgt dat lange lussen minder zwaar wegen (en dat de reeks convergeert).

**Bijdragen per lengte:**
- $n=0$: triviale bijdrage 1
- $n=1$: 0 (geen zelfverbindingen)
- $n=2$: aantal directe buren
- $n=3$: driehoeken door $i$
- $n=4$: vierhoeken door $i$, enzovoort

**De meest centrale knoop** heeft de grootste $(e^A)_{ii}$.

---

### Berekening van de centraliteit

**Hoe bereken je $e^A$?**

Gebruik de eigenwaardenontbinding $A = X\Lambda X^T$ (reëel en symmetrisch → $X$ orthogonaal):
$$e^A = X\begin{pmatrix}e^{\lambda_1}&&\\&\ddots&\\&&e^{\lambda_N}\end{pmatrix}X^T.$$

Dit is het **spectraal mapping theorema**: $f(A) = Xf(\Lambda)X^{-1}$ voor elke functie $f$, waarbij $f$ componentgewijs op de eigenwaarden wordt toegepast.

Lees dan $(e^A)_{ii}$ af. Het QR-algoritme berekent de volledige eigenwaardenontbinding in $O(N^3)$ bewerkingen.

> ⚠️ **Examentip:** Ken de definitie $(e^A)_{ii}$ en de berekening via eigenwaardenontbinding. Formule: $e^A = Xe^\Lambda X^T$.

---

## V8 – Het algoritme van Arnoldi

### (a) Voor welk soort problemen?

**Het probleem:**

Je hebt een grote matrix $A$ ($m = 10^6$ of meer) en je wil enkele **extremale eigenwaarden** berekenen. Het volledige QR-algoritme kost $O(m^3)$ bewerkingen en vereist de volledige matrix in geheugen — beide onhaalbaar.

**Arnoldi's oplossing:**

Gebruik alleen **matrix-vector producten** $A\mathbf{v}$ (voor ijle matrices: $O(\text{nnz})$ per product). Bouw stapsgewijs een kleine ruimte op waarop je het eigenwaardenprobleem projecteert. Na $k \ll m$ stappen heb je $k$ goede benaderende eigenwaarden.

**Concreet voorbeeld: PageRank.**

De linkmatrix van het web heeft $N \sim 10^{10}$ pagina's maar is ijl ($\sim 10$ links per pagina). We zoeken de dominante eigenvector. Arnoldi (of zijn variant: methode van de machten) lost dit op met enkel matrix-vector producten.

> ⚠️ **Examentip:** Twee sleutelwoorden: "grote ijle matrix" en "slechts een paar extremale eigenwaarden". Het voordeel: enkel matrix-vector producten nodig.

---

### (b) Orthonormale basis van de Krylov-ruimte

**Wat is de Krylov-ruimte?**

Herhaaldelijk toepassen van $A$ op een startvector $\mathbf{b}$ geeft vectoren $\mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots$ Na $k$ stappen beslaan die een $k$-dimensionale ruimte:
$$\mathcal{K}_k(A,\mathbf{b}) = \text{span}\{\mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}.$$

**Het probleem:** de Krylov-vectoren worden snel evenredig met de dominante eigenvector — gevaarlijke numerieke afhankelijkheid.

**Arnoldi's aanpak:** pas Gram-Schmidt toe op de opeenvolgend gegenereerde vectoren. Het algoritme bouwt zo een orthonormale basis $\{q_1,\ldots,q_k\}$ van $\mathcal{K}_k$ op.

**Arnoldi Algoritme:**
```
q₁ = b / ||b||₂
voor j = 1, 2, 3, ...:
    v = A qⱼ                      ← matrix-vector product
    voor i = 1 tot j:
        h_{ij} = qᵢ* v            ← scalair product (Gram-Schmidt)
        v -= h_{ij} qᵢ
    h_{j+1,j} = ||v||₂
    q_{j+1} = v / h_{j+1,j}       ← normaliseer
```

**Arnoldi-recursievergelijking** (na $k$ stappen):
$$\boxed{AQ_k = Q_kH_k + h_{k+1,k}\,\mathbf{q}_{k+1}\mathbf{e}_k^T,}$$
met:
- $Q_k = [\mathbf{q}_1\cdots\mathbf{q}_k]$: orthonormale kolommen
- $H_k = Q_k^*AQ_k$: **Hessenberg-matrix** — de projectie van $A$ op de Krylov-ruimte

> ⚠️ **Examentip:** Ken de recursievergelijking en het feit dat $H_k = Q_k^*AQ_k$.

---

### (c) Eigenwaarden benaderen: Ritz-waarden

**Wat zijn Ritz-waarden?**

In stap $k$ zoek je een benadering van een eigenpaar van $A$ in de Krylov-ruimte $\mathcal{K}_k$. Eis dat het residu loodrecht staat op $\mathcal{K}_k$ (Galerkin-eis):
$$Q_k^*(A - \theta I)Q_k\mathbf{y} = \mathbf{0} \quad\Rightarrow\quad H_k\mathbf{y} = \theta\mathbf{y}.$$

De **Ritz-waarden** $\theta_1^{(k)},\ldots,\theta_k^{(k)}$ zijn de eigenwaarden van de kleine $k\times k$ matrix $H_k$.

**Het voordeel:**
- Eigenwaarden van $A$: $m\times m$ probleem (onhaalbaar)
- Ritz-waarden: $k\times k$ probleem ($k \ll m$) — oplosbaar met gewoon QR-algoritme

**Vroegtijdig stoppen:** als $h_{k+1,k} = 0$ zijn de Ritz-waarden **exact** eigenwaarden van $A$.

> ⚠️ **Examentip:** Ritz-waarden = eigenwaarden van $H_k = Q_k^*AQ_k$. Dit is de sleutelformule.

---

### (d) Variante met gewijzigd Gram-Schmidt?

**Antwoord: Ja, zinvol en standaard gebruikt.**

MGS is numeriek stabieler dan klassiek GS zonder extra kost ($O(km)$ per stap in beide gevallen). Bij ernstig slechte conditie kan ook MGS zijn orthogonaliteit verliezen — dan is herorthogonalisatie nodig.

**Concreet in Arnoldi (MGS-variant):**
```
v = A qⱼ
voor i = 1 tot j:
    h_{ij} = qᵢ* v          ← inproduct met bijgewerkte v (MGS)
    v -= h_{ij} qᵢ           ← onmiddellijk verwijderen
```

**Standaardimplementatie:** Arnoldi met MGS (+ herorthogonalisatie indien nodig).

> ⚠️ **Examentip:** "Ja, zinvol" — MGS is stabieler, zelfde kost.

---

## V9 – Convergentiegedrag van Arnoldi

### Eigenschap 1 — Verschuivingsinvariantie

De Arnoldi-recursie voor $(A - \sigma I, \mathbf{b})$ geeft dezelfde $Q_k$ en $H_k - \sigma I$ als Hessenberg-projectie. Daardoor convergeert Arnoldi zowel naar $\lambda_{\min}$ als $\lambda_{\max}$ — in tegenstelling tot de methode van de machten (die alleen naar de grootste convergeert).

De stoorterm $h_{k+1,k}$ meet hoe goed $\mathcal{K}_k$ de gezochte eigenruimte benadert: **klein $h_{k+1,k}$ $\iff$ Ritz-waarden zijn geconvergeerd**.

---

### Eigenschap 2 — Veeltermbenadering en lokalisatie

**Welke eigenwaarden convergeren eerst?**

De optimale veelterm $\hat{p}$ (karakteristiek van $H_k$) minimaliseert $\|p(A)\mathbf{b}\|_2$. Via de eigenwaardenontbinding $A = X\Lambda X^{-1}$ plaatst $\hat{p}$ zijn nulpunten zo dat $|p(\lambda_i)|$ klein is voor de eigenwaarden die sterk bijdragen aan $\mathbf{b}$.

**Gevolg:** De **geïsoleerde eigenwaarden aan de rand van het spectrum** convergeren het snelst — een nulpunt van $p$ kan dicht bij $\lambda_j$ gezet worden zonder andere eigenwaarden te verstoren. Geclusterde eigenwaarden in het midden vereisen veel iteraties.

---

## V10 – Twee methodes voor lage-rangbenadering

**Het probleem:** Gegeven $A \in \mathbb{C}^{m \times n}$, zoek de beste benadering van rang $\leq k$.

De verzameling lage-rangmatrices is **geen deelruimte** (niet gesloten onder optelling) — orthogonale projectietheorie werkt hier niet.

### Methode 1 — QR met kolomverwisselingen

**Idee:** bij elke stap pivoteer je de kolom met de grootste norm naar voor:
$$AP = QR, \qquad |r_{11}| \geq |r_{22}| \geq \cdots$$

De rang-$k$ benadering: $A_k^{qr} = Q_k R_k P^T$.

- **Rekenkost:** $O(mnk)$
- **Kwaliteit:** goed maar **niet gegarandeerd optimaal**

---

### Methode 2 — Getrunceerde SVD

**Idee:** de Singulierwaarden-ontbinding $A = U\Sigma V^*$ sorteert de "richtingen" van grootste naar kleinste bijdrage. Houd gewoon de $k$ grootste:
$$A_k^{svo} = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^* = U_k\Sigma_k V_k^*.$$

**Fout:**
$$\|A - A_k^{svo}\|_2 = \sigma_{k+1}, \qquad \|A - A_k^{svo}\|_F = \sqrt{\sigma_{k+1}^2 + \cdots + \sigma_r^2}.$$

- **Rekenkost:** $O(mn^2)$ (duurder)
- **Kwaliteit:** **optimaal** (Stelling Eckart–Young–Mirsky)

---

### Stelling — Beste rang-$k$ benadering (Eckart–Young–Mirsky)

**Wat zegt het?** Er bestaat geen rang-$k$ matrix die $A$ beter benadert dan $A_k^{svo}$.

**Bewijs:** Voor elke rang-$k$ matrix $B = XY^*$ zoek je een $\mathbf{w}$ met $\|\mathbf{w}\|_2 = 1$ in zowel $\text{span}\{\mathbf{v}_1,\ldots,\mathbf{v}_{k+1}\}$ als $\mathcal{N}(Y^*)$. Zo'n $\mathbf{w}$ bestaat (dimensies overlappen: $(k+1)+(n-k) = n+1 > n$).

Schrijf $\mathbf{w} = \sum_{i=1}^{k+1} c_i \mathbf{v}_i$ met $\sum c_i^2 = 1$. Dan $B\mathbf{w} = 0$, dus:
$$\|B - A\|_2 \geq \|A\mathbf{w}\|_2 = \sqrt{\sum c_i^2 \sigma_i^2} \geq \sigma_{k+1}. \quad\square$$

---

### Vergelijking

| | QR met kolomverwisselingen | Getrunceerde SVD |
|---|---|---|
| **Kwaliteit** | Suboptimaal | Optimaal |
| **Rekenkost** | $O(mnk)$ | $O(mn^2)$ |
| **Fout gekend?** | Niet precies | $\sigma_{k+1}$ exact |
| **Wanneer** | Snelheid primeert | Beste kwaliteit |

---

### Geval: gigantische ijle matrix ($10^7 \times 10^7$)

**Beide methodes zijn onhaalbaar:** volledige SVD kost $O(m^3) = O(10^{21})$ bewerkingen en de matrix opslaan vereist $10^{14}$ getallen.

**Oplossing: Krylov-methoden** — enkel matrix-vector producten $A\mathbf{v}$ nodig.

Voor **symmetrische** matrices: **Lanczos-algoritme** (= Arnoldi voor symmetrische matrices). De Hessenbergmatrix $H_k$ vereenvoudigt tot een **tridiagonale** matrix $T_k$. De singuliere waarden van $A$ volgen via $\sigma_i(A) = \sqrt{\lambda_i(A^TA)}$.

| Methode | Rekenkost | Geheugen | Haalbaar bij $10^7\times 10^7$? |
|---|---|---|---|
| Volledige SVD | $O(mn^2)$ | $O(mn)$ | Nee |
| QR met kolomverwisseling | $O(mnk)$ | $O(mn)$ | Nee |
| Lanczos/Arnoldi | $O(k\cdot\text{nnz})$ | $O(kn)$ | **Ja** |

> ⚠️ **Examentip:** Bij matrices van orde $10^7$ of meer: het antwoord is altijd **Krylov-methode (Lanczos voor SVD, Arnoldi voor niet-symmetrisch)**. Sleutelzin: "enkel matrix-vector producten, matrix hoeft niet in geheugen."

---

## Snel-referentie: Sleutelformules

| Formule | Betekenis |
|---|---|
| $P^2 = P$ | Definitie projector |
| $P_{S_1,S_2} = A(B^*A)^{-1}B^*$ | Scheve projector op $S_1$ langs $S_2^\perp$ |
| $P_D = A(A^*A)^{-1}A^* = QQ^*$ | Orthogonale projector op $\mathcal{D}$ |
| $P = P^* \iff P$ orthogonaal | Karakterisering orthogonale projector |
| $F_m = I - 2\tfrac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}$ | Householder-reflector |
| $\hat{R}\hat{x} = b_1 = \hat{Q}^T b$ | KK-oplossing via QR |
| $\kappa(A^TA) = \kappa(A)^2$ | Reden voor instabiliteit normaalstelsel |
| $AQ_k = Q_kH_k + h_{k+1,k}\mathbf{q}_{k+1}\mathbf{e}_k^T$ | Arnoldi-recursie |
| $H_k = Q_k^*AQ_k$ | Ritz-matrix = geprojecteerde $A$ |
| $A_k^{svo} = U_k\Sigma_kV_k^*$ | Getrunceerde SVD |
| $\|A - A_k^{svo}\|_2 = \sigma_{k+1}$ | Fout van getrunceerde SVD |
| $(e^A)_{ii}$ | Centraliteit van knoop $i$ |
