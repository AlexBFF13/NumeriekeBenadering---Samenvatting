# Modelexamen Wim Michiels — Uitgewerkte antwoorden

> Alle vragen komen uit `modelvragen_deel_WM(1).pdf`.  
> Bronnen: cursustekst + samenvattingen H2, H3, H4, H8, H9, H12.  
> Formularium is beschikbaar op het examen — bewijzen en afleidsels worden hieronder volledig uitgewerkt.

---

## V1 – Projector in $\mathbb{C}^m$

**Vraag:** Bespreek het begrip projector in de ruimte $\mathbb{C}^m$.
(a) Geef definitie. Waarop wordt geprojecteerd? In welke richting?
(b) Leid een expliciete uitdrukking af voor $P_{S_1,S_2}$ die projecteert op $S_1 = \langle a_1,\ldots,a_n\rangle$ in de richting van $S_2 = \langle b_1,\ldots,b_n\rangle^\perp$.
(c) Bewijs dat $P$ een orthogonale projector is indien $P = P^*$.
(d) Beschouw de uitmiddelingsmatrix $T$. Is $T$ een projector? Een orthogonale projector?

**Kernidee:** Een projector is een lineaire afbeelding die "stopt na één keer": dubbel projecteren geeft hetzelfde als één keer. Orthogonale projectoren projecteren loodrecht, en dat is precies wat we nodig hebben voor beste benaderingen.

---

### (a) Definitie

### Definitie — Projector
Een matrix $P \in \mathbb{C}^{m \times m}$ is een **projector** als zij **idempotent** is:
$$P^2 = P.$$

**Intuïtie:** Als je een vector eenmaal projecteert ($P\mathbf{v}$) en daarna opnieuw ($P(P\mathbf{v}) = P^2\mathbf{v}$), kom je op hetzelfde punt uit — projecteren is een "eenmalige" operatie.

**Waarop wordt geprojecteerd?** Op de **beeldruimte** $\mathcal{R}(P)$.

**In welke richting?** Langs de **nulruimte** $\mathcal{N}(P)$.

**Structuur:** Elke vector $\mathbf{v} \in \mathbb{C}^m$ ontbindt uniek als:
$$\mathbf{v} = \underbrace{P\mathbf{v}}_{\in\,\mathcal{R}(P)} + \underbrace{(I-P)\mathbf{v}}_{\in\,\mathcal{N}(P)}.$$

**Bijkomende eigenschappen:**
- $\mathcal{R}(P) \cap \mathcal{N}(P) = \{\mathbf{0}\}$
- $\dim(\mathcal{R}(P)) + \dim(\mathcal{N}(P)) = m$
- Als $\mathbf{v} \in \mathcal{R}(P)$: $P\mathbf{v} = \mathbf{v}$ (projectie op eigen subspace = identiteit)
- $\tilde{P} = I - P$ is ook een projector (de complementaire projector)

> ⚠️ **Examentip:** Weet de definitie en de unieke ontbinding vanbuiten. Het idempotentie-bewijs gaat altijd via $P^2 = P$.

---

### (b) Expliciete uitdrukking voor $P_{S_1, S_2}$

**Opstelling:** We projecteren op $S_1 = \langle a_1, \ldots, a_n \rangle$ in de richting van $S_2^\perp$ (d.w.z. de projectie van $\mathbf{v}$ op $S_1$ ligt in $S_1$, en het residu $\mathbf{v} - P\mathbf{v}$ is loodrecht op $S_2 = \langle b_1,\ldots,b_n\rangle$, dus het residu zit in $S_2^\perp$).

Stel $A = [a_1 \cdots a_n]$ en $B = [b_1 \cdots b_n]$.

**Wat willen we:** Schrijf $P\mathbf{v} = A\mathbf{x}$ voor een onbekend $\mathbf{x}$. Het residu $\mathbf{v} - A\mathbf{x}$ moet loodrecht staan op $S_2$, d.w.z. loodrecht op elke kolom van $B$:

$$B^*(\mathbf{v} - A\mathbf{x}) = \mathbf{0}.$$

**Stap 1:** Uitschrijven geeft:
$$B^* A\mathbf{x} = B^* \mathbf{v}.$$

**Stap 2:** Los op voor $\mathbf{x}$ (mits $B^*A$ inverteerbaar):
$$\mathbf{x} = (B^*A)^{-1}B^*\mathbf{v}.$$

**Stap 3:** De projectie wordt:
$$P_{S_1,S_2}\,\mathbf{v} = A\mathbf{x} = A(B^*A)^{-1}B^*\mathbf{v}.$$

**Conclusie:**
$$\boxed{P_{S_1,S_2} = A(B^*A)^{-1}B^*.} \tag{1}$$

**Verificatie (idempotentie):**
$$P^2 = A(B^*A)^{-1}B^* \cdot A(B^*A)^{-1}B^* = A(B^*A)^{-1}(B^*A)(B^*A)^{-1}B^* = A(B^*A)^{-1}B^* = P. \checkmark$$

**Speciaal geval — orthogonale projector:** Als $B = A$ (de nulruimte is $S_1^\perp$):
$$P_{S_1} = A(A^*A)^{-1}A^*. \tag{2}$$

> ⚠️ **Examentip:** Formule (1) is de algemene scheve projector; (2) is de orthogonale. Ken het verschil en de afleiding.

---

### (c) Bewijs: $P = P^* \Rightarrow$ orthogonale projector (en omgekeerd)

**Te bewijzen:** Een projector $P$ is orthogonaal $\iff$ $P = P^*$.

**Definitie orthogonale projector:** $\mathcal{R}(P) \perp \mathcal{N}(P)$.

#### $(\Rightarrow)$ Als $P$ orthogonaal is, dan $P = P^*$

**Stap 1:** Kies een orthonormale basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$ van $\mathcal{R}(P)$ en $\{\mathbf{q}_{n+1}, \ldots, \mathbf{q}_m\}$ van $\mathcal{N}(P)$.

**Stap 2:** Omdat $\mathcal{R}(P) \perp \mathcal{N}(P)$ is de samengevoegde matrix $Q = [\mathbf{q}_1 \cdots \mathbf{q}_m]$ **unitair**.

**Stap 3:** Bereken $Q^*PQ$. Voor $i \leq n$: $P\mathbf{q}_i = \mathbf{q}_i$ (zit in $\mathcal{R}(P)$). Voor $j > n$: $P\mathbf{q}_j = \mathbf{0}$ (zit in $\mathcal{N}(P)$). Dus:
$$Q^*PQ = \text{diag}(\underbrace{1, \ldots, 1}_{n}, \underbrace{0, \ldots, 0}_{m-n}).$$

**Stap 4:** Deze diagonaalmatrix is reëel en symmetrisch, dus $Q^*PQ = (Q^*PQ)^* = Q^*P^*Q$. Omdat $Q$ unitair: $P = P^*$. $\checkmark$

#### $(\Leftarrow)$ Als $P = P^*$, dan is $P$ orthogonaal

**Neem** $\mathbf{x} = P\mathbf{u} \in \mathcal{R}(P)$ en $\mathbf{y} \in \mathcal{N}(P)$ willekeurig.

$$\langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^*\mathbf{y} = (P\mathbf{u})^*\mathbf{y} = \mathbf{u}^*P^*\mathbf{y} \overset{P^*=P}{=} \mathbf{u}^*(P\mathbf{y}) = \mathbf{u}^* \cdot \mathbf{0} = 0.$$

Dus $\mathcal{R}(P) \perp \mathcal{N}(P)$: $P$ is orthogonaal. $\square$

**Dit zegt ons:** De "zelfgeadjungeerdheid" $P = P^*$ is de algebraïsche vertaling van "de projectie staat loodrecht op het residu."

> ⚠️ **Examentip:** Dit bewijs wordt bijna altijd gevraagd. De kernstap is de berekening $\mathbf{u}^*P^*\mathbf{y} = \mathbf{u}^*P\mathbf{y} = 0$.

---

### (d) De uitmiddelingsmatrix $T$

De matrix $T$ werkt als:
$$T\begin{pmatrix}x_1\\\vdots\\x_m\end{pmatrix} = \frac{x_1 + x_2 + \cdots + x_m}{m}\begin{pmatrix}1\\\vdots\\1\end{pmatrix}.$$

In matrixvorm: $T = \dfrac{1}{m}\mathbf{1}\mathbf{1}^T$, met $\mathbf{1} = (1,1,\ldots,1)^T \in \mathbb{R}^m$.

**Is $T$ een projector?**

Bereken $T^2$:
$$T^2 = \frac{1}{m^2}\mathbf{1}\mathbf{1}^T\mathbf{1}\mathbf{1}^T = \frac{1}{m^2}\mathbf{1}(\mathbf{1}^T\mathbf{1})\mathbf{1}^T = \frac{1}{m^2}\mathbf{1} \cdot m \cdot \mathbf{1}^T = \frac{1}{m}\mathbf{1}\mathbf{1}^T = T.$$

**Ja**, $T^2 = T$, dus $T$ is een projector. ✓

**Waarop projecteert $T$?** Op $\mathcal{R}(T) = \langle \mathbf{1} \rangle$ — de lijn van constante vectoren.

**Is $T$ een orthogonale projector?**

Bereken $T^*$ (voor reëel: $T^T$):
$$T^T = \left(\frac{1}{m}\mathbf{1}\mathbf{1}^T\right)^T = \frac{1}{m}\mathbf{1}\mathbf{1}^T = T.$$

Omdat $T = T^*$ (en $T^2 = T$), is $T$ een **orthogonale projector**. ✓

**Intuïtie:** $T$ projecteert elke vector loodrecht op de deelruimte van constante vectoren — het gemiddelde van de componenten. Het residu $(I-T)\mathbf{x}$ heeft som nul en staat loodrecht op $\mathbf{1}$.

> ⚠️ **Examentip:** De twee berekeningen zijn kort: $T^2 = T$ (projector) en $T = T^T$ (orthogonaal). Verwar niet de twee.

---

## V2 – Householder-transformaties en QR-factorisatie

**Vraag:** Bespreek hoe je met Householder-transformaties een QR-factorisatie van $A \in \mathbb{C}^{m \times n}$ berekent.

**Kernidee:** Householder-transformaties zijn unitaire reflecties die in één stap een volledige kolom op één getal brengen. Door opeenvolgend kolommen te reduceren, brengen we $A$ in bovendriehoekse vorm $R$.

---

### (a) Betekenis van $F_m = I_m - 2\dfrac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}$ — Afleiding

**Vertrekpunt:** De orthogonale projector op het orthogonale complement van $\langle \mathbf{v} \rangle$ is:
$$P_{\langle \mathbf{v}\rangle^\perp} = I - \frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}.$$

Deze projector beeldt een vector af op zijn component loodrecht op $\mathbf{v}$.

**Van projector naar reflector:** Een reflectie t.o.v. een hypervlak loodrecht op $\mathbf{v}$ doet dit: de component langs $\mathbf{v}$ wordt omgekeerd, de rest blijft. Dit is equivalent met twee keer de projectie op $\langle\mathbf{v}\rangle$ aftrekken:
$$F_m = I - 2\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}} = 2P_{\langle\mathbf{v}\rangle^\perp} - I.$$

**Verificatie — unitair:**
$$F_m^* F_m = \left(I - 2\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}}\right)^2 = I - 4\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}} + 4\frac{\mathbf{v}\mathbf{v}^*\mathbf{v}\mathbf{v}^*}{(\mathbf{v}^*\mathbf{v})^2} = I - 4\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}} + 4\frac{\mathbf{v}\mathbf{v}^*}{\mathbf{v}^*\mathbf{v}} = I. \checkmark$$

**Wat doet $F_m$ precies?** Met de keuze $\mathbf{v} = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2\,\mathbf{e}_1$:

**Stap 1:** We willen $F_m \mathbf{x} = -\text{sign}(x_1)\|\mathbf{x}\|_2 \mathbf{e}_1$, d.w.z. alle componenten op nul brengen behalve de eerste.

**Stap 2:** Bereken:
$$F_m\mathbf{x} = \mathbf{x} - 2\frac{\mathbf{v}^*\mathbf{x}}{\mathbf{v}^*\mathbf{v}}\mathbf{v}.$$

**Stap 3:** Met $\mathbf{v} = \mathbf{x} + \alpha\mathbf{e}_1$ (waarbij $\alpha = \text{sign}(x_1)\|\mathbf{x}\|_2$):
$$\mathbf{v}^*\mathbf{x} = \|\mathbf{x}\|_2^2 + \bar{\alpha}x_1, \qquad \mathbf{v}^*\mathbf{v} = \|\mathbf{x}\|_2^2 + 2\text{Re}(\bar{\alpha}x_1) + |\alpha|^2 = 2(\|\mathbf{x}\|_2^2 + \text{Re}(\bar{\alpha}x_1)).$$

**Stap 4:** Na invullen: $F_m\mathbf{x} = -\alpha\mathbf{e}_1 = -\text{sign}(x_1)\|\mathbf{x}\|_2\mathbf{e}_1$. ✓

**Geometrische interpretatie:** $F_m$ spiegelt $\mathbf{x}$ t.o.v. het hypervlak dat loodrecht staat op $\mathbf{v}$, zodat $\mathbf{x}$ precies op de $e_1$-as terechtkomt.

> ⚠️ **Examentip:** De Householder-transformatie is een *reflector* (unitair en zelfinvers: $F_m^2 = I$). Ze brengt een volledige vector op de eerste standaardbasisvector.

---

### (b) Waarom rekening houden met het teken van $x_1$?

**Probleem:** Zonder tekenconventie zijn er twee keuzen: $\mathbf{v} = \mathbf{x} + \|\mathbf{x}\|_2\mathbf{e}_1$ of $\mathbf{v} = \mathbf{x} - \|\mathbf{x}\|_2\mathbf{e}_1$.

**Numeriek probleem:** Kies $\alpha = -\|\mathbf{x}\|_2$ (tegengesteld teken aan $x_1$). Dan:
$$\mathbf{v} = \mathbf{x} - \text{sign}(x_1)\|\mathbf{x}\|_2\mathbf{e}_1.$$

Als $x_1$ bijna gelijk is aan $\|\mathbf{x}\|_2$ (d.w.z. $\mathbf{x}$ wijst bijna in de richting van $\mathbf{e}_1$), dan is $v_1 = x_1 - \|\mathbf{x}\|_2 \approx 0$: **gevaarlijke aftrekking** van bijna gelijke getallen. Dit leidt tot verlies van significante cijfers.

**Oplossing:** Kies het teken zodat $v_1 = x_1 + \text{sign}(x_1)\|\mathbf{x}\|_2$ **optelling** is van gelijknamige getallen — nooit aftrekking van bijna gelijke grootheden. Concreet:
$$\mathbf{v} = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2\mathbf{e}_1.$$

Dan is $v_1 = x_1 + |x_1| \geq |x_1|$ nooit klein in vergelijking met $\|\mathbf{x}\|_2$.

**💡 Inzicht:** Het teken van het resultaat ($\mp\|\mathbf{x}\|_2\mathbf{e}_1$) speelt niet mee voor de QR-factorisatie (enkel de absolute waarde van $r_{11}$ is voorgeschreven); het teken wordt gekozen voor **numerieke stabiliteit**.

> ⚠️ **Examentip:** De reden is altijd "gevaarlijke aftrekking vermijden" — niet meer en niet minder.

---

### (c) Algoritme voor $R$ en $Q$: stap voor stap

**Doel:** Bereken de QR-factorisatie $A = QR$ met $Q$ unitair en $R$ bovendriehoeks.

#### Berekening van $R$

Kolom per kolom brengen we nullen onder de diagonaal.

**Stap $j$ ($j = 1, \ldots, n$):**

1. Beschouw de subkolom $\mathbf{x} = A^{(j-1)}_{j:m, j}$ (elementen $j$ tot $m$ van kolom $j$ van de actuele matrix).
2. Bereken $\mathbf{v}_j = \mathbf{x} + \text{sign}(x_1)\|\mathbf{x}\|_2\mathbf{e}_1$.
3. Definieer $Q_{j,m} = \text{diag}(I_{j-1},\, F_{m-j+1})$ met $F_{m-j+1} = I - 2\dfrac{\mathbf{v}_j\mathbf{v}_j^*}{\mathbf{v}_j^*\mathbf{v}_j}$.
4. Update: $A^{(j)} = Q_{j,m} A^{(j-1)}$.

Na $n$ stappen: $A^{(n)} = Q_{n,m} \cdots Q_{1,m} A = R$ (bovendriehoeks).

**Pseudocode voor $R$:**
```
Input: A ∈ ℂ^{m×n}
for j = 1 to n:
    x = A[j:m, j]                              ← subkolom
    v = x + sign(x[1]) * ||x||₂ * e₁           ← Householder-vector
    A[j:m, j:n] = A[j:m, j:n] - 2*(v*v')/(v'*v) * A[j:m, j:n]
                                                ← efficiënte toepassing
Output: R = A (bovendriehoeks)
```

**Efficiënte toepassing (vermijd expliciete matrix):**
$$Q_{j,m}\mathbf{q} = \mathbf{q} - 2\frac{\mathbf{v}_j^*\mathbf{q}}{\mathbf{v}_j^*\mathbf{v}_j}\mathbf{v}_j \qquad \text{(kost slechts } O(m) \text{ bewerkingen)}$$

#### Berekening van $Q$

Omdat $Q_{n,m} \cdots Q_{1,m} A = R$ en elke $Q_{j,m}$ unitair (en zelfsinvers: $Q_{j,m}^{-1} = Q_{j,m}$):
$$A = Q_{1,m}Q_{2,m} \cdots Q_{n,m}\, R = QR,$$
dus $Q = Q_{1,m}Q_{2,m} \cdots Q_{n,m}$.

In de praktijk bouwt men $Q$ op door achtereenvolgens de inverses toe te passen op $I_m$:

```
Q = I_m
for j = n downto 1:
    Q[j:m, :] = Q[j:m, :] - 2*(v_j * v_j')/(v_j'*v_j) * Q[j:m, :]
```

Of equivalenter: sla de vectoren $\mathbf{v}_j$ op en pas ze later toe.

**Rekenkost:**
- Enkel $R$: $\sim 2mn^2 - \tfrac{2}{3}n^3$ bewerkingen
- Inclusief $Q$: $\sim 4m^2n - \tfrac{2}{3}n^3$ bewerkingen

**Vergelijking:** Householder is numeriek stabiel (achterwaarts stabiel) en goedkoper dan Givens voor algemene matrices.

> ⚠️ **Examentip:** Ken het algoritme in pseudocode en begrijp de efficiënte toepassing via $F_m\mathbf{q} = \mathbf{q} - 2\tfrac{\mathbf{v}^*\mathbf{q}}{\mathbf{v}^*\mathbf{v}}\mathbf{v}$. Dit vermijdt het expliciet opbouwen van de $m \times m$ matrix $F_m$.

---

## V3 – Overgedetermineerd stelsel $Ax \approx b$ via QR

**Vraag:** Beschouw het overgedetermineerde stelsel $Ax \approx b$, $A \in \mathbb{R}^{m \times n}$, $m > n$.

**Kernidee:** Het kleinste-kwadratenprobleem $\min_x \|Ax - b\|_2$ heeft als geometrische interpretatie: zoek de vector $\hat{y} = A\hat{x}$ in $\mathcal{R}(A)$ die het dichtst bij $b$ ligt. De QR-factorisatie geeft een stabielere aanpak dan het normaalstelsel.

---

### (a) Kleinste-kwadratenoplossing via QR-factorisatie

**Vertrekpunt:** We willen $\hat{x}$ zodat $\|Ax - b\|_2$ minimaal is.

**Stap 1:** Bereken de volledige QR-factorisatie $A = QR$ met $Q \in \mathbb{R}^{m \times m}$ orthogonaal en

$$R = \begin{pmatrix}\hat{R} \\ 0\end{pmatrix}, \quad \hat{R} \in \mathbb{R}^{n \times n} \text{ bovendriehoeks, inverteerbaar (als } A \text{ volle kolomrang heeft)}.$$

**Stap 2:** Gebruik dat orthogonale transformaties de 2-norm bewaren ($\|Q^Tz\|_2 = \|z\|_2$):

$$\|Ax - b\|_2^2 = \|Q^T(Ax - b)\|_2^2 = \left\|Q^T A x - Q^T b\right\|_2^2 = \left\|\begin{pmatrix}\hat{R}x \\ 0\end{pmatrix} - \begin{pmatrix}b_1 \\ b_2\end{pmatrix}\right\|_2^2$$

waarbij $\begin{pmatrix}b_1 \\ b_2\end{pmatrix} = Q^T b$ met $b_1 \in \mathbb{R}^n$, $b_2 \in \mathbb{R}^{m-n}$.

**Stap 3:** Uitschrijven:
$$\|Ax - b\|_2^2 = \|\hat{R}x - b_1\|_2^2 + \|b_2\|_2^2.$$

Enkel de eerste term hangt af van $x$. De minimum wordt bereikt wanneer:
$$\hat{R}\hat{x} = b_1. \tag{3}$$

**Stap 4:** Los het bovendriehoekse stelsel (3) op via achterwaartse substitutie.

**De minimale fout is $\|b_2\|_2$** (de restterm die we nooit kunnen elimineren).

**Algoritme:**
1. Bereken onvolledige QR: $A = \hat{Q}\hat{R}$ ($\hat{Q} \in \mathbb{R}^{m \times n}$, $\hat{R} \in \mathbb{R}^{n \times n}$)
2. Bereken $b_1 = \hat{Q}^T b$
3. Los $\hat{R}\hat{x} = b_1$ op via achterwaartse substitutie

> ⚠️ **Examentip:** Stap 2 (norm-invariantie onder orthogonale transformatie) is de sleutel. De splitsing in $b_1$ en $b_2$ maakt het probleem direct oplosbaar.

---

### (b) Onvolledige of volledige QR-factorisatie nodig?

**Antwoord:** De **onvolledige (dunne) QR-factorisatie** $A = \hat{Q}\hat{R}$ volstaat.

**Uitleg:**
- De oplossing hangt enkel af van $b_1 = \hat{Q}^T b$ (de eerste $n$ componenten van $Q^T b$).
- De rest-$b_2$ geeft de minimale fout maar we hebben $b_2$ niet nodig om $\hat{x}$ te berekenen. **Geometrisch:** $b_2 = (I - \hat{Q}\hat{Q}^T)b$ is de projectie van $b$ op het **orthogonale complement van $\mathcal{R}(A)$** — geen keuze van $x$ kan dit deel elimineren, want $Ax \in \mathcal{R}(A)$ per definitie. Het is de intrinsieke, onvermijdelijke restfout van het benaderingsprobleem.
- Bijgevolg: men berekent enkel $\hat{Q} \in \mathbb{R}^{m \times n}$ (de eerste $n$ kolommen van $Q$) en $\hat{R} \in \mathbb{R}^{n \times n}$.

**Voordeel onvolledige QR:** Rekenkost $O(mn^2)$ i.p.v. $O(m^2 n)$ voor de volledige factorisatie.

💡 **Inzicht:** De volledige QR bevat extra informatie (de restfout $b_2$) die we wel kunnen berekenen als we de minimale fout ook willen weten, maar die strikt niet nodig is voor $\hat{x}$.

---

### (c) $A$ is $10^7 \times 10$ — welk algoritme?

**Antwoord:** **Gram-Schmidt** (bij voorkeur het gewijzigd Gram-Schmidt algoritme).

**Motivering:**

| Eigenschap | Gram-Schmidt | Householder | Givens |
|---|---|---|---|
| Rekenkost voor $R$ | $O(2mn^2)$ | $O(2mn^2)$ | $O(3mn^2)$ |
| Rekenkost voor $Q$ expliciet | $O(2mn^2)$ | $O(4m^2n)$ | $O(6m^2n)$ |
| Kolommen op voorhand nodig? | Nee | Ja | Ja |
| Geheugen | $O(mn)$ | $O(m^2)$ (voor volle $Q$) | $O(m^2)$ |

**Voor $m = 10^7$, $n = 10$:**
- Householder: $Q \in \mathbb{R}^{10^7 \times 10^7}$ — **totaal onhaalbaar** qua geheugen.
- Gram-Schmidt: bouwt enkel $\hat{Q} \in \mathbb{R}^{10^7 \times 10}$ op — $10^7 \times 10$ getallen, dat is haalbaar.
- Bovendien: voor $m \gg n$ (zwaar overgedetermineerd) heeft Gram-Schmidt precies dezelfde asymptotische kost $O(2mn^2)$ voor de onvolledige factorisatie.

**Rekenkost vergelijking:**
- Gram-Schmidt ($\hat{Q}$ + $\hat{R}$): $\approx 2 \times 10^7 \times 100 = 2 \times 10^9$ bewerkingen
- Householder (volledige $Q$): $\approx 4 \times (10^7)^2 \times 10 = 4 \times 10^{15}$ bewerkingen — factoor $10^6$ duurder!

> ⚠️ **Examentip:** Het antwoord is altijd Gram-Schmidt wanneer $m \gg n$, met als reden het geheugenverbruik en de kost van de volledige $Q$ bij Householder.

---

### (d) Waarom beter dan rechtstreeks het normaalstelsel oplossen?

Het normaalstelsel is $A^T A \hat{x} = A^T b$.

**Probleem met het normaalstelsel:**

$$\kappa(A^T A) = \kappa(A)^2.$$

**Gevolg:** Als $\kappa(A) = 100$, dan is $\kappa(A^TA) = 10^4$. Een kleine fout op $b$ of $A$ wordt $10^4 \times$ vergroot in de oplossing.

**Relatieve fout bij normaalstelsel:**
$$\frac{\|\tilde{x} - \hat{x}\|_2}{\|\hat{x}\|_2} = O\!\left(\kappa(A)^2 \cdot \epsilon_{\text{mach}}\right).$$

**Relatieve fout bij QR-aanpak:**
$$\frac{\|\tilde{x} - \hat{x}\|_2}{\|\hat{x}\|_2} = O\!\left(\kappa(A) \cdot \epsilon_{\text{mach}}\right).$$

De QR-aanpak is een factor $\kappa(A)$ nauwkeuriger — voor slecht geconditioneerde matrices kan dit het verschil zijn tussen een zinvolle en een volledig foutieve oplossing.

**Intuïtie:** $A^T A$ kwadreert de conditioneringsverhouding tussen grootste en kleinste singuliere waarden van $A$. De QR-factorisatie werkt rechtstreeks op $A$ en vermijdt dit kwadreren.

> ⚠️ **Examentip:** Onthoud de formule $\kappa(A^TA) = \kappa(A)^2$ en de vergelijking van de fouten. Dit is een standaardargument dat ook buiten dit vak terugkomt.

---

## V4 – De maximumnorm op $C[a,b]$

**Vraag:** Beschouw de functionaal $I(f) = \max_{x \in [a,b]} |f(x)|$.

**Kernidee:** De maximumnorm is de meest intuïtieve norm op functies. Ze meet de "grootste afwijking" — relevant voor minimax-benaderingen. Ze is echter niet strict en voldoet niet aan de parallellogramgelijkheid, dus ze komt niet van een scalair product.

---

### (a) Bewijs dat $I(f)$ een norm is op $C[a,b]$

We moeten de vier normeigenschappen verifiëren.

**Eigenschap 1 — Niet-negativiteit:** $|f(x)| \geq 0$ voor alle $x$, dus $\max|f(x)| \geq 0$. ✓

**Eigenschap 2 — Nulvector:** $I(f) = 0 \iff \max_x|f(x)| = 0 \iff |f(x)| = 0$ voor alle $x \iff f \equiv 0$.
→ Wél: het maximum van een *continue* functie op een *gesloten begrensd interval* is bereikbaar (stelling van Weierstrass). ✓

**Eigenschap 3 — Homogeniteit:** Voor $\alpha \in \mathbb{R}$:
$$I(\alpha f) = \max_x |\alpha f(x)| = |\alpha|\max_x|f(x)| = |\alpha|\,I(f). \checkmark$$

**Eigenschap 4 — Driehoeksongelijkheid:** Voor alle $x \in [a,b]$:
$$|f(x) + g(x)| \leq |f(x)| + |g(x)| \leq \max_x|f(x)| + \max_x|g(x)| = I(f) + I(g).$$

Dit geldt voor alle $x$, dus ook voor de maximumplaats:
$$I(f + g) = \max_x|f(x) + g(x)| \leq I(f) + I(g). \checkmark$$

**Conclusie:** $I(f) = \|f\|_\infty$ is een norm op $C[a,b]$. $\square$

> ⚠️ **Examentip:** Het bewijsje van eigenschap 4 is de meest gevraagde stap. Schrijf het expliciet uit.

---

### (b) Is het een strikte norm?

**Definitie strikte norm:** Een genormeerde ruimte is *strikt* als voor $\|f\| = \|g\| = 1$ en $f \neq g$:
$$\|f + g\| < 2.$$

Equivalenter: gelijkheid in de driehoeksongelijkheid geldt alleen als één van de twee een niet-negatief veelvoud is van de ander.

**Antwoord: Nee**, de maximumnorm is **niet** strikt.

**Tegenvoorbeeld:** Kies $[a,b] = [0,1]$, $f(x) = 1$ en $g(x) = x$. Dan:
$$\|f\|_\infty = 1, \quad \|g\|_\infty = 1, \quad \|f + g\|_\infty = \max_{x \in [0,1]}|1 + x| = 2.$$

Maar $f$ is geen veelvoud van $g$, toch geldt $\|f + g\|_\infty = \|f\|_\infty + \|g\|_\infty = 2$.

**Gevolg:** De beste benadering in de maximumnorm is niet altijd uniek (de existentie-uniciteitsresultaten die *strikte* normen vereisen gelden hier niet).

💡 **Inzicht:** De maximumnorm heeft "rechte stukken" op haar eenheidsbol — overal waar twee functies tegelijk hun maximum bereiken. In tegenstelling tot de $L^2$-norm (die stikt convex is).

---

### (c) Voldoet deze norm aan de parallellogramgelijkheid?

**Parallellogramgelijkheid:** $\|f + g\|^2 + \|f - g\|^2 = 2(\|f\|^2 + \|g\|^2)$.

**Antwoord: Nee.**

**Tegenvoorbeeld:** $[a,b] = [0,1]$, $f(x) = 1$, $g(x) = x$.
$$\|f+g\|_\infty = 2, \quad \|f-g\|_\infty = 1, \quad \|f\|_\infty = 1, \quad \|g\|_\infty = 1.$$

Linkerlid: $4 + 1 = 5$. Rechterlid: $2(1 + 1) = 4$. Ongelijk! $5 \neq 4$.

**Gevolg:** De maximumnorm op $C[a,b]$ komt **niet voort uit een scalair product**. Er bestaat geen scalair product $(\cdot,\cdot)$ op $C[a,b]$ zodat $\sqrt{(f,f)} = \|f\|_\infty$.

**Vergelijk:** De $L^2$-norm $\|f\|_2 = \sqrt{\int_a^b f^2\,dx}$ voldoet wel aan de parallellogramgelijkheid — zij is de geïnduceerde norm van het scalair product $(f,g) = \int_a^b f(x)g(x)\,dx$.

> ⚠️ **Examentip:** De parallellogramgelijkheid is het criterium om te beslissen of een norm van een scalair product afkomt (stelling van Jordan-von Neumann). Het tegenvoorbeeld met $f=1$, $g=x$ werkt altijd voor de maximumnorm.

---

## V5 – Gram-Schmidt in unitaire ruimte

**Vraag:** Beschouw een unitaire ruimte met scalair product $(\cdot,\cdot)$.

**Kernidee:** Gram-Schmidt veralgemening van de QR-factorisatie naar abstracte ruimten van functies. Inwendige producten worden scalaire producten; Euclidische norm wordt geïnduceerde norm. Maar Householder/Givens werken hier niet meer.

---

### (a) Gram-Schmidt algoritme + verband met orthogonale projecties

**Opgave:** Gegeven een basis $\{a_1,\ldots,a_n\}$ van een deelruimte $\mathcal{D} \subset V$, construeer een **geneste orthonormale basis** $\{q_1,\ldots,q_n\}$:
$$\langle q_1,\ldots,q_k\rangle = \langle a_1,\ldots,a_k\rangle, \quad k = 1,\ldots,n.$$

**Algoritme (Klassiek Gram-Schmidt):**
```
voor j = 1 tot n:
    v_j = a_j
    voor i = 1 tot j-1:
        r_{ij} = (q_i, a_j)         ← scalair product
        v_j = v_j - r_{ij} · q_i   ← verwijder component langs q_i
    r_{jj} = ||v_j||                ← geïnduceerde norm
    q_j = v_j / r_{jj}             ← normaliseer
```

**Verband met orthogonale projecties:** In stap $j$ berekenen we:
$$v_j = a_j - P_{\langle q_1,\ldots,q_{j-1}\rangle}\,a_j,$$

waarbij $P_{\langle q_1,\ldots,q_{j-1}\rangle}$ de orthogonale projector is op de reeds gevonden deelruimte:
$$P_{\langle q_1,\ldots,q_{j-1}\rangle}\,a_j = \sum_{i=1}^{j-1}(q_i, a_j)\,q_i.$$

We trekken de **loodrechte projectie** van $a_j$ op de reeds gebouwde ruimte af. Het resultaat $v_j$ is dan orthogonaal op alle $q_1,\ldots,q_{j-1}$. Na normalisatie geeft dit $q_j$.

**Relatie met QR-factorisatie:** De coëfficiënten $r_{ij}$ vormen een bovendriehoekse matrix $\hat{R}$, en de constructie geeft de abstracte analogon van $A = \hat{Q}\hat{R}$:
$$a_j = \sum_{i=1}^j r_{ij}\,q_i.$$

> ⚠️ **Examentip:** De geometrische interpretatie ("projecteer en trek af") is essentieel. $v_j$ staat loodrecht op de vorige ruimte precies omdat we de projectie aftrekken.

---

### (b) Hoeveelheid rekenwerk: aantal scalaire producten

**Stap $j$** vereist:
- $j-1$ scalaire producten $(q_i, a_j)$ voor $i = 1,\ldots,j-1$
- $1$ scalaire product voor de norm $\|v_j\|^2 = (v_j, v_j)$

**Totaal voor $n$ stappen:**
$$\sum_{j=1}^n j = \frac{n(n+1)}{2} = O(n^2) \text{ scalaire producten.}$$

**Vergelijk met drie directe benaderingsmethodes:**
- **Normaalstelsel $G\mathbf{c} = \mathbf{b}$:** ook $O(n^2)$ scalaire producten voor het vullen van $G$ en $\mathbf{b}$
- Maar het normaalstelsel is numeriek gevaarlijker door $\kappa(G) = \kappa(A)^2$

> ⚠️ **Examentip:** De drietermsrecursiebetrekking voor orthogonale veeltermen reduceert dit tot $O(n)$ scalaire producten — een enorme besparing bij de klasse van veeltermen (zie V6).

---

### (c) Herorthogonalisatie

**Wat is het?** Na afloop van het Gram-Schmidt algoritme past men het algoritme **een tweede keer** toe op de bekomen vectoren $\{q_1,\ldots,q_n\}$. Men orthogonaliseert de uitvoer tegenover zichzelf.

**Wanneer nodig?** Gram-Schmidt is wiskundig correct maar **numeriek niet achterwaarts stabiel**. Wanneer de originele basis $\{a_1,\ldots,a_n\}$ slecht geconditioneerd is (vectoren liggen bijna in dezelfde richting), accumuleren afrondingsfouten en raken de berekende $q_j$'s hun orthogonaliteit kwijt:
$$(q_i, q_j) \approx \epsilon_{\text{mach}} \cdot \kappa(A) \quad \text{(verlies van orthogonaliteit)}.$$

**Effect van herorthogonalisatie:** Na één extra ronde Gram-Schmidt worden de afwijkingen van de orde $\epsilon_{\text{mach}}^2 \cdot \kappa(A)^2$, wat in de meeste praktische situaties te verwaarlozen is.

**Vergelijk:** Householder-transformaties (voor vectoren in $\mathbb{R}^m$) zijn achterwaarts stabiel en vereisen nooit herorthogonalisatie. Maar voor abstracte unitaire ruimten is Householder niet toepasbaar (zie §d).

> ⚠️ **Examentip:** "Herorthogonalisatie" = "Gram-Schmidt nogmaals uitvoeren op de uitvoer." De reden: Gram-Schmidt is niet achterwaarts stabiel bij slecht geconditioneerde bases.

#### Gewijzigd Gram-Schmidt (Modified Gram-Schmidt, MGS) — numeriek stabielere variant

**Verschil met klassiek GS:** Bij klassiek GS worden alle componenten $r_{ij} = (q_i, a_j)$ berekend t.o.v. de **originele** vector $a_j$. Bij MGS worden de componenten **stapsgewijs** verwijderd: na het aftrekken van de $q_i$-component, wordt de **bijgewerkte** $v_j$ gebruikt voor de volgende orthogonalisatie.

**Algoritme (Gewijzigd Gram-Schmidt):**

```text
voor j = 1 tot n:
    v_j = a_j
    voor i = 1 tot j-1:
        r_{ij} = (q_i, v_j)         ← inproduct met bijgewerkte v_j (niet a_j!)
        v_j = v_j - r_{ij} · q_i   ← onmiddellijk verwijderen
    r_{jj} = ||v_j||
    q_j = v_j / r_{jj}
```

**Waarom stabieler?** Bij KGS propageren kleine afrondingsfouten in elke stap naar alle volgende stappen (de fout op $a_j$ beïnvloedt alle $r_{ij}$ tegelijk). Bij MGS wordt de fout na elke deelstap "opgeschoond" — de gecumuleerde fout op de orthogonaliteit is van orde $\epsilon_{\text{mach}} \cdot \kappa(A)$ in plaats van $\epsilon_{\text{mach}} \cdot \kappa(A)^2$.

**Wiskundig identiek, numeriek beter:** KGS en MGS geven in exacte rekenkunde identieke resultaten. Bij drijvende-kommarekenkunde is MGS significant robuuster, en vervangt in de praktijk steeds het klassieke Gram-Schmidt.

> ⚠️ **Examentip:** Michiels vraagt dit altijd: het woord **Modified Gram-Schmidt (MGS)** moet vallen. Het sleutelwoord is "stapsgewijs verwijderen van componenten t.o.v. de bijgewerkte vector", in tegenstelling tot "alles t.o.v. de originele vector" bij KGS.

---

### (d) Orthogonale veeltermen opstellen met Gram-Schmidt: efficiënt?

**Vraag:** Stel dat we orthogonale veeltermen willen t.o.v. $(f,g) = \int_a^b w(x)f(x)g(x)\,dx$, vertrekkende van $\{1, x, x^2, \ldots\}$.

**Antwoord: Nee, dit is niet de computationeel meest efficiënte methode.**

**Waarom niet?**

Gram-Schmidt kost $O(n^2)$ scalaire producten. Elk scalair product is een integraal:
$$(q_i, x^j) = \int_a^b w(x)\,q_i(x)\,x^j\,dx.$$

**Betere aanpak: drietermsrecursiebetrekking (zie V6).**

Voor orthogonale veeltermen geldt dat elke nieuwe veelterm $\phi_k$ enkel afhangt van $\phi_{k-1}$ en $\phi_{k-2}$:
$$\phi_k(x) = \lambda_k\!\left[(x - \alpha_k)\phi_{k-1}(x) - \beta_k\phi_{k-2}(x)\right].$$

Dit vraagt slechts **2 scalaire producten per stap** (voor $\alpha_k$ en $\beta_k$): in totaal $O(n)$ scalaire producten voor de volledige rij.

**Verklaring:** De speciale structuur van veeltermen (scalair product met $x\phi_{k-1}$ is enkel niet-nul voor twee buurpunten in de rij) maakt de $O(n^2)$ matrix $G$ effectief $O(n)$ — een eigenschap die Gram-Schmidt niet benut.

💡 **Inzicht:** Gram-Schmidt is een generiek algoritme. Voor de speciale klasse van veeltermen is een slimmere, op de structuur gebaseerde methode beschikbaar.

---

## V6 – Orthogonale rij veeltermen $\{\phi_0, \phi_1, \phi_2, \ldots\}$

**Vraag:** Beschouw een orthogonale rij veeltermen t.o.v. $(f,g) = \int_a^b w(x)f(x)g(x)\,dx$, $w(x) > 0$.

**Kernidee:** Orthogonale veeltermen hebben fascinerende eigenschappen: hun nulpunten zijn enkelvoudig, reëel en liggen in het open interval. Die nulpunten zijn tegelijk eigenwaarden van een tridiagonale matrix — wat ze numeriek berekenbaar maakt. Interpolatie in die nulpunten is bijna optimaal.

---

### (a) Nulpunten van $\phi_n$

**Stelling:** De veelterm $\phi_n(x)$ heeft precies **$n$ enkelvoudige reële nulpunten**, alle gelegen in het **open interval** $(a, b)$.

**Argumenten (geen bewijs vereist):**
1. **Reëel:** Alle coëfficiënten van $\phi_n$ zijn reëel (het scalair product is reëel).
2. **Precies $n$:** $\phi_n$ is een veelterm van graad $n$, dus hoogstens $n$ nulpunten. De orthogonaliteit t.o.v. constanten ($(\phi_n, 1) = 0$) dwingt minstens één nulpunt in $(a,b)$. Een vollediger argument toont dat er minstens $n$ tekenwisselingen zijn.
3. **Enkelvoudig:** Geen meervoudige nulpunten (anders zou het product $\phi_n \cdot q$ voor geschikte $q$ niet van teken wisselen maar het integraal toch nul moeten zijn — contradictie).
4. **Open interval:** De randpunten $a$ en $b$ zijn nooit nulpunten (volgt uit de normalisatiekeuze en de positieve gewichtsfunctie).

> ⚠️ **Examentip:** Weet dit te formuleren. Het bewijs hoeft niet volledig, maar weet dat het via een tegenspraakargument gaat met het teken van $\phi_n \cdot q$.

---

### (b) Bewijs: nulpunten zijn eigenwaarden van een tridiagonale matrix

**Te bewijzen:** De $n$ nulpunten van $\phi_n$ zijn de eigenwaarden van de $n \times n$ matrix:
$$A = \begin{pmatrix}
\alpha_1 & \mu_1 & & \\
\nu_1 & \alpha_2 & \mu_2 & \\
& \nu_2 & \alpha_3 & \ddots \\
& & \ddots & \ddots & \mu_{n-1} \\
& & & \nu_{n-1} & \alpha_n
\end{pmatrix}.$$

**Vertrekpunt:** De **drietermsrecursiebetrekking** voor de orthogonale rij luidt (herschreven zoals gegeven):
$$\nu_{k-1}\phi_{k-2}(x) + \alpha_k\phi_{k-1}(x) + \mu_k\phi_k(x) = x\phi_{k-1}(x), \quad k = 1,\ldots,n.$$

*(Noot: vergelijkt met de standaardvorm $\phi_k(x) = \lambda_k[(x-\alpha_k)\phi_{k-1}(x) - \beta_k\phi_{k-2}(x)]$, zijn $\mu_k, \nu_k, \alpha_k$ herscalings van de recursiecoëfficiënten.)*

**Stap 1:** Definieer de vector $\Phi(x) = [\phi_0(x), \phi_1(x), \ldots, \phi_{n-1}(x)]^T$.

**Stap 2:** Schrijf de recursiebetrekking voor $k = 1, \ldots, n-1$ als een matrixvergelijking:
$$A\Phi(x) = x\Phi(x) - \mu_n\phi_n(x)\,\mathbf{e}_n,$$

waarbij $\mathbf{e}_n$ de $n$-de standaardbasisvector is.

**Verificatie per component:**
- Rij $k$ ($1 \leq k \leq n-1$): geeft $\nu_{k-1}\phi_{k-2}(x) + \alpha_k\phi_{k-1}(x) + \mu_k\phi_k(x) = x\phi_{k-1}(x)$ ✓
- Rij $n$: geeft $\nu_{n-1}\phi_{n-2}(x) + \alpha_n\phi_{n-1}(x) = x\phi_{n-1}(x) - \mu_n\phi_n(x)$ ✓

**Stap 3:** Evalueer in een nulpunt $x_j$ van $\phi_n$, d.w.z. $\phi_n(x_j) = 0$. Dan valt de stoorterm weg:
$$A\Phi(x_j) = x_j\Phi(x_j).$$

**Conclusie:** $x_j$ is een eigenwaarde van $A$ met eigenvector $\Phi(x_j) = [\phi_0(x_j), \phi_1(x_j), \ldots, \phi_{n-1}(x_j)]^T$.

**Stap 4:** Er zijn $n$ nulpunten van $\phi_n$ (Stelling uit §a), en $A$ heeft hoogstens $n$ eigenwaarden. Elke eigenwaarde van $A$ is een nulpunt (omgekeerde richting volgt door rang-argumenten). Dus de nulpunten en eigenwaarden zijn precies dezelfde $n$ getallen. $\square$

**Praktisch belang:** De nulpunten van $\phi_n$ worden berekend door het eigenwaardenprobleem van de $n \times n$ tridiagonale matrix $A$ op te lossen via het QR-algoritme — efficiënt, nauwkeurig, en in $O(n^2)$ bewerkingen.

> ⚠️ **Examentip:** Stap 3 is de kern: door te evalueren in een nulpunt van $\phi_n$ verdwijnt de stoorterm en krijg je een eigenwaardevergelijking. Het eigenvectorbewijs is elegant maar kort.

---

### (c) Interpolatie in nulpunten van Chebyshev: goede aanpak?

**Context:** We willen een functie $f$ benaderen door een veelterm van graad $n$. We bepalen de coëfficiënten via **interpolatie in de nulpunten van de Chebyshev-veelterm** $T_{n+1}$.

**Antwoord: Ja, dit is een uitstekende aanpak.**

**Argumenten:**

**1. Verband met kleinste-kwadratenbenadering:**
De fout van de kleinste-kwadratenbenadering van graad $n$ met Chebyshev-gewicht $w(x) = 1/\sqrt{1-x^2}$ wordt benaderd door:
$$r_n(x) \approx -a_{n+1}T_{n+1}(x).$$

De nulpunten van $T_{n+1}$ zijn precies de plaatsen waar de fout nul is. Interpolatie in die punten geeft een benadering die bijna identiek is aan de kleinste-kwadratenbenadering.

**2. Nabijheid aan de minimax-benadering:**
Chebyshev-veeltermen $T_k$ oscilleren equi-amplitudinaal tussen $\pm 1$ over $[-1,1]$ ($n+2$ extrema). De fout van interpolatie in Chebyshev-punten voldoet daarmee bijna aan de equioscillatiestelling — het is bijna een minimax-benadering.

**3. Vermijding van Runge-fenomeen:**
Bij equidistante interpolatiepunten kan de fout exploderen voor hoge graden (Runge-fenomeen). Chebyshev-punten concentreren zich aan de randen van het interval, wat dit vermijdt.

**4. Concreet voorbeeld ($e^x$ op $[-1,1]$, graad 4):**

| Methode | Max. fout |
|---|---|
| Taylor (Maclaurin) | $\approx 10^{-2}$ |
| Equidistante interpolatie | $\approx 10^{-3}$ |
| Interpolatie in Chebyshev-punten | $\approx 3 \times 10^{-4}$ |
| Kleinste kwadraten (Chebyshev-gewicht) | $\approx 3 \times 10^{-4}$ |

**De Chebyshev-punten zijn:**
$$x_k = \cos\!\left(\frac{(2k-1)\pi}{2(n+1)}\right), \quad k = 1,\ldots,n+1.$$

> ⚠️ **Examentip:** Drie kernpunten: (1) bijna gelijk aan kleinste-kwadratenbenadering, (2) bijna minimax, (3) vermijdt Runge-fenomeen. Alle drie komen van de equi-oscillatorische eigenschap van Chebyshev-veeltermen.

---

## V7 – Grafen: directe verbindingen en meest centrale knoop

**Vraag:** Beschouw een netwerk beschreven door een enkelvoudige graaf. Leid een uitdrukking af voor het aantal directe verbindingen vanuit de knopen, in termen van de verbindingsmatrix. Geef een definitie van de meest centrale knoop die zich niet beperkt tot directe verbindingen. Hoe kan men deze berekenen?

**Kernidee:** De verbindingsmatrix $A$ codeert de directe verbindingen. Machten van $A$ tellen lussen van grotere lengtes — en de matrixexponentiaal $e^A$ som al die bijdragen samen tot een definitie van centraliteit die het volledige netwerkperspectief omvat.

---

### Directe verbindingen via de verbindingsmatrix

**Enkelvoudige ongerichte graaf:** $N$ knopen, maximaal één boog per paar, geen zelfverbindingen.

**Verbindingsmatrix** $A \in \{0,1\}^{N \times N}$: symmetrisch, $A_{ij} = 1$ als knopen $i$ en $j$ verbonden zijn, 0 anders. Diagonaal is nul.

**Aantal directe verbindingen vanuit knoop $i$:**
$$\text{graad}(i) = \sum_{j=1}^N A_{ij}.$$

**Uitdrukking via $A^2$:** Merk op dat:
$$\sum_j A_{ij}A_{ji} = (A^2)_{ii}.$$

Omdat $A$ symmetrisch is ($A_{ji} = A_{ij}$): $(A^2)_{ii} = \sum_j A_{ij}^2 = \sum_j A_{ij}$ = graad$(i)$.

Dus: **het aantal directe verbindingen vanuit knoop $i$ is $(A^2)_{ii}$**.

**Veralgemeening — lussen van lengte $n$:**

### Stelling
Het aantal lussen van lengte $n$ die starten en eindigen in knoop $i$ is $(A^n)_{ii}$.

**Bewijs (inductie):** $(A^n)_{ij}$ telt het aantal paden van lengte $n$ van knoop $i$ naar knoop $j$. Voor $i = j$ geeft dit het aantal lussen. $\square$

**Intuïtie:** Een knoop $i$ met veel lussen van lengte 2, 3, 4, ... is goed ingebed in het netwerk — niet alleen directe buren, maar ook indirecte verbindingen tellen.

---

### Definitie meest centrale knoop

**Probleem met enkel directe verbindingen:** Een knoop kan veel buren hebben die onderling weinig verbonden zijn. Een ander knooppunt met minder directe buren maar verbonden met goed-verbonden knopen is in de praktijk centraler.

**Oplossing: combineer lussen van alle lengtes** met afnemend gewicht voor grotere lengtes:

$$\text{centraliteit}(i) = \sum_{n=0}^\infty \frac{1}{n!}(A^n)_{ii} = (e^A)_{ii}.$$

De factor $1/n!$ zorgt dat lange lussen minder wegen (en dat de reeks convergeert).

### Definitie 8.2.1 — Centraliteit
$$\text{centraliteit van knoop } i = (e^A)_{ii}.$$

**De meest centrale knoop** is de knoop waarvoor $(e^A)_{ii}$ maximaal is.

**Toelichting op de bijdragen:**
- $n=0$: $(A^0)_{ii} = (I)_{ii} = 1$ — triviale bijdrage
- $n=1$: $(A)_{ii} = 0$ — geen zelfverbindingen
- $n=2$: $(A^2)_{ii} = $ graad$(i)$ — directe buren
- $n=3$: $(A^3)_{ii}$ — driehoeken door $i$
- $n=4$: $(A^4)_{ii}$ — vierhoeken door $i$, enz.

---

### Berekening van de centraliteit

**Stap 1:** Bereken de **eigenwaardenontbinding** van $A$:
$$A = X\Lambda X^{-1}, \quad \Lambda = \text{diag}(\lambda_1,\ldots,\lambda_N).$$

Omdat $A$ symmetrisch is, zijn de eigenwaarden reëel en is $X$ orthogonaal ($X^{-1} = X^T$).

**Stap 2:** Bereken de matrixexponentiaal via de eigenwaardenontbinding (**spectraal mapping theorema**: als $A = X\Lambda X^{-1}$ diagonaliseerbaar is, dan geldt $f(A) = Xf(\Lambda)X^{-1}$ voor elke functie $f$, waarbij $f$ componentgewijs op de eigenwaarden wordt toegepast):
$$e^A = Xe^\Lambda X^{-1} = X\begin{pmatrix}e^{\lambda_1}&&\\&\ddots&\\&&e^{\lambda_N}\end{pmatrix}X^T.$$

**Stap 3:** Lees de diagonaalelementen $(e^A)_{ii}$ af.

**Algoritme:** Het QR-algoritme (Hoofdstuk 9) berekent de volledige eigenwaardenontbinding van $A$ in $O(N^3)$ bewerkingen.

**Voorbeeld:** Voor een graaf van 8 knopen met verbindingsmatrix $A$:
- $(A^2)_{ii}$: knoop 5 heeft 4 directe buren (maximum).
- $(e^A)_{ii}$: knopen 2 en 4 zijn het meest centraal — zij zijn verbonden met andere goed-verbonden knopen; de lussen van lengte 3, 4, ... geven de doorslag.

> ⚠️ **Examentip:** Ken de definitie $(e^A)_{ii}$ en de berekening via eigenwaardenontbinding. Weet dat $e^A = Xe^\Lambda X^{-1}$ — dit is de matrix-analogon van de scalaire exponentiaalfunctie.

---

## V8 – Het algoritme van Arnoldi

**Vraag:** Bespreek het algoritme van Arnoldi.

**Kernidee:** Arnoldi is een Krylov-methode voor het benaderen van extremale eigenwaarden van grote, ijle matrices. Door Gram-Schmidt te combineren met herhaalde matrix-vectorvermenigvuldiging, bouwt het een orthonormale basis van de Krylov-ruimte op en projecteert het eigenwaardenprobleem op die lage-dimensionale ruimte.

---

### (a) Voor welk soort problemen? Waarom? Voorbeeld in context van grafen.

**Geschikt voor:** Het berekenen van **een paar extremale eigenwaarden** van grote, ijle matrices $A \in \mathbb{C}^{m \times m}$, met $m$ groot ($\sim 10^3$–$10^9$).

**Waarom niet gewoon QR-algoritme?** Het volledige QR-algoritme berekent alle $m$ eigenwaarden in $O(m^3)$ bewerkingen — bij $m = 10^6$ volledig onhaalbaar. Arnoldi berekent slechts $k \ll m$ eigenwaarden in $O(km^2)$ bewerkingen (gedomineerd door $k$ matrix-vector producten).

**Cruciale vereiste:** Alleen het product $A\mathbf{v}$ (matrix-vector) hoeft beschikbaar te zijn — de matrix $A$ zelf hoeft nooit volledig in geheugen te zitten. Voor grote ijle matrices is $A\mathbf{v}$ in $O(nnz)$ bewerkingen (waar $nnz$ het aantal niet-nul elementen is).

**Voorbeeld: grafen en PageRank.**
- De linkmatrix $\hat{A} \in \mathbb{R}^{N \times N}$ voor het web heeft $N \sim 10^{10}$ knopen maar is ijl (elke pagina heeft gemiddeld $\sim 10$ links).
- We zoeken slechts de **dominante eigenvector** (eigenwaarde $\lambda = 1$): de PageRank-vector.
- Het Arnoldi-algoritme (of zijn variante: de methode van de machten) lost dit op met slechts matrix-vector producten — de matrix zelf hoeft niet volledig in geheugen.

**Ander voorbeeld:** Graafpartitionering vereist de Fiedler-vector $\mathbf{v}_{N-1}$ van de Laplaciaan $L$ — opnieuw een extreme eigenvector van een grote ijle matrix.

> ⚠️ **Examentip:** Twee sleutelwoorden: "grote ijle matrix" en "slechts een paar extreme eigenwaarden." Het voordeel is dat enkel matrix-vector producten nodig zijn.

---

### (b) Orthonormale basis van de Krylov-ruimte

**Krylov-ruimte:**
$$\mathcal{K}_k(A,\mathbf{b}) = \langle \mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\rangle.$$

**Probleem met de Krylov-vectoren als basis:** Voor grote $k$ worden $A^{k-1}\mathbf{b}$ bijna allemaal evenredig met de dominante eigenvector — gevaarlijke lineaire afhankelijkheid.

**Oplossing:** Bouw via Gram-Schmidt een orthonormale basis $\{\mathbf{q}_1,\ldots,\mathbf{q}_k\}$ van $\mathcal{K}_k$ op.

**Arnoldi Algoritme:**
```
q₁ = b / ||b||₂
voor j = 1, 2, 3, ...:
    v_j = A q_j                           ← matrix-vector product
    voor i = 1 tot j:
        h_{ij} = q_i* v_j                 ← scalair product
        v_j = v_j - h_{ij} q_i            ← Gram-Schmidt
    h_{j+1,j} = ||v_j||₂
    q_{j+1} = v_j / h_{j+1,j}            ← normaliseer
```

**Arnoldi-recursievergelijking:** Na $k$ stappen:
$$\boxed{AQ_k = Q_k H_k + h_{k+1,k}\,\mathbf{q}_{k+1}\mathbf{e}_k^T} \tag{9.8}$$

met:
- $Q_k = [\mathbf{q}_1 \cdots \mathbf{q}_k]$: matrix met orthonormale kolommen
- $H_k \in \mathbb{C}^{k \times k}$: **Hessenberg-matrix** (bovendriehoeks + één subdiagonaal), $H_k = Q_k^* A Q_k$
- $h_{k+1,k}$: de subdiagonaalcoëfficiënt (residu na Gram-Schmidt)

**Afleiding van $H_k = Q_k^* AQ_k$:** Links vermenigvuldigen van (9.8) met $Q_k^*$:
$$Q_k^* AQ_k = Q_k^* Q_k H_k + h_{k+1,k}(Q_k^*\mathbf{q}_{k+1})\mathbf{e}_k^T = H_k + \mathbf{0} = H_k. \checkmark$$

> ⚠️ **Examentip:** Ken de recursievergelijking (9.8) vanbuiten en weet dat $H_k = Q_k^* AQ_k$. De Hessenberg-structuur van $H_k$ is een directe consequentie van het Gram-Schmidt-proces.

---

### (c) Eigenwaarden benaderen: Ritz-waarden

**Idee:** In stap $k$ zoeken we een benadering $(\theta, \mathbf{q})$ van een eigenwaardepaar van $A$, met $\mathbf{q} \in \mathcal{K}_k$, d.w.z. $\mathbf{q} = Q_k\mathbf{y}$.

**Eis:** Het residu $A\mathbf{q} - \theta\mathbf{q}$ moet orthogonaal zijn op $\mathcal{K}_k$:
$$Q_k^*(A - \theta I)Q_k \mathbf{y} = \mathbf{0} \implies H_k\mathbf{y} = \theta\mathbf{y}.$$

**Definitie:** De **Ritz-waarden** $\theta_1^{(k)},\ldots,\theta_k^{(k)}$ zijn de eigenwaarden van de $k \times k$ Hessenberg-matrix $H_k$.

**Verband met orthogonale projectie:** De Ritz-waarden zijn de eigenwaarden van de **orthogonale projectie** $P_k A P_k$ (geprojecteerd op $\mathcal{K}_k$), d.w.z. $H_k = Q_k^* A Q_k$. Dit is dezelfde constructie als de Galerkin-benadering: eis dat het residu loodrecht staat op de benaderingsruimte.

**Bewijs:** De orthogonale projector op $\mathcal{K}_k$ is $P_k = Q_k Q_k^*$ (want $Q_k$ heeft orthonormale kolommen). Dus:
$$P_k A P_k = Q_k Q_k^* A Q_k Q_k^* = Q_k H_k Q_k^*.$$
De eigenwaarden van $Q_k H_k Q_k^*$ zijn gelijk aan die van $H_k$ (unitaire gelijkenis: $Q_k^*(Q_k H_k Q_k^*)Q_k = H_k$). $\square$

**Voordeel in stap $k$:**
- Eigenwaarden van $A$: $m \times m$ probleem ($m$ groot)
- Ritz-waarden: $k \times k$ probleem ($k \ll m$) → oplosbaar met het volledige QR-algoritme

**Prematuur afbreken ($h_{k+1,k} = 0$):** Als de subdiagonaalcoëfficiënt nul wordt:
$$AQ_k = Q_k H_k \implies \mathcal{K}_k \text{ is A-invariant.}$$
De Ritz-waarden zijn dan **exact** de eigenwaarden van $A$ die beginvector $\mathbf{b}$ kan bereiken.

> ⚠️ **Examentip:** De Ritz-waarden zijn eigenwaarden van $H_k = Q_k^*AQ_k$ — de projectie van $A$ op de Krylov-ruimte. Dit is de sleutelformule.

---

### (d) Variante met gewijzigd Gram-Schmidt (Modified Gram-Schmidt, MGS)?

**Antwoord: Ja, zinvol en in de praktijk gebruikt.**

**Motivering:**

| Eigenschap | Klassiek GS (KGS) | Gewijzigd GS (MGS) |
|---|---|---|
| Wiskundig equivalent? | Ja | Ja |
| Numerieke stabiliteit | Matig | Beter |
| Kost per stap | $O(km)$ | $O(km)$ |

**Waarom MGS beter?** Bij KGS worden alle componenten $r_{ij}$ berekend t.o.v. de **originele** $A\mathbf{q}_j$. Bij MGS worden ze stapsgewijs verwijderd: na het aftrekken van de $q_i$-component, wordt de bijgewerkte $\mathbf{v}_j$ gebruikt voor de volgende orthogonalisatie. Kleine afrondingsfouten propageren minder ver.

**Concreet in Arnoldi:**
```
v = A q_j
voor i = 1 tot j:
    h_{ij} = q_i* v         ← inproduct met bijgewerkte v (MGS)
    v = v - h_{ij} q_i      ← onmiddellijk verwijderen
```

**Nadeel:** Bij ernstig slechte conditie kan ook MGS zijn orthogonaliteit verliezen. Dan is **herorthogonalisatie** nodig.

**Conclusie:** Het Arnoldi-algoritme met MGS (of MGS + herorthogonalisatie) is de standaardimplementatie voor grote eigenwaardesystemen.

> ⚠️ **Examentip:** Het antwoord is "ja, zinvol" met als reden: MGS is numeriek stabieler dan klassiek GS, zonder extra kost.

---

## V9 – Convergentiegedrag van Arnoldi en eigenwaardenlokalisatie

**Vraag:** Bespreek het convergentiegedrag en hoe Arnoldi eigenwaarden lokaliseert, gebruik makend van $(A - \sigma I)Q_k = Q_k(H_k - \sigma I) + \mathbf{q}_{k+1}\mathbf{e}_k^T h_{k+1,k}$ en $\hat{p} = \arg\min_{p \in \mathcal{M}_k}\|p(A)\mathbf{b}\|_2 \iff \hat{p}(z) = \det(zI - H_k)$. Bewijzen niet nodig.

---

### Eigenschap 1 — Verschuivingsinvariantie

De formule $(A - \sigma I)Q_k = Q_k(H_k - \sigma I) + \mathbf{q}_{k+1}\mathbf{e}_k^T h_{k+1,k}$ is een geldige Arnoldi-recursie voor $(A - \sigma I, \mathbf{b})$, met dezelfde $Q_k$ en $H_k - \sigma I$ als Hessenberg-projectie. Omdat $\mathcal{K}_k(A-\sigma I,\mathbf{b}) = \mathcal{K}_k(A,\mathbf{b})$, geeft Arnoldi op $(A-\sigma I)$ Ritz-waarden $\{\theta_i - \sigma\}$.

**Gevolg:** Arnoldi convergeert zowel naar $\lambda_{\min}$ als $\lambda_{\max}$ — in tegenstelling tot de methode van de machten. De stoorterm $h_{k+1,k}$ meet hoe goed $\mathcal{K}_k$ de gezochte eigenruimte benadert: klein $h_{k+1,k}$ $\iff$ Ritz-waarden zijn geconvergeerd.

---

### Eigenschap 2 — Veeltermbenadering en lokalisatie

De karakteristieke veelterm $\hat{p}$ van $H_k$ minimaliseert $\|p(A)\mathbf{b}\|_2$ over alle monische veeltermen van graad $k$. Via de eigenwaardenontbinding $A = X\Lambda X^{-1}$:

$$\min_{p \in \mathcal{M}_k}\left\|X\begin{pmatrix}p(\lambda_1)&&\\&\ddots&\\&&p(\lambda_m)\end{pmatrix}X^{-1}\mathbf{b}\right\|_2.$$

De optimale veelterm plaatst zijn nulpunten zo dat $|p(\lambda_i)|$ klein is voor eigenwaarden die sterk bijdragen aan $\mathbf{b}$.

**Welke eigenwaarden convergeren eerst?** De **geïsoleerde eigenwaarden aan de rand van het spectrum**. Als $\lambda_j$ ver van de rest staat, kan één nulpunt van $p$ dicht bij $\lambda_j$ geplaatst worden zonder andere termen te verstoren $\Rightarrow$ snelle convergentie. Geclusterde eigenwaarden in het midden van het spectrum vereisen veel iteraties.

---

## V10 – Twee methodes voor lage-rangbenadering van een matrix

**Vraag:** Bespreek twee methodes, vergelijk ze, met aandacht voor rekenkost en kwaliteit.

**Probleemstelling:** Gegeven $A \in \mathbb{C}^{m \times n}$ rang $r$, zoek $\min_{\text{rang}(B) \leq k} \|B - A\|$. De verzameling $\mathcal{M}_k = \{B : \text{rang}(B) \leq k\}$ is **geen deelruimte** (niet gesloten onder optelling) $\Rightarrow$ orthogonale projectietheorie niet toepasbaar.

---

### Methode 1 — QR met kolomverwisselingen

Bij elke stap wordt de kolom met grootste norm naar voor gepivoteerd:
$$AP = QR, \qquad |r_{11}| \geq |r_{22}| \geq \cdots$$

**Rang-$k$ benadering:** $A_k^{qr} = Q_k R_k P^T$ ($Q_k$: eerste $k$ kolommen van $Q$, $R_k$: eerste $k$ rijen van $R$).

- **Rekenkost:** $O(mnk)$
- **Kwaliteit:** goed maar **niet gegarandeerd optimaal**

---

### Methode 2 — Getrunceerde SVD

$$A_k^{svo} = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^* = U_k \Sigma_k V_k^*$$

**Gemaakte fout:**
$$\|A - A_k^{svo}\|_2 = \sigma_{k+1}, \qquad \|A - A_k^{svo}\|_F = \sqrt{\sigma_{k+1}^2 + \cdots + \sigma_r^2}$$

- **Rekenkost:** $O(mn^2)$ (duurder)
- **Kwaliteit:** **optimaal** (Stelling Eckart–Young–Mirsky)

---

### Stelling — Beste rang-$k$ benadering (Eckart–Young–Mirsky)

$$\|A_k^{svo} - A\|_2 = \min_{B:\,\text{rang}(B)\leq k} \|B - A\|_2$$

**Bewijs (spectraalnorm):** Te tonen: voor elke $B = XY^*$ met rang $\leq k$ geldt $\|B - A\|_2 \geq \sigma_{k+1}$.

Zoek $\mathbf{w}$ met $\|\mathbf{w}\|_2 = 1$ in $\text{span}\{\mathbf{v}_1,\ldots,\mathbf{v}_{k+1}\} \cap \mathcal{N}(Y^*)$. Zo'n $\mathbf{w}$ bestaat want $(k+1)+(n-k) = n+1 > n$ (dimensies overlappen).

Schrijf $\mathbf{w} = \sum_{i=1}^{k+1} c_i \mathbf{v}_i$ met $\sum c_i^2 = 1$. Dan $B\mathbf{w} = XY^*\mathbf{w} = \mathbf{0}$, dus:
$$\|B - A\|_2 \geq \|A\mathbf{w}\|_2 = \left\|\sum_{i=1}^{k+1} c_i \sigma_i \mathbf{u}_i\right\|_2 = \sqrt{\sum c_i^2 \sigma_i^2} \geq \sigma_{k+1}\underbrace{\sqrt{\sum c_i^2}}_{=1} = \sigma_{k+1}. \quad\square$$

---

### Vergelijking

| | QR met kolomverwisselingen | Getrunceerde SVD |
|---|---|---|
| **Kwaliteit** | Suboptimaal | Optimaal (EYM) |
| **Rekenkost** | $O(mnk)$ | $O(mn^2)$ |
| **Fout** | Niet precies kwantificeerbaar | $\|A - A_k^{svo}\|_2 = \sigma_{k+1}$ |
| **Wanneer** | Snelheid primeert | Beste kwaliteit vereist |

---

### Geval: gigantische ijle matrix ($A \in \mathbb{R}^{10^7 \times 10^7}$)

Voor matrices van deze omvang zijn **beide bovenstaande methodes onhaalbaar**:
- Volledige SVD kost $O(mn^2)$ — bij $m = n = 10^7$ is dit $O(10^{21})$ bewerkingen: volstrekt onhaalbaar.
- Zelfs de matrix in geheugen opslaan vereist $10^{14}$ getallen — niet mogelijk.

**Oplossing: iteratieve Krylov-methoden** die enkel matrix-vector producten $A\mathbf{v}$ vereisen.

**Specifieke aanpak voor symmetrische matrices — Lanczos-algoritme:**

Het Lanczos-algoritme is de variant van Arnoldi voor **symmetrische** matrices. Omdat $A = A^T$, vereenvoudigt de Hessenbergmatrix $H_k$ tot een **tridiagonale** matrix $T_k$. Na $k$ iteraties:
$$T_k = Q_k^T A Q_k, \qquad T_k \text{ tridiagonaal.}$$

De $k$ grootste singuliere waarden van $A$ worden benaderd via de **bidiagonalisatie-aanpak** (Golub-Kahan): men past Lanczos toe op $A^T A$ of $AA^T$ via de relatie $\sigma_i(A) = \sqrt{\lambda_i(A^TA)}$.

**Rekenkost:** slechts $k$ matrix-vector producten $A\mathbf{v}$ (elk $O(nnz)$ voor een ijle matrix) + $O(k^2 n)$ voor de eigenwaarden van $T_k$.

**Vergelijking voor grote ijle matrices:**

| Methode | Rekenkost | Geheugen | Toepassbaar bij $10^7 \times 10^7$? |
|---|---|---|---|
| Volledige SVD | $O(mn^2)$ | $O(mn)$ | Nee |
| QR met kolomverwisseling | $O(mnk)$ | $O(mn)$ | Nee |
| Lanczos/Arnoldi (Krylov) | $O(k \cdot nnz)$ | $O(kn)$ | **Ja** |

> ⚠️ **Examentip:** Zodra Michiels een matrix van orde $10^7$ of meer noemt in de context van SVD of eigenwaarden, is het verwachte antwoord altijd **Krylov-methode (Lanczos voor SVD, Arnoldi voor niet-symmetrisch)**. De sleutelzin: "enkel matrix-vector producten nodig, de matrix hoeft niet expliciet in geheugen."

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
| $L = D - A$ | Graaf-Laplaciaan |
| Fiedler-vector $\mathbf{v}_{N-1}$ | Optimale graafpartitionering |
