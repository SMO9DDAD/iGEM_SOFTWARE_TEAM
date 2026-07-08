# SCENTINEL CODE: Canvis V2

![Example Image](1.png)

## Resum

La V1 puntuava cada OBP amb un sol joc de pesos fix. La V2 afegeix **dues maneres noves** de tractar la incertesa dels pesos, tres taules noves per llegir el resultat, dues correccions importants als criteris i un mode nou que genera un informe PDF complet. Aquest document resumeix tot el que ha canviat.

---

## Recordatori dels 4 criteris i el seu ordre

El programa avalua cada OBP amb afinitat (s1), selectivitat (s2), estabilitat (s4) i promiscuïtat (s5). L'ordre d'importància per defecte és:

> **afinitat (0,45) > selectivitat (0,25) > estabilitat (0,15) = promiscuïtat (0,15)**

---

## Dos nous modes!

La idea de fons: sovint sabem quin criteri és més important que un altre, però no sabem exactament **quant**. Els dos modes nous serveixen per a aquesta situació, cadascun amb un grau de certesa diferent.

### Mode 1: Pesos deterministes (recordatori)

El de sempre: poses un joc de pesos i obtens un únic rànquing. Ràpid i transparent. Serveix de base i és el punt de partida dels altres dos modes.

### Mode 2: Pertorbació (c%)

**Què fa:** parteix dels teus pesos i els deixa "ballar" una mica al voltant. Repeteix el càlcul moltes vegades (per defecte 200.000; ajustable) amb pesos lleugerament diferents cada cop i compta quantes vegades acaba guanyant cada OBP.

**Per a què serveix:** quan estàs força segur dels pesos — vols, per exemple, que la selectivitat valgui ~el doble que l'estabilitat — però no al 100%. Aquest mode et diu si el guanyador aguanta quan mous una mica els pesos o si només era qüestió de sort.

**La matemàtica, simplificada:**

- Imagina els pesos ordenats de gran a petit com una escala amb graons.
- El programa no mou els pesos directament: mou els **graons** (els salts entre un criteri i el següent) i la **base** (el valor comú dels dos últims criteris).
- Cada graó i la base poden encongir-se o estirar-se com a màxim un c% del seu valor.
- `c = 0` → pesos fixos (idèntic al Mode 1). `c = 100` → el graó més petit es pot tancar del tot (els dos criteris més propers queden quasi empatats).
- Després es torna a muntar l'escala i es normalitza perquè els pesos sumin 1.

Així mai es trenca l'ordre (un criteri menys important no pot superar-ne un de més important) i l'estabilitat i la promiscuïtat sempre queden iguals. A la pràctica: `c = 10` és molt concentrat (pesos quasi fixos), `c = 50` és el recomanat (mig salt) i `c = 100` és el màxim.

### Mode 3: SMAA pur

**Què fa:** no fas servir cap valor de pes concret. Només dius l'ordre d'importància i el programa prova **tots** els pesos possibles que respecten aquest ordre, comptant amb quina freqüència guanya cada OBP.

**Per a què serveix:** quan saps que l'afinitat importa més que la selectivitat, però no tens ni idea de si és 2× o 5× més important. Deixes que el mètode exploti tot l'espai i et digui quin candidat surt bé "gairebé passi el que passi".

**La matemàtica, simplificada:**

- Es generen pesos aleatoris **uniformes** sobre la regió permesa: w₁ > w₂ > … > (estabilitat = promiscuïtat).
- Tècnicament es fa amb una distribució de Dirichlet sobre els graons (els salts entre criteris) més la base, de manera que cada combinació vàlida té exactament la mateixa probabilitat de sortir.
- Resultat: explores honestament tot l'univers de pesos compatible amb l'ordre, sense afegir-hi cap opinió sobre els valors.

### Mode 2 vs Mode 3 d'un cop d'ull

|  | **Mode 2 — Pertorbació** | **Mode 3 — SMAA pur** |
|---|---|---|
| Què saps | Els pesos aproximats (vols donar-los un valor concret). | Només l'ordre d'importància, no els valors. |
| Què explora | Un entorn al voltant dels TEUS pesos, amb marge c%. | TOT l'espai de pesos compatible amb l'ordre. |
| Paràmetre clau | c% (dispersió). | Cap: només l'ordre dels criteris. |
| Perfil de pesos | No es calcula. | Sí: en quin "món de pesos" guanya cada OBP. |

---

## Tres maneres noves de veure el resultat!

Les tres taules surten en els tres modes. Distingeixen els candidats segons si tenim o no totes les dades dels interferents (i, per tant, si la selectivitat és de fiar).

### Taula 1: Selectivitat fiable

Només hi surten els OBPs que tenen **tots** els interferents mesurats. Per a aquests, la selectivitat (s2) és de fiar i es fa la comparació honesta amb els 4 criteris (s1+s2+s4+s5).

> **Cas especial:** si cap OBP té tots els interferents complets, el programa no s'inventa res. T'avisa i et demana tornar a posar els pesos **sense** selectivitat (es bloqueja a 0 i el seu pes es reparteix entre la resta de criteris). Després salta directament a la Taula 2.

### Taula 2: Prometedors (ignorant la selectivitat)

Hi surten els OBPs als quals falta ≥ 1 interferent. Com que el que falta podria amagar el mínim real, s'elimina el criteri 2 (selectivitat) i es puntua només amb la resta (s1+s4+s5, renormalitzats). Són candidats bons, però els queda confirmar la selectivitat.

### Taula 3: Llista de feina (Boltz / Tanimoto)

Agafa les parelles OBP/interferent que falten i les ordena segons la prioritat de la Taula 2 (primer les dels OBPs més prometedors). És la teva **to-do list**: quines afinitats has de calcular abans —per Boltz o Tanimoto— perquè val més la pena.

### Sortida extra dels Modes 2 i 3

A més de les taules, els modes probabilístics mostren, per a cada OBP, la **probabilitat de ser el 1r, el 2n, el 3r…** (p(rang 1), p(rang 2)…), amb una barra visual.

Només el **Mode 3** afegeix el **perfil de pesos típics**: en quin "món de pesos" (quina combinació concreta) guanya cada candidat. És la resposta a "quina és la situació ideal que fa que aquest OBP surti primer?".

---

## Correccions d'errors

### Afinitat (s1) — nova normalització

**Abans:** s'escalava amb el mín/màx **dins de cada compost**. Això separava molt els candidats dins d'un mateix VOC, però perdia el sentit absolut: un binder mediocre de 18,79 μM treia s1 = 1 només per ser el millor d'aquell VOC. Resultat: impossible comparar entre VOCs diferents.

**Ara:** àncores fixes de la literatura → sostre 1 μM / terra 40 μM (Cui et al.):

- **40 μM** = llindar de literatura per a unió feble → s1 = 0.
- **1 μM** = sostre pràctic per a biosensors: per sota d'1 μM ja és excel·lent i millorar més aporta poc → s1 = 1.
- **Escala independent del dataset:** afegir o treure una proteïna no canvia cap valor de s1.
- Bona discriminació en el rang crític 1–40 μM (~2 ordres de magnitud, ben resolts).

### Selectivitat (s2) — τ personalitzable

Ara pots triar el valor de τ (tau). La fórmula és:

```
s2 = min( 1,  Ki_interferent / (Ki_diana × τ) )
```

| Valor | Descripció |
|---|---|
| `τ = 1` | Lax — l'interferent només ha de lligar una mica menys que la diana. |
| `τ = 10` | **Recomanat** — l'interferent ha de lligar 10× més fluix que la diana. |
| `τ = 50–100` | Exigent — per a dianes i interferents molt semblants. |

Rang vàlid: [1, 100]. Per defecte, **τ = 10**.

---

## NOU MODE! *(en procés…)*

El fitxer executa els 3 mètodes de cop i genera un **informe PDF descarregable** a cada execució. Inclou:

- Portada amb tota la configuració (VOC diana, interferents, pesos, τ, iteracions, c%).
- Mode 1: gràfic de scores amb la descomposició dels subscores i taula detallada.
- Taules 1 / 2 i Taula 3 (llista Boltz), si hi ha dades incompletes.
- Mode 2 (Pertorbació) i Mode 3 (SMAA pur, amb els perfils de pesos).
- Comparació final dels 3 mètodes, un al costat de l'altre.

Així, en lloc de llegir-ho tot per consola, tens un document gràfic per guardar, comentar amb el PI i adjuntar a la documentació del projecte.
