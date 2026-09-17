# aesplot — précision le long d'un AES homomorphe

`TestAES` écrit une ligne par opération dans un CSV ; `plot.py` en fait des courbes. Une couleur par
run, `prec_min` (le pire slot sur **tous** les blocs) en ordonnée, et en abscisse le temps de
pipeline cumulé — `--x seq` remet l'indice de l'opération.

Toutes les commandes de ce README se lancent **depuis la racine du dépôt**.

## Produire les traces

Le flag `-csv` de `TestAES` **ajoute** au fichier : pointe tous tes runs sur le même et ils
s'empilent.

```sh
go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -csv runs/aes.csv
go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -cleanextract -csv runs/aes2.csv ; go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -cleanextract -xor sq -csv runs/aes2.csv
```

`-csv` est relatif au **répertoire du package**, pas au tien : `go test` y place le binaire de test,
donc `runs/aes.csv` écrit dans `nathanfau/transciphering/runs/`. `plot.py`, lui, résout depuis ton
cwd — d'où les chemins complets dans tout ce qui suit.

Les flags qui changent la configuration : `-subbytes`, `-xor`, `-clean`, `-place`, `-cleanextract`,
`-seed`, et pour la chaîne `-logn` (défaut 11), `-stc` (primes de SlotsToCoeffs, défaut `2x30`) et
`-logqi` (taille des primes q_i de la chaîne et de l'échelle, défaut 38 ; q0 suit, P et `-stc` non).
La chaîne se relit dans les colonnes `logn`, `logscale`, `q_sizes` et `s2c_levels`.

Chaque run porte son horodatage de démarrage (`run_ts`) et **toutes** ses colonnes de configuration,
répétées sur chaque ligne : le fichier se suffit à lui-même et deux runs se concatènent sans
histoire d'en-tête. Un run coûte une dizaine de minutes ; le fichier est validé au démarrage, donc
un en-tête incompatible échoue tout de suite et pas à la fin.

## Tracer

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv
```

C'est tout ce qu'il faut pour le cas courant : un PNG par valeur de `logn`, écrit à côté du CSV —
`nathanfau/transciphering/runs/aes-logn11.png`, `…-logn12.png`… Ajoute un run à logN 13 et le
prochain appel sortira `aes-logn13.png` sans que tu touches à quoi que ce soit.

Le découpage par défaut est `--split logn` parce que deux degrés d'anneau ne se comparent pas sur un
même axe : le bruit de décodage CKKS croît en `σ√N`, donc l'écart entre deux logN n'est pas une
propriété du circuit.

| flag | effet |
|---|---|
| `--x time\|seq` | abscisse : le temps de pipeline cumulé (défaut) ou le rang de l'opération |
| `--where COL=VAL` | ne garder que ces runs. Répétable |
| `--split COL` | un PNG par valeur de `COL`. Défaut `logn`, `none` pour tout réunir |
| `--avg` | ajoute `prec_avg` en pointillés (par défaut, seul le pire slot) |
| `--all` | garder plusieurs passes d'une même configuration |
| `--out` | un nom de fichier précis (incompatible avec `--split`) |
| `--title` | remplace le titre, filtres compris |
| `--dpi` | 160 par défaut |

### Exemples

**Un seul degré d'anneau.** Le suffixe n'est pas doublé : un `--where` qui épingle déjà la colonne
de découpage la désactive. Sort `runs/aes-logn11.png`.

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where logn=11
```

**Isoler un circuit XOR**, aux deux degrés. Les filtres se composent avec le découpage, donc deux
PNG : `runs/aes-xorsq-logn11.png` et `…-logn12.png`.

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where xor=sq
```

**Croiser deux filtres.** Répète `--where` :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv \
    --where logn=11 --where xor=nosq
```

**Filtrer sur l'extraction sans taper son libellé.** `extract` vaut `half (demi-spectre, k lv)`,
pénible à citer dans un shell ; `extractlv` vaut 4 ou 5 et dit la même chose :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where extractlv=5
```

**Une seule opération, tour par tour.** `--where step=...` ne garde que ces lignes : dix points par
run au lieu de cinquante, et on lit directement ce que le refresh rend à chaque tour.

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where step=Refresh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where step=Cleaning
```

**Découper sur autre chose que `logn`.** Un PNG par extraction, ou par circuit XOR :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --split extract
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --split xor
```

Attention : `--split` remplace le découpage par logN, il ne s'y ajoute pas. Ces deux commandes
mélangent donc les degrés d'anneau sur un même axe — c'est justement ce que le défaut évite.

**Ajouter la moyenne.** Trait plein pour le pire slot, pointillés pour la moyenne, même couleur :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where logn=11 --avg
```

**Une figure ponctuelle**, nom et titre choisis :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where logn=11 \
    --out /tmp/figure.png --title "Effet du cleaning fusionne" --dpi 220
```

**Voir la dispersion entre deux passes** d'une même configuration, que le script écarte par défaut :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv --where logn=11 --all
```

**Tout dans un seul graphique**, si tu sais pourquoi :

```sh
python3 nathanfau/transciphering/aesplot/plot.py --csv nathanfau/transciphering/runs/aes.csv \
    --split none --out /tmp/tout.png
```

Le nom du PNG découle du CSV et des filtres, les caractères non alphanumériques retirés : `--where
xor=sq` donne `nathanfau/transciphering/runs/aes-xorsq-logn11.png`. Les colonnes filtrables sont
toutes celles du CSV, listées plus bas.

Le script imprime aussi le tableau des moyennes par opération. C'est la même information en texte,
lisible quand la couleur ne l'est pas.

## Trois choses à savoir pour lire les tracés

**Un run par configuration.** Si tu relances le même test, la nouvelle passe **remplace** l'ancienne
dans le tracé plutôt que d'ajouter une courbe superposée. Le script dit à l'écran ce qu'il écarte.
`--all` pour tout voir. La configuration est définie **sans la seed** : sans `-seed`, elle est tirée
de l'horloge et changerait à chaque run sans que rien n'ait bougé.

**Une courbe collée à 0 est morte, pas mauvaise.** Sous 0 bit l'erreur dépasse 1 : le bit n'existe
plus, et de combien il n'existe plus ne veut rien dire. Un run qui plonge à −110 écraserait tous les
autres, donc le tracé plafonne à 0 — **au tracé seulement**, le CSV garde la valeur réelle. Une note
sous l'axe le rappelle.

**Une précision ne prouve rien si l'oracle a décroché.** La colonne `blocks_wrong` dit combien de
blocs ne correspondent plus à l'AES en clair. Au-delà du premier point où elle est non nulle, la
courbe mesure une dérive, pas une précision :

```sh
python3 -c "import pandas as pd; d=pd.read_csv('nathanfau/transciphering/runs/aes.csv'); print(d[d.blocks_wrong>0][['run_ts','round','step','blocks_wrong']].to_string())"
```

## Colonnes du CSV

Tout ce qui suit est **répété sur chaque ligne**, pour qu'aucune ne soit ambiguë sur ce qui l'a
produite.

*Le run.* `run_ts` (démarrage), `host`, `gomaxprocs`, `go_version`, `git_commit` (avec un `+` final
si l'arbre était sale). Un fichier qui s'accumule sur des semaines n'est pas relisible sans eux —
et les temps, eux, dépendent de la machine.

*Ce qui a été lancé.* `logn`, `k`, `slots`, `blocks`, `seed`, `rounds`, `subbytes`, `xor`, `clean`,
`cleandepth`, `place`, `extract`, `extractlv`.

*La chaîne.* `logscale`, `primes` (leur nombre), `logq0`, `logq`, `logp`, `logqp`, `dnum`,
`digit_bits` (le plus large paquet de `#P` primes consécutives), `marge` (`logp − digit_bits`,
ce qui reste au-dessus du bruit de key-switch), `q_sizes` et `p_sizes` (les tailles, arrondies),
et `chain_id`. Ce dernier hache les primes **elles-mêmes** : changer la taille d'une seule prime
fait servir au générateur une autre prime à tous les niveaux en dessous, et ça vaut des bits de
précision. Deux lignes partagent un `chain_id` si et seulement si elles ont tourné sur la même
chaîne.

*Le secret et le bootstrap.* `h` (poids de Hamming du secret), `h_tilde` (celui du secret
éphémère), `s2c_levels` et `c2s_levels` (le découpage des deux DFT), `mod1_type`, `mod1_deg`,
`mod1_k`, `logmsgratio`.

*Les clefs.* `galois_keys` (leur nombre, que le découpage du SlotsToCoeffs fait varier d'un facteur
3), `key_gb`, `keygen_ms`, `dft_ms` (l'encodage des matrices DFT, qui n'est pas de la génération de
clefs et coûte plus qu'elle).

*Où le pipeline se pose.* `refresh_lv`, `ark_lv`.

Position et coût : `seq` (l'ordre d'exécution), `round`, `step`, `ms` (la durée de l'opération) et
`elapsed_ms` (leur somme courante). Ces deux abscisses ne disent pas la même chose : `seq` compte
des opérations, et deux configurations qui n'en font pas le même nombre — `XorClean` en fusionne
deux — ne s'y comparent pas. `elapsed_ms` est le coût du **pipeline seul** : le chronomètre du run
porte en plus les traces d'oracle, qui ne sont pas du calcul homomorphe.

Résultat : `level`, `prec_avg`, `prec_min`, `worst_err`, `slots_pooled`, `bit_err`, `blocks_wrong`.

`prec_avg` et `prec_min` sont poolées sur **les 64 chiffrés et tous leurs slots** — pas par octet.
`prec_min` est le `−log2` de la plus grande distance du lot, donc la marge qui reste avant qu'un bit
bascule ; `slots_pooled` dit sur combien de valeurs.

## Couleurs

Palette catégorielle validée, dans un **ordre fixe** : le run *n* prend toujours la couleur *n*, et
elle n'est jamais recyclée. Au-delà de 8 runs le script s'arrête plutôt que de générer une 9ᵉ teinte
qui ne serait plus distinguable — filtre le fichier. Chaque courbe porte aussi son étiquette en bout
de tracé, pour que l'identité ne repose pas sur la seule couleur.
