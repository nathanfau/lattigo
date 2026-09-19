#!/usr/bin/env python3
"""Précision le long d'un run d'AES homomorphe, une courbe par run.

Lit le CSV que `go test -run '^TestAES$' -csv <fichier>` écrit (une ligne par
opération, plusieurs runs empilés dans le même fichier) et trace, sur UN seul
graphique :

  - en abscisse  : le temps de pipeline cumulé (colonne `elapsed_ms`), ou le rang de
                   l'opération avec --x seq
  - en ordonnée  : la précision en bits
  - une courbe   : `prec_min`, le pire slot de tous les blocs (--avg ajoute la moyenne
                   en pointillés)
  - une couleur  : un run

    go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -csv runs/2026-09-07_aes.csv
    go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -cleanextract -csv runs/2026-09-07_aes.csv
    python3 nathanfau/transciphering/aesplot/plot.py --csv runs/2026-09-07_aes.csv
"""

import argparse
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

# Palette catégorielle de référence, dans son ORDRE FIXE : une couleur par run, jamais recyclée.
# Validée pour la liste de paires adjacentes (le cas des courbes) en mode clair — ΔE CVD 9.1,
# vision normale 19.6. Ne pas réordonner ni intercaler.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e3e2df"

# L'axe des ordonnées ne descend pas plus bas, quoi que fassent les données : voir set_ylim.
FLOOR = -5.0

# Ce qui définit une CONFIGURATION. `seed` en est exclue : sans -seed elle est tirée de l'horloge,
# donc elle diffère à chaque run sans que rien de la configuration n'ait changé. Le label de légende
# est construit sur les colonnes qui DIFFÈRENT entre les configurations retenues.
# `logscale` (la taille des q_i, flag -logqi) et `s2c_levels` (la forme de -stc) passent AVANT
# `primes` et `chain_id` : quand ils changent, c'est eux que l'etiquette courte doit nommer.
# `chain_id` identifie les primes elles-memes : deux chaines de memes TAILLES mais dont une prime
# a change plus bas ne portent pas les memes primes, et ca vaut des bits de precision.
CONFIG = ["logn", "k", "blocks", "rounds", "subbytes", "xor", "clean",
          "cleandepth", "clean_scale", "refresh_scale", "place", "extract", "extractlv", "logscale", "s2c_levels", "primes", "chain_id"]


def pick_runs(df, keep_all):
    """Un run par configuration, le plus récent. Le même test relancé n'ajoute pas une courbe :
    il en remplace une, sinon deux tracés identiques mangent deux couleurs pour rien."""
    per_run = df.groupby("run_ts", sort=True)[CONFIG].first()
    if keep_all:
        return list(per_run.index), []
    kept, dropped = {}, []
    for ts, row in per_run.iterrows():
        key = tuple(row)
        if key in kept:
            dropped.append(kept[key])  # les precedents cedent la place au plus recent
        kept[key] = ts
    return sorted(kept.values()), sorted(dropped)


def run_labels(df):
    """Deux libellés par run : le long pour la légende, le court pour le bout de courbe."""
    per_run = df.groupby("run_ts", sort=True)[CONFIG].first()
    varying = [c for c in CONFIG if c in per_run and per_run[c].nunique() > 1]
    # Le nom du fichier ne se lit que s'il separe des runs que les colonnes de configuration ne
    # separent pas deja : sinon il allonge la legende sans rien apprendre.
    if "source" in varying:
        rest = [c for c in varying if c != "source"]
        if rest and len(per_run[rest].drop_duplicates()) == len(per_run[varying].drop_duplicates()):
            varying = rest
    # L'etiquette courte prend le plus petit prefixe de `varying` qui suffise a distinguer les
    # runs : avec un seul axe qui change elle vaut "clean", avec deux elle vaut "11/clean".
    def tag(row, cols):
        return "/".join(str(row[c]).split()[0] for c in cols)

    # On n'ajoute une colonne que si elle separe VRAIMENT : `blocks` est determine par `logn`,
    # `extractlv` et `primes` par `extract`, et les repeter n'apprend rien a la lecture.
    cols, seen = [], 1
    for c in varying:
        n = per_run.apply(lambda r: tag(r, cols + [c]), axis=1).nunique()
        if n > seen:
            cols.append(c)
            seen = n
        if seen == len(per_run):
            break

    full, short = {}, {}
    for i, (ts, row) in enumerate(per_run.iterrows()):
        if varying:
            full[ts] = ", ".join(f"{c}={row[c]}" for c in varying)
            short[ts] = tag(row, cols)
        else:
            full[ts] = str(ts)  # runs identiques : seul l'horodatage les separe
            short[ts] = f"#{i + 1}"
    return full, short


def add_x(df, mode):
    """Pose la colonne d'abscisse `_x`, le début de l'opération `_x0`, et rend le libellé de l'axe.

    En temps, l'abscisse est le cumul des durées d'opération, donc le coût du PIPELINE seul : le
    chronomètre du run porte en plus les traces d'oracle, qui ne sont pas du calcul homomorphe.
    C'est la seule abscisse sur laquelle deux configurations qui ne font pas le même nombre
    d'opérations se comparent."""
    if mode == "seq":
        df["_x"] = df.seq.astype(float)
        df["_x0"] = df._x - 0.5  # le separateur de tour tombe ENTRE deux operations
        return "n — n-ème opération du run"
    scale, unit = (1 / 60_000, "min") if df.elapsed_ms.max() > 900_000 else (1 / 1000, "s")
    df["_x"] = df.elapsed_ms * scale
    df["_x0"] = df._x - df.ms * scale  # le debut de l'operation, i.e. la fin de la precedente
    return f"temps de pipeline cumulé ({unit})"


def round_marks(df, gap=0.04):
    """L'abscisse où chaque tour commence. Le tour 0 n'a qu'une opération, donc son séparateur
    tomberait sur celui du tour 1 : on ne garde que les marques assez espacées pour être lues,
    le seuil étant une fraction de la largeur du graphique."""
    one = df[df.run_ts == df.run_ts.iloc[0]].sort_values("seq")
    starts, prev_round = [], None
    for _, r in one.iterrows():
        if r["round"] != prev_round:
            starts.append((float(r["_x0"]), int(r["round"])))
            prev_round = r["round"]
    # De deux marques trop proches on garde la SECONDE : le tour 0 n'a qu'une operation, et c'est
    # le debut du tour 1 qui interesse.
    span = float(one._x.max() - one._x.min()) or 1.0
    return [m for i, m in enumerate(starts)
            if i + 1 == len(starts) or starts[i + 1][0] - m[0] >= gap * span]


def draw(df, args, where):
    """Un graphique pour le sous-ensemble deja filtre. `where` sert au titre et au nom du PNG."""
    runs, dropped = pick_runs(df, args.all)
    if dropped:
        print(f"{len(dropped)} run(s) ecarte(s), meme configuration qu'un run plus recent "
              f"(--all pour tout tracer) : {', '.join(dropped)}")
    df = df[df.run_ts.isin(runs)]
    if len(runs) > len(SERIES):
        sys.exit(f"{len(runs)} runs pour {len(SERIES)} couleurs : filtre le fichier, "
                 "une 9e couleur ne serait plus distinguable")

    full, short = run_labels(df)
    if args.label:
        # La legende vient d'une colonne choisie : utile quand les runs different sur trop de
        # colonnes de configuration pour qu'une liste reste lisible.
        per_run = df.groupby("run_ts")[args.label].first().astype(str)
        full, short = per_run.to_dict(), per_run.to_dict()
    fig, ax = plt.subplots(figsize=(11, 5.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    # Séparateurs de tours, derrière les données et volontairement discrets : c'est par tour qu'on
    # lit une dégradation. En temps ils ne valent que pour UN run : deux configurations n'arrivent
    # pas au tour 5 au même instant, et une grille tiree du premier run mentirait sur les autres.
    if args.x == "seq" or len(runs) == 1:
        for x0, rnd in round_marks(df):
            ax.axvline(x0, color=GRID, lw=1, zorder=0)
            ax.annotate(f"T{rnd}", xy=(x0, 1.005), xycoords=("data", "axes fraction"),
                        color=INK_MUTED, fontsize=7.5, ha="left", va="bottom")

    ends = []
    for i, ts in enumerate(runs):
        color = SERIES[i]
        d = df[df.run_ts == ts].sort_values("seq")
        ax.plot(d._x, d.prec_min, color=color, lw=2, zorder=3)
        if args.avg:
            ax.plot(d._x, d.prec_avg, color=color, lw=2, ls=(0, (2, 2)), zorder=3)

        # Étiquette directe en bout de courbe : l'identité ne repose jamais sur la seule couleur,
        # et trois des huit teintes passent sous 3:1 sur fond clair.
        last = d.dropna(subset=["prec_min"]).iloc[-1]
        ends.append((float(last.prec_min), float(last._x), color, short[ts]))

    ax.set_xlabel(args.xlabel, color=INK_MUTED, fontsize=9.5)
    ax.set_ylabel("précision du pire slot (bits)" if not args.avg else "précision (bits)",
                  color=INK_MUTED, fontsize=9.5)
    title = args.title or "Précision le long d'un AES-128 homomorphe"
    if where and not args.title:
        title += "  —  " + ", ".join(where)
    ax.set_title(title, color=INK, fontsize=12, pad=22, loc="left")

    ax.grid(axis="y", color=GRID, lw=1, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=8.5)

    # Les valeurs sont tracées telles quelles, mais l'axe ne descend pas sous FLOOR : sous 0 bit
    # l'erreur dépasse 1, le bit n'existe plus, et de COMBIEN il n'existe plus ne veut rien dire.
    # Un run qui plonge à -110 écraserait tous les autres dans une bande illisible.
    lo = float(min(df.prec_min.min(), df.prec_avg.min() if args.avg else 0))
    ax.set_ylim(bottom=max(FLOOR, lo - 0.5))
    if lo < FLOOR:
        ax.annotate(f"une courbe qui sort par le bas est passée sous {FLOOR:.0f} bits : "
                    "l'erreur y dépasse 1 depuis longtemps, le bit n'existe plus",
                    xy=(0, -0.13), xycoords="axes fraction",
                    color=INK_MUTED, fontsize=8, ha="left", va="top")

    # Étiquettes de bout de courbe, écartées quand deux fins se touchent : plusieurs runs morts
    # finissent tous sur 0, et leurs étiquettes se superposeraient exactement.
    y0, y1 = ax.get_ylim()
    gap = 0.05 * (y1 - y0)
    placed = None
    for y, x, color, tag in sorted(ends):
        ax.plot([x], [y], marker="o", ms=5, color=color, zorder=4)
        y = max(y0, min(y1, y))  # une courbe finie hors cadre garde son etiquette sur le bord
        if placed is not None and y - placed < gap:
            y = placed + gap
        placed = y
        ax.annotate(f"  {tag}", xy=(x, y), color=INK_MUTED, fontsize=9,
                    va="center", ha="left", zorder=4, annotation_clip=False)

    # Place pour les etiquettes de bout de courbe, mesuree sur la plus longue et rapportee a la
    # largeur des donnees : en secondes comme en rangs d'operation. Les graduations, elles,
    # s'arretent aux donnees : une marge n'est pas une plage de valeurs.
    lo, hi = float(df._x.min()), float(df._x.max())
    span = (hi - lo) or 1.0
    pad = span * (0.012 * max(len(s) for s in short.values()) + 0.02)
    ax.set_xlim(lo - 0.02 * span, hi + pad)
    ax.set_xticks([t for t in ax.get_xticks() if lo <= t <= hi])

    # Deux légendes, SOUS les axes : la couleur dit quel run, le style de trait quelle statistique.
    # Hors du cadre, elles ne peuvent pas recouvrir une courbe, quelle que soit l'allure des données.
    # La legende des couleurs suffit quand une seule statistique est tracee : le style de trait
    # n'encode plus rien, donc sa legende n'aurait rien a dire.
    rows = len(runs) + (1 if args.avg else 0)
    fig.subplots_adjust(bottom=0.10 + 0.045 * rows)
    if args.avg:
        style_legend = fig.legend(
            handles=[Line2D([], [], color=INK_MUTED, lw=2, label="prec_min — pire slot, trait plein"),
                     Line2D([], [], color=INK_MUTED, lw=2, ls=(0, (2, 2)), label="prec_avg — moyenne, pointillés")],
            loc="lower left", bbox_to_anchor=(0.08, 0.005), frameon=False, fontsize=8.5,
            labelcolor=INK_MUTED, ncol=2)
        fig.add_artist(style_legend)
    fig.legend(handles=[Line2D([], [], color=SERIES[i], lw=2, label=full[ts])
                        for i, ts in enumerate(runs)],
               loc="lower left", bbox_to_anchor=(0.08, 0.005 + (0.045 if args.avg else 0)),
               frameon=False, fontsize=8.5, labelcolor=INK_MUTED, ncol=1)

    out = args.out
    if not out:
        # Les valeurs des colonnes portent des espaces et des parentheses -- "sq ((x-y)^2)" --
        # qu'un nom de fichier ne doit pas heriter.
        suffix = "".join("-" + re.sub(r"[^0-9A-Za-z]+", "", w) for w in where)
        out = args.stem + suffix + ".png"
    fig.savefig(out, dpi=args.dpi, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"{len(runs)} run(s), {len(df)} lignes -> {out}")

    # La même chose en texte : une courbe sous 3:1 de contraste doit rester lisible autrement.
    piv = df.pivot_table(index="step", columns="run_ts", values=["prec_min", "prec_avg"], aggfunc="mean")
    print(piv.round(2).to_string())
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", action="append", default=None, metavar="FICHIER",
                    help="repetable : plusieurs CSV sont reunis, et une colonne `source` (le nom "
                         "du fichier) distingue leurs runs. Defaut runs/2026-09-07_aes.csv")
    ap.add_argument("--out", default=None,
                    help="PNG de sortie (defaut : le nom du CSV, plus les filtres). Incompatible avec --split")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--where", action="append", default=[], metavar="COL=VAL",
                    help="ne garder que les runs dont COL vaut VAL, repetable (ex: --where logn=11)")
    ap.add_argument("--split", default="logn", metavar="COL",
                    help="un PNG par valeur de COL trouvee dans le CSV. Defaut logn : deux degres "
                         "d'anneau ne se comparent pas sur un meme axe. 'none' pour tout reunir")
    ap.add_argument("--x", choices=("time", "seq"), default="time",
                    help="abscisse : 'time' le temps de pipeline cumule (defaut), 'seq' le rang "
                         "de l'operation. Deux configurations qui ne font pas le meme nombre "
                         "d'operations ne se comparent qu'en temps")
    ap.add_argument("--title", default=None, help="titre du graphique")
    ap.add_argument("--avg", action="store_true",
                    help="tracer aussi prec_avg, en pointilles (par defaut seul le pire slot)")
    ap.add_argument("--label", default=None, metavar="COL",
                    help="legende et etiquettes de bout de courbe tirees de la colonne COL, au lieu "
                         "de la liste des colonnes de configuration qui different")
    ap.add_argument("--all", action="store_true",
                    help="tracer TOUS les runs, y compris plusieurs passes d'une meme configuration")
    args = ap.parse_args()

    if args.split not in ("", "none") and args.out:
        sys.exit("--split ecrit plusieurs PNG, il ne peut pas partager un seul --out")

    args.csv = args.csv or ["runs/2026-09-07_aes.csv"]
    args.stem = args.csv[0].rsplit(".", 1)[0]

    # Plusieurs fichiers se reunissent sur l'union de leurs colonnes : les jeux de colonnes ont
    # change dans le temps, et un CSV ancien n'a pas a etre exclu pour autant.
    parts = []
    for path in args.csv:
        d = pd.read_csv(path)
        if d.empty:
            sys.exit(f"{path} : aucune ligne")
        if len(args.csv) > 1:
            d["source"] = path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    # En TETE de CONFIG : l'etiquette courte prend les premieres colonnes qui separent, et le nom
    # du fichier se lit ou un chain_id en hexadecimal ne se lit pas.
    if len(args.csv) > 1 and "source" not in CONFIG:
        CONFIG.insert(0, "source")
    # Les jeux de colonnes ont change dans le temps. Une colonne de CONFIG que le fichier n'a pas
    # du tout fait echouer le groupby de pick_runs ; presente dans un fichier et pas dans l'autre,
    # elle vaut NaN, et un groupby laisse alors tomber ces lignes sans un mot. Dans les deux cas
    # elle vaut "?", ce qui ne separe rien et ne casse rien.
    # Exception : un run sans `clean_scale` date d'avant le flag -cleanfixed, donc d'un Cleaning
    # qui gardait l'echelle de son entree -- c'est un fait, pas une inconnue.
    if "clean_scale" in df.columns:
        df["clean_scale"] = df.clean_scale.fillna("input")
    else:
        df["clean_scale"] = "input"
    # Idem pour `refresh_scale` : avant -refreshcanon, le refresh reetiquetait toujours.
    if "refresh_scale" in df.columns:
        df["refresh_scale"] = df.refresh_scale.fillna("relabel")
    else:
        df["refresh_scale"] = "relabel"
    for c in CONFIG:
        df[c] = df[c].fillna("?") if c in df.columns else "?"

    # Les CSV ecrits avant la colonne `elapsed_ms` la retrouvent exactement : elle EST la somme
    # courante des durees d'operation du run.
    if "elapsed_ms" not in df.columns:
        df["elapsed_ms"] = float("nan")
    df["elapsed_ms"] = df.elapsed_ms.fillna(df.groupby("run_ts").ms.cumsum())
    args.xlabel = add_x(df, args.x)

    for w in args.where:
        col, _, val = w.partition("=")
        if col not in df.columns:
            sys.exit(f"--where {w} : pas de colonne {col!r}")
        df = df[df[col].astype(str) == val]
        if df.empty:
            sys.exit(f"--where {w} : aucun run ne correspond")

    # Un --where qui epingle deja la colonne de decoupage ne laisse qu'une valeur : decouper
    # dessus n'ajouterait qu'un suffixe en double au nom du fichier.
    pinned = {w.partition("=")[0] for w in args.where}
    if args.split in ("", "none") or args.split in pinned:
        draw(df, args, args.where)
        return

    if args.split not in df.columns:
        sys.exit(f"--split {args.split} : pas de colonne {args.split!r}")
    for v in sorted(df[args.split].astype(str).unique()):
        draw(df[df[args.split].astype(str) == v], args, args.where + [f"{args.split}={v}"])


if __name__ == "__main__":
    main()
