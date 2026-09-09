#!/usr/bin/env python3
"""Précision le long d'un run d'AES homomorphe, une courbe par run.

Lit le CSV que `go test -run '^TestAES$' -csv <fichier>` écrit (une ligne par
opération, plusieurs runs empilés dans le même fichier) et trace, sur UN seul
graphique :

  - en abscisse  : n, la n-ème opération du run (colonne `seq`)
  - en ordonnée  : la précision en bits
  - une courbe   : `prec_min`, le pire slot de tous les blocs (--avg ajoute la moyenne
                   en pointillés)
  - une couleur  : un run

    go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -csv runs/aes.csv
    go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -cleanextract -csv runs/aes.csv
    python3 nathanfau/transciphering/aesplot/plot.py --csv runs/aes.csv
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
CONFIG = ["logn", "k", "blocks", "rounds", "subbytes", "xor", "clean",
          "cleandepth", "place", "extract", "extractlv", "primes"]


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


def round_marks(df, gap=2):
    """Le `seq` où chaque tour commence. Le tour 0 n'a qu'une opération, donc son séparateur
    tomberait sur celui du tour 1 : on ne garde que les marques assez espacées pour être lues."""
    one = df[df.run_ts == df.run_ts.iloc[0]].sort_values("seq")
    starts, prev_round = [], None
    for _, r in one.iterrows():
        if r["round"] != prev_round:
            starts.append((int(r["seq"]), int(r["round"])))
            prev_round = r["round"]
    # De deux marques trop proches on garde la SECONDE : le tour 0 n'a qu'une operation, et c'est
    # le debut du tour 1 qui interesse.
    return [m for i, m in enumerate(starts)
            if i + 1 == len(starts) or starts[i + 1][0] - m[0] >= gap]


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
    fig, ax = plt.subplots(figsize=(11, 5.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    # Séparateurs de tours, derrière les données et volontairement discrets : l'axe des x compte
    # des opérations, mais c'est par tour qu'on lit une dégradation.
    for seq, rnd in round_marks(df):
        ax.axvline(seq - 0.5, color=GRID, lw=1, zorder=0)
        ax.annotate(f"T{rnd}", xy=(seq - 0.5, 1.005), xycoords=("data", "axes fraction"),
                    color=INK_MUTED, fontsize=7.5, ha="left", va="bottom")

    ends = []
    for i, ts in enumerate(runs):
        color = SERIES[i]
        d = df[df.run_ts == ts].sort_values("seq")
        ax.plot(d.seq, d.prec_min, color=color, lw=2, zorder=3)
        if args.avg:
            ax.plot(d.seq, d.prec_avg, color=color, lw=2, ls=(0, (2, 2)), zorder=3)

        # Étiquette directe en bout de courbe : l'identité ne repose jamais sur la seule couleur,
        # et trois des huit teintes passent sous 3:1 sur fond clair.
        last = d.dropna(subset=["prec_min"]).iloc[-1]
        ends.append((float(last.prec_min), float(last.seq), color, short[ts]))

    ax.set_xlabel("n — n-ème opération du run", color=INK_MUTED, fontsize=9.5)
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

    # Place pour les etiquettes de bout de courbe, mesuree sur la plus longue. Les graduations,
    # elles, s'arretent aux donnees : une marge n'est pas une plage de valeurs.
    pad = 0.75 * max(len(s) for s in short.values()) + 2
    ax.set_xlim(df.seq.min() - 1, df.seq.max() + pad)
    ax.set_xticks([t for t in ax.get_xticks() if df.seq.min() <= t <= df.seq.max()])

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
        out = args.csv.rsplit(".", 1)[0] + suffix + ".png"
    fig.savefig(out, dpi=args.dpi, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"{len(runs)} run(s), {len(df)} lignes -> {out}")

    # La même chose en texte : une courbe sous 3:1 de contraste doit rester lisible autrement.
    piv = df.pivot_table(index="step", columns="run_ts", values=["prec_min", "prec_avg"], aggfunc="mean")
    print(piv.round(2).to_string())
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default="runs/aes.csv")
    ap.add_argument("--out", default=None,
                    help="PNG de sortie (defaut : le nom du CSV, plus les filtres). Incompatible avec --split")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--where", action="append", default=[], metavar="COL=VAL",
                    help="ne garder que les runs dont COL vaut VAL, repetable (ex: --where logn=11)")
    ap.add_argument("--split", default="logn", metavar="COL",
                    help="un PNG par valeur de COL trouvee dans le CSV. Defaut logn : deux degres "
                         "d'anneau ne se comparent pas sur un meme axe. 'none' pour tout reunir")
    ap.add_argument("--title", default=None, help="titre du graphique")
    ap.add_argument("--avg", action="store_true",
                    help="tracer aussi prec_avg, en pointilles (par defaut seul le pire slot)")
    ap.add_argument("--all", action="store_true",
                    help="tracer TOUS les runs, y compris plusieurs passes d'une meme configuration")
    args = ap.parse_args()

    if args.split not in ("", "none") and args.out:
        sys.exit("--split ecrit plusieurs PNG, il ne peut pas partager un seul --out")

    df = pd.read_csv(args.csv)
    if df.empty:
        sys.exit(f"{args.csv} : aucune ligne")

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
