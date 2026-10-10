import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns
import os
import re
import glob
import networkx as nx
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc as sk_auc

# --- MASTER PUBLICATION STYLE GUIDE ---
# Every figure is drawn at its FINAL printed width (FULL_W below) rather than at some
# large working size that is later scaled down by the word processor. Scaling is what
# silently destroys figure legibility: an earlier version of this file drew figures
# 6-20 in wide and let the .docx shrink them to 5.85 in, which left body text at
# 3.7-4.4 pt on the page. Drawing at final size means the point sizes below ARE the
# point sizes a reader sees, and the smallest is 6.5 pt (journal floor is ~6 pt).
FULL_W = 5.83  # inches; measured from the .docx reference template's actual text
               # column (page width minus its margins), figures embedded at 100%

# The type scale below is set for legibility at final size, which means each figure
# must now budget its canvas for that type rather than the reverse. Two rules follow,
# and both were learned by measuring rendered output rather than by eye:
#   1. Vertical space is allocated from the text, not guessed. fig1_text_height()
#      computes how many points a stacked text block needs, and the figure is sized
#      from that, so labels cannot silently overlap when the canvas shrinks.
#   2. Adding type costs canvas. Ink coverage is the check: figures above ~12% non-
#      white pixels read as crowded once they are only 5.83 in wide. Where coverage
#      rose, content was cut or the canvas grown rather than the type shrunk.
FIG1_MIN_H = 2.9   # inches; clears the 5-box stack with padding
HEATMAP_MAX_ANNOT = 40   # above this, drop cell text and keep the colour scale
IMPORTANCE_TOP_N = 12    # bars shown; full ranking stays in the CSV

# Okabe-Ito: colorblind-safe by construction, so red/green confusion and
# deuteranopia cannot merge two series. Used everywhere instead of ad-hoc hues.
COLOR_CASE = '#D55E00'    # vermillion - IBD cases
COLOR_CTRL = '#0072B2'    # blue - controls
COLOR_ACCENT = '#009E73'  # bluish green - highlights
# Darkened from #999999: that grey only reached 2.85:1 on white, below the 4.5:1
# floor, which made the Figure 3 chance diagonal and gridlines hard to see.
COLOR_MUTED = '#5F6B76'   # 5.9:1 on white
COLOR_TEXT = '#1A1A1A'
PALETTE_CYCLE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9']

# Type scale, in true final-size points.
FS_PANEL = 8.5    # panel / subplot titles
FS_AXIS = 8.0     # axis labels
FS_TICK = 7.5     # tick labels
FS_ANNOT = 6.5    # in-plot annotations (lowest permitted)

plt.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'serif'],
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.titlesize': FS_PANEL,
    'axes.labelsize': FS_AXIS,
    'xtick.labelsize': FS_TICK,
    'ytick.labelsize': FS_TICK,
    'legend.fontsize': FS_TICK,
    'text.color': COLOR_TEXT,
    'axes.labelcolor': COLOR_TEXT,
    'axes.edgecolor': COLOR_TEXT,
    'xtick.color': COLOR_TEXT,
    'ytick.color': COLOR_TEXT,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'axes.linewidth': 0.7,
    'pdf.fonttype': 42,
})


def panel_label(ax, letter, dx=-0.14, dy=1.06):
    """Bold (a)/(b) marker so multi-panel figures can be cited unambiguously."""
    ax.text(dx, dy, f"({letter})", transform=ax.transAxes, fontsize=FS_PANEL + 1.5,
            fontweight='bold', va='bottom', ha='left')


def text_block_height(lines, fontsize, linespacing=1.2):
    """Height in points of a stacked n-line text block at a given size."""
    return lines * fontsize * linespacing


def text_width_pt(s, size, weight='normal', style='normal'):
    """Rendered width of a string in points, for real font metrics.

    Used instead of a character-count estimate because proportional fonts make
    'Harmonization' 25% wider than 'Validation' at the same size, and guessing is
    what let Figure 1's titles overrun their boxes.
    """
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    return TextPath((0, 0), s,
                    prop=FontProperties(size=size, weight=weight, style=style)
                    ).get_extents().width


def wrap_to_width(s, size, max_w, weight='normal', style='normal'):
    """Greedy word wrap so no rendered line exceeds max_w points.

    Words longer than max_w are returned on their own line rather than split, since
    a hyphenless mid-token break is less readable than a slightly wide line.
    """
    # Existing hard breaks are honoured: split() treats '\n' as whitespace, so a
    # flat split would silently reflow a caption the author had already broken and
    # then join the fragments with spaces, producing wider lines than intended.
    out = []
    for para in s.split('\n'):
        words, lines, cur = para.split(), [], ''
        for word in words:
            trial = f"{cur} {word}".strip()
            if cur and text_width_pt(trial, size, weight, style) > max_w:
                lines.append(cur)
                cur = word
            else:
                cur = trial
        if cur:
            lines.append(cur)
        out.append('\n'.join(lines))
    return '\n'.join(out)


def ink_coverage(path):
    """Split coverage into solid fills and glyph strokes.

    Total non-white fraction is a misleading crowding signal: a bar chart or a
    heatmap is mostly solid colour by construction (Figure 8 is 18.7% solid fill but
    only 4.1% text), so flagging on total ink condemns healthy figures. Glyph density
    is the part that competes for attention and predicts legibility, so that is what
    gets measured and thresholded.
    """
    from PIL import Image
    from scipy.ndimage import binary_erosion
    a = np.asarray(Image.open(path).convert('L'))
    ink = a < 240
    solid = binary_erosion(ink, np.ones((5, 5)))
    return {'solid': float(solid.mean()), 'glyph': float((ink & ~solid).mean())}


def report_coverage(path):
    """Print glyph density per figure; flag only genuinely text-heavy panels."""
    cov = ink_coverage(path)
    pct = cov['glyph'] * 100
    flag = "TEXT-CROWDED" if pct > 6.5 else "ok"
    print(f"    fill {cov['solid']*100:5.2f}%  glyph {pct:5.2f}%  {flag}")
    return pct

OUT_DIR = "/Users/edkwami/Desktop/new/figures_final"
os.makedirs(OUT_DIR, exist_ok=True)

SCORECARD = "out/replication_scorecard_species.csv"
COHORTS = ["Nielsen2014", "Ijaz2017", "HMP2"]

# All figures are derived from the analysis outputs. Nothing below is hard-coded:
# consensus taxa come from the meta-analysis scorecard, effect sizes from the same
# file, and ROC curves from the per-sample predictions written by ml_crosscohort.py.


def short_name(taxon):
    """Strip the Greengenes-style lineage down to the genus/species label."""
    if '.g__' in taxon:
        genus = taxon.split('.g__')[-1]
        if '.s__' in genus:
            genus = genus.split('.s__')[0]
        return genus
    return taxon.split('.s__')[-1].replace('_', ' ')


def full_label(taxon):
    """Human-readable 'Genus species' for a lineage."""
    genus = short_name(taxon)
    if '.s__' in taxon:
        sp = taxon.split('.s__')[-1].replace('_', ' ')
        return f"{genus} {sp}"
    return genus


def unique_labels(names):
    """Disambiguate a list of taxon labels so no two read identically.

    short_name() returns the genus alone, which is fine for prose but wrong on an
    axis: the 8 consensus taxa collapse to just five distinct genus strings
    (Roseburia x3, Lachnoclostridium x3), so a reader could not tell the bars, heatmap
    rows or network nodes apart. Where two labels collide after trimming, the full
    'Genus species' form is used for that label only, which keeps the short form
    wherever it is already unambiguous.
    """
    from collections import Counter
    counts = Counter(short_name(n) for n in names)
    return [full_label(n) if counts[short_name(n)] > 1 else short_name(n) for n in names]


def compact_label(taxon):
    """Abbreviated but unique label for tight spaces such as network nodes.

    Returns 'L. bolteae', not 'L. Clostridium bolteae'. Greengenes species strings
    already embed a genus name, so the lineage reads
    'g__Lachnoclostridium;s__Clostridium_bolteae': the genus is Lachnoclostridium but
    the species field is literally 'Clostridium bolteae'. Prefixing the genus initial
    onto that whole field produced 'L. Clostridium bolteae', which repeats a genus,
    runs 22-34 characters, and in a 2-inch network panel pushes labels into each
    other. The trailing epithet is the informative part, so it is used alone; any
    trailing clade/strain qualifier ('CAG 59', 'sp CAG 182') is kept because it is
    what distinguishes two otherwise identical binomials.
    """
    genus = short_name(taxon)
    if '.s__' not in taxon:
        return genus
    sp = taxon.split('.s__')[-1].replace('_', ' ').strip()
    words = sp.split()
    # Drop leading words that merely restate the genus. Suffix matching rather than
    # equality, because the species field usually carries a *coarser* genus than the
    # lineage does: g__Lachnoclostridium pairs with s__Clostridium_bolteae, so
    # 'Lachnoclostridium' != 'Clostridium' even though the latter is contained in it.
    while len(words) > 1 and genus.lower().endswith(words[0].lower()):
        words.pop(0)
    epithet = ' '.join(words)
    return f"{genus[0]}. {epithet}"


def load_scorecard():
    if not os.path.exists(SCORECARD):
        raise SystemExit(f"Missing {SCORECARD}; run ./run_all.sh analyze first.")
    s = pd.read_csv(SCORECARD)
    return s.rename(columns={'species': 'taxon'})


def _norm_key(name):
    """Alphanumeric-only key for matching taxon labels across sanitization layers.

    Maaslin2 rewrites feature names (spaces and dots become underscores), so the
    scorecard labels are not byte-identical to the row labels in the exported
    tables. Both sides are reduced to alphanumerics to join them reliably. All
    cohorts pass through the same sanitizer inside the analysis, so the
    meta-analysis itself is unaffected; only this figure-level join needs this.
    """
    return re.sub(r'[^A-Za-z0-9]', '', str(name))


def _resolve_taxa(taxa, cohort):
    """Map scorecard taxon labels onto the row labels actually present in a table.

    taxa=None means "keep every row", which is how Figure 8 obtains the full
    feature space rather than just the consensus taxa.
    """
    tax = pd.read_csv(f"data/taxon/{cohort}/combined_counts.csv", index_col=0)
    if taxa is None:
        return tax, None, None
    lookup = {}
    for raw in tax.index:
        lookup.setdefault(_norm_key(raw), raw)
    return tax, [_norm_key(t) for t in taxa], lookup


def load_consolidated_data(taxa=None):
    """Sample x taxon abundance matrix for the real cohorts, with status labels.

    taxa=None loads every taxon present in the cohort tables. That wider frame is
    what Figure 8 needs: the classifier's own feature set is defined by a >=20%
    prevalence rule over all taxa, so restricting the frame to the 8 consensus taxa
    would silently reduce the feature space the reported model was fitted on.
    """
    frames = []
    for cohort in COHORTS:
        meta_path = f"data/metadata/{cohort}.csv"
        if not os.path.exists(meta_path):
            continue
        tax, keys, lookup = _resolve_taxa(taxa, cohort)
        if keys is None:
            cols = [c for c in tax.index
                    if c not in ('unclassified', 'UNFRAC')]
        else:
            cols = [lookup[k] for k in keys if k in lookup]
        if not cols:
            continue
        meta = pd.read_csv(meta_path, index_col='sample_name')
        # Taxa are the INDEX of this table and samples are its columns, so select
        # rows with .loc[row_labels, :] (not .loc[:, cols], which is column-axis).
        df = tax.loc[cols, :].T.join(meta[['status']])
        # A taxon absent from a cohort table is a structural gap, not a zero count;
        # leaving it NaN lets each cohort contribute the taxa it actually measured
        # instead of imputing an absent taxon as undetected.
        df = df.loc[:, df.notna().any(axis=0)]
        df['cohort'] = cohort
        frames.append(df)
    if not frames:
        return None
    merged = pd.concat(frames)
    # Keep only taxa measured in every cohort, so a taxon is compared on the same
    # samples everywhere rather than being silently dropped from the cohorts that
    # lack it (which would bias cohort-wise means and correlations).
    common = [c for c in merged.columns
              if c not in ('status', 'cohort')
              and merged.groupby('cohort')[c].apply(lambda s: s.notna().any()).all()]
    return merged[common + ['status', 'cohort']]


def library_sizes():
    """Per-sample total counts across every taxon measured in that cohort.

    Computed from the cohort tables directly so it is independent of which subset of
    taxa a figure happens to display, and keyed by the same sample index the merged
    frame uses so the two align exactly.
    """
    sizes = {}
    for cohort in COHORTS:
        path = f"data/taxon/{cohort}/combined_counts.csv"
        if not os.path.exists(path):
            continue
        tax = pd.read_csv(path, index_col=0)
        tax = tax.loc[[c for c in tax.index if c not in ('unclassified', 'UNFRAC')]]
        sizes.update(tax.sum(axis=0).to_dict())
    return pd.Series(sizes)


def to_relative_abundance(df, taxa, library=None):
    """Convert counts to percent relative abundance, scaled by total library size.

    combined_counts.csv holds RAW COUNTS, not proportions (row sums differ by an
    order of magnitude between samples). The differential-abundance stage handles
    this correctly with Maaslin2 TSS+LOG, but this figure module used to read the
    counts directly, so any panel that compared samples or cohorts was really
    plotting sequencing depth, and the axis label claiming "relative abundance (%)"
    did not describe the data.

    The divisor must be the sample's TOTAL library size across every taxon measured
    in that cohort, not the sum over the handful of taxa being displayed. Dividing by
    a subset sum is wrong twice over: a taxon shown is then rescaled to look dominant
    if the other taxa are omitted, and samples whose displayed taxa are all absent
    produce 0/0 -> NaN. With only the 8 consensus taxa, 14 of 498 samples summed to
    zero and were silently dropped from every panel. So `library` is the per-sample
    total over all taxa, and `taxa` only selects which columns to return.
    """
    if library is None:
        library = df[[c for c in df.columns
                      if c not in ('status', 'cohort')]].sum(axis=1)
    # copy=True: asarray on a Series can hand back a read-only view, and the
    # sanitising assignment below would then raise.
    library = np.array(library, dtype=float, copy=True)
    library[~np.isfinite(library) | (library <= 0)] = np.nan
    counts = df[taxa].to_numpy(dtype=float)
    rel = pd.DataFrame(counts / library[:, None] * 100.0, index=df.index, columns=taxa)
    return rel


def consensus_taxa(scorecard, n=8):
    """Top meta-analytic features that replicated in >=2 cohorts.

    Replication is taken from the scorecard's k column (set by the REML stage), not
    from presence in every cohort table: with real data almost no feature is
    estimated in all three cohorts, so demanding that would leave a single hit and
    misrepresent how much of the signature actually generalizes.
    """
    rep = scorecard[scorecard['replicates'] == True].copy()
    rep['neglogp'] = -np.log10(rep['pval'].clip(lower=1e-300))
    top = rep.sort_values('neglogp', ascending=False)['taxon'].head(n)
    print(f"Consensus taxa (replicated in >=2 cohorts): {len(top)}")
    return list(top)


# -----------------------------------------------------------------------------
# FIGURE 1: THE PIPELINE (Process Flow)
# -----------------------------------------------------------------------------
def generate_figure_1():
    # Every box dimension and every text position below is derived from measured
    # rendered text width. Previously the boxes were a fixed 3.0 x 1.9 data units
    # and the type was laid on top, so 'Data Acquisition' (77.8 pt) rendered wider
    # than the 54.8 pt box it sat in, and the in/out captions (up to 91.6 pt for
    # 'out: Generalizable signature') spilled into the neighbouring stages. The
    # layout now sizes itself from the text and wraps anything that will not fit.
    # Full names are kept; the flow is set on two rows so each box gets ~2/3 of the
    # column instead of 1/5. A single row gave each stage ~68 pt, and real labels
    # like 'curatedMetagenomicData' (85.5 pt) are single words that cannot be
    # wrapped narrower, so a one-row layout can only fit them by abbreviating the
    # method names or the type, neither of which is worth the loss.
    stages = [
        ("Data\nAcquisition", "curatedMetagenomicData\nMetaPhlAn3 profiles",
         "3 IBD cohorts", "Harmonized CSVs"),
        ("Harmonization", "Lineage normalization\n& schema alignment",
         "Harmonized CSVs", "Standardized tables"),
        ("Differential\nAbundance", "Maaslin2 & ANCOM-BC2\n(dual-algorithm, TSS)",
         "Standardized tables", "Coef tables"),
        ("Meta-Analysis", "REML aggregation\n(random effects)",
         "Coef tables", "Replication scorecard"),
        ("Validation", "Random forest\nleave-one-cohort-out",
         "Replication scorecard", "Generalizable signature"),
    ]

    # Three boxes on the top row, two on the bottom, so each gets ~130 pt of the
    # 419.8 pt column. Wrapping still runs afterwards, but it now has room to work
    # and does not have to mangle method names to fit.
    per_row = [3, 2]
    arrow_pt, row_gap_pt = 12.0, 16.0
    col_pt = FULL_W * 72
    n = len(stages)

    def row_layout(count):
        return (col_pt - (count - 1) * arrow_pt) / count

    widths = [row_layout(c) for c in per_row]
    inner_pts = [wd - 8.0 for wd in widths]      # 4 pt padding each side

    row_of = [0, 0, 0, 1, 1]
    wrapped = []
    for (title, detail, inp, outp), ri in zip(stages, row_of):
        inner = inner_pts[ri]
        wrapped.append((
            wrap_to_width(title, FS_PANEL, inner, weight='bold'),
            wrap_to_width(detail, FS_ANNOT, inner, style='italic'),
            wrap_to_width(f"in: {inp}", FS_ANNOT, inner),
            wrap_to_width(f"out: {outp}", FS_ANNOT, inner),
        ))

    # Verify against real font metrics instead of trusting the wrapper. The row
    # index must come from row_of, not from `i < 3`, or the two rows get each
    # other's widths and the check is meaningless.
    over = []
    for i, grp in enumerate(wrapped):
        inner = inner_pts[row_of[i]]
        for tag, seg, size, wt, st in (
                ('title', grp[0], FS_PANEL, 'bold', 'normal'),
                ('detail', grp[1], FS_ANNOT, 'normal', 'italic'),
                ('in', grp[2], FS_ANNOT, 'normal', 'normal'),
                ('out', grp[3], FS_ANNOT, 'bold', 'normal')):
            for line in seg.split('\n'):
                wpt = text_width_pt(line, size, wt, st)
                if wpt > inner:
                    over.append((i, tag, line, wpt, inner))
    if over:
        for i, tag, line, wpt, inner in over:
            print(f"  OVERWIDE fig1 stage{i} {tag}: {wpt:.1f}pt > {inner:.1f}pt  {line!r}")
        raise SystemExit("Figure 1: text exceeds its box width after wrapping.")

    box_pad_pt, cap_gap_pt, title_gap_pt = 9.0, 5.0, 3.0
    box_h = (box_pad_pt
             + text_block_height(max(len(g[0].split('\n')) for g in wrapped),
                                 FS_PANEL, 1.15)
             + title_gap_pt
             + text_block_height(max(len(g[1].split('\n')) for g in wrapped),
                                 FS_ANNOT, 1.25)
             + box_pad_pt)
    in_h = text_block_height(max(len(g[2].split('\n')) for g in wrapped),
                             FS_ANNOT, 1.25)
    out_h = text_block_height(max(len(g[3].split('\n')) for g in wrapped),
                              FS_ANNOT, 1.25)
    # Captions live under their own box, so each row is box + gap + in + gap + out.
    # Caption stacks below a box occupy [yb - cap_gap - in_h - 3 - out_h, yb - cap_gap];
    # the previous code put the lower row's box at y = 3, so its captions landed at
    # negative y and were clipped. Set the lower row's box high enough for its
    # captions to clear the bottom margin, and the upper row's box high enough for
    # its captions to clear the lower row's box top.
    margin = 4.0
    row1_y0 = margin + out_h + 3.0 + in_h + cap_gap_pt
    caption_band = cap_gap_pt + in_h + 3.0 + out_h
    row0_y0 = row1_y0 + box_h + caption_band + 6.0
    height_in = (row0_y0 + box_h + margin) / 72.0

    # The axes must fill the figure exactly. plt.subplots leaves default
    # margins, so the axes is smaller than the figure and one data unit is
    # worth fewer than one point. Text renders at an absolute point size,
    # but box dimensions are data units, so every box came out scaled down
    # while its text stayed full size -- the text then overflowed and
    # overlapped no matter how the widths were tuned. Filling the figure
    # makes 1 data unit == 1 point, so the point-based layout below holds.
    fig = plt.figure(figsize=(FULL_W, height_in))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, col_pt)
    ax.set_ylim(0, height_in * 72)
    ax.axis('off')

    colors = ['#F2F6FA', '#DCE9F5', '#BBD5EC', '#8FBADD', '#2E6DA4']
    text_colors = ['#1A1A1A', '#1A1A1A', '#141414', '#0D0D0D', '#FFFFFF']

    title_h = text_block_height(max(len(g[0].split('\n')) for g in wrapped),
                                FS_PANEL, 1.15)
    detail_h = text_block_height(max(len(g[1].split('\n')) for g in wrapped),
                                 FS_ANNOT, 1.25)

    def place(i, cx, box_y0):
        title, detail, inp, outp = wrapped[i]
        ax.add_patch(patches.FancyBboxPatch(
            (cx - widths[row_of[i]] / 2, box_y0),
            widths[row_of[i]], box_h, boxstyle="round,pad=0.1",
            facecolor=colors[i], edgecolor='#2A2A2A', linewidth=0.9, zorder=2))
        ax.text(cx, box_y0 + box_h - box_pad_pt - title_h / 2, title,
                ha='center', va='center', fontsize=FS_PANEL, fontweight='bold',
                color=text_colors[i], zorder=3, linespacing=1.15)
        ax.text(cx, box_y0 + box_pad_pt + detail_h / 2, detail,
                ha='center', va='center', fontsize=FS_ANNOT, style='italic',
                color=text_colors[i], zorder=3, linespacing=1.25)
        ax.text(cx, box_y0 - cap_gap_pt - in_h / 2, inp, ha='center', va='center',
                fontsize=FS_ANNOT, color='#3A3A3A', zorder=3, linespacing=1.25)
        ax.text(cx, box_y0 - cap_gap_pt - in_h - 3.0 - out_h / 2, outp,
                ha='center', va='center', fontsize=FS_ANNOT, color='#2A2A2A',
                fontweight='bold', zorder=3, linespacing=1.25)

    centres = []
    k_in_row = [0, 1, 2, 0, 1]
    for i in range(n):
        if row_of[i] == 0:
            wd, y0, k = widths[0], row0_y0, k_in_row[i]
        else:
            wd, y0, k = widths[1], row1_y0, k_in_row[i]
        cx = k * (wd + arrow_pt) + wd / 2
        centres.append((cx, y0 + box_h / 2))
        place(i, cx, y0)

    # Arrows: across each row, then down from the end of row 1 to the start of row 2.
    for i in range(2):
        ax.annotate('', xy=(centres[i + 1][0] - widths[row_of[i]] / 2 - 1.5, centres[i][1]),
                    xytext=(centres[i][0] + widths[row_of[i]] / 2 + 1.5, centres[i][1]),
                    arrowprops=dict(arrowstyle='-|>', color='#3A3A3A', lw=1.0),
                    zorder=1)
    ax.annotate('', xy=(centres[3][0], centres[3][1] + box_h / 2 + 1.5),
                xytext=(centres[2][0], centres[2][1] - box_h / 2 - 1.5),
                arrowprops=dict(arrowstyle='-|>', color='#3A3A3A', lw=1.0,
                                connectionstyle="arc3,rad=-0.25"),
                zorder=1)

    plt.savefig(f"{OUT_DIR}/figure_7_flowchart.png", bbox_inches='tight',
                facecolor='white')
    plt.close()
    print("Figure 1 saved."); report_coverage(f"{OUT_DIR}/figure_7_flowchart.png")


# -----------------------------------------------------------------------------
# FIGURE 2: FOREST PLOT (effect sizes read from the meta-analysis scorecard)
# -----------------------------------------------------------------------------
def generate_figure_2(scorecard, taxa):
    df = scorecard[scorecard['taxon'].isin(taxa)].copy()
    df['label'] = df['taxon'].apply(full_label)
    df = df.sort_values('est', ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(FULL_W, 0.34 * len(df) + 1.5))
    ax.axvspan(-0.2, 0.2, color='#F4F4F4', zorder=0, label='Null zone ($\\beta \\approx 0$)')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.9, zorder=1)

    for i in range(len(df)):
        color = COLOR_CASE if df['est'][i] > 0 else COLOR_CTRL
        lo = df['est'][i] - 1.96 * df['se'][i]
        hi = df['est'][i] + 1.96 * df['se'][i]
        ax.errorbar(df['est'][i], i, xerr=[[df['est'][i] - lo], [hi - df['est'][i]]],
                    fmt='o', color=color, capsize=2.5, markersize=5.5, linewidth=1.4, zorder=3)
        p = df['pval'][i]
        ptxt = f"{p:.1e}" if p >= 1e-4 else "<1e-4"
        # Reserve the right-hand annotation strip in axes coordinates instead of
        # letting the text run past the last datum. With bbox_inches='tight' the
        # saved width is driven by whatever extends furthest right, and these
        # strings pushed the figure to 7.8 in against a 5.83 in column, so Word
        # scaled it to 0.75x and the 6.5 pt type came out at 4.9 pt.
        ax.text(1.02, i, f"{ptxt}  $I^2$={df['I2'][i]:.0f}%  k={int(df['k'][i])}",
                    transform=ax.get_yaxis_transform(), va='center',
                    fontsize=FS_ANNOT, color='#3A3A3A', clip_on=False)

    ax.set_yticks(range(len(df)))
    # Abbreviated taxa labels; the binomial is given in full in Figure 4 and the
    # scorecard, and the forest plot's job is comparing effect sizes.
    ax.set_yticklabels([compact_label(t) for t in df['taxon']], fontsize=FS_TICK)
    ax.set_xlabel(r'Pooled effect size ($\beta$, log scale; >0 = enriched in IBD)',
                  fontsize=FS_AXIS)
    # Colour carries meaning, so name it explicitly instead of relying on hue alone.
    ax.plot([], [], 'o', color=COLOR_CASE, markersize=5.5, label='Enriched in IBD')
    ax.plot([], [], 'o', color=COLOR_CTRL, markersize=5.5, label='Depleted in IBD')
    # The legend sat inside the axes at 'lower right', covering the three most negative
    # rows and their right-hand annotations. It is now a figure-level legend on its
    # own reserved strip below the plot, so it cannot overlap data by construction.
    # tight_layout is deliberately not used: it fights figure legends over the same
    # margin and silently pulls the legend back up into the axes.
    leg = fig.legend(loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=3,
                     fontsize=FS_ANNOT, frameon=False,
                     handletextpad=0.4, columnspacing=1.8, borderpad=0.0)
    ax.grid(axis='x', linestyle='--', alpha=0.35, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    sns.despine()
    # Reserve a fixed strip for the legend instead of letting tight_layout guess.
    fig.subplots_adjust(left=0.235, right=0.755, top=0.97, bottom=0.13)
    plt.savefig(f"{OUT_DIR}/figure_8_forest_plot.png", bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 2 saved."); report_coverage(f"{OUT_DIR}/figure_8_forest_plot.png")


# -----------------------------------------------------------------------------
# FIGURE 3: ROC CURVES from observed out-of-cohort predictions
# -----------------------------------------------------------------------------
def generate_figure_3():
    files = sorted(glob.glob("out/ml_results/predictions_*__*__*.csv"))
    if not files:
        print("Figure 3 skipped: no prediction files found.")
        return

    fig, ax = plt.subplots(figsize=(FULL_W, FULL_W * 0.86))
    ax.plot([0, 1], [0, 1], color=COLOR_MUTED, linestyle='--', linewidth=1.0,
            label='Random chance (AUC = 0.50)', zorder=1)

    summary = []
    dashes = [(None, None), (5, 2), (1.5, 1.5), (6, 2, 1, 2), (3, 1.5), (7, 2)]
    for i, f in enumerate(files):
        p = pd.read_csv(f)
        if p['y_true'].nunique() < 2:
            continue
        fpr, tpr, _ = roc_curve(p['y_true'], p['prob_ibd'])
        a = sk_auc(fpr, tpr)
        label = f"{p['train_cohort'].iloc[0]} $\\rightarrow$ {p['test_cohort'].iloc[0]}"
        # Colourblind-safe hues plus a distinct dash per curve, so the series stay
        # separable in greyscale and for any colour vision deficiency.
        ax.plot(fpr, tpr, label=f"{label}  (AUC {a:.3f}, n={len(p)})",
                color=PALETTE_CYCLE[i % len(PALETTE_CYCLE)], linewidth=1.5,
                dashes=dashes[i % len(dashes)] or (None, None), zorder=2 + i)
        summary.append({'train': p['train_cohort'].iloc[0], 'test': p['test_cohort'].iloc[0],
                        'auc': a, 'n_test': len(p)})

    ax.set_xlim([0, 1]); ax.set_ylim([0, 1.02])
    ax.set_xlabel('False positive rate', fontsize=FS_AXIS)
    ax.set_ylabel('True positive rate', fontsize=FS_AXIS)
    ax.legend(loc='lower right', fontsize=FS_ANNOT, frameon=True, framealpha=0.95,
              edgecolor='#CCCCCC', borderpad=0.4, labelspacing=0.35, handlelength=2.2)
    ax.grid(alpha=0.25, linestyle=':', linewidth=0.5)
    ax.set_axisbelow(True)
    sns.despine()
    fig.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_9_roc_curves.png", bbox_inches='tight', facecolor='white')
    plt.close()

    if summary:
        sdf = pd.DataFrame(summary).sort_values('auc', ascending=False)
        sdf.to_csv("out/ml_results/roc_auc_summary.csv", index=False)
        print(f"Figure 3 saved. Mean AUC = {sdf['auc'].mean():.3f} "
              f"(range {sdf['auc'].min():.3f}-{sdf['auc'].max():.3f}).")


# -----------------------------------------------------------------------------
# FIGURES 4-9
# -----------------------------------------------------------------------------
def generate_complex_analysis(taxa):
    df = load_consolidated_data(taxa)
    if df is None:
        print("Figures 4-9 skipped: no cohort data matched the consensus taxa.")
        return
    df_all = df
    # Figure 8 needs the classifier's own feature space (all taxa passing the >=20%
    # prevalence rule), which is wider than the 8 consensus taxa the other panels
    # show. Load it separately so Figures 4-7 and 9 keep their focused taxon set.
    try:
        df_all = load_consolidated_data(None)
    except Exception as e:  # noqa: BLE001 - degrade to the narrow frame
        print(f"Figure 8: full feature set unavailable ({e}); using consensus taxa.")
    # Use the column labels actually present in the merged frame (the raw table
    # row names), not the scorecard labels they were resolved from. The merged
    # frame also carries the 'status'/'cohort' annotation columns, which must be
    # excluded so groupby().mean() never sees string data.
    taxa = [c for c in df.columns if c not in ('status', 'cohort')]

# Normalise once, here. Everything below that compares samples or cohorts now
    # works on percent relative abundance, which is what the axis labels claim and
    # what makes the panels comparable across cohorts with different library sizes.
    # Scale by each sample's full library (all taxa measured in its cohort), not by
    # the consensus-taxon subset, so the percentages are true and no sample is lost.
    lib = library_sizes().reindex(df.index)
    rel = to_relative_abundance(df, taxa, library=lib)
    rel['status'] = df['status']
    rel['cohort'] = df['cohort']
    print(f"Figures 4-7/9: {len(rel)} samples normalized "
          f"({int(rel[taxa].isna().all(axis=1).sum())} unusable).")

    # Figure 4: composition, split by cohort so the per-taxon pattern is visible
    # instead of being averaged into one bar per cohort/status pair.
    cohorts_present = [c for c in COHORTS if c in set(rel['cohort'])]
    # Only the first panel carries taxon names. `sharey=True` means the other panels
    # draw ticks at the same positions, and matplotlib still reserves the label
    # gutter in every panel, so three identical copies of 'Lachnoclostridium
    # Clostridium clostridioforme' were squeezed into one third of the column and
    # overlapped vertically. Hiding them on panels b/c returns that width to the
    # bars; the rows line up because the y positions are shared.
    #
    # Row pitch is then set from measured line height rather than guessed: 8 taxa
    # in a 2.76 in figure gave 0.345 in = 24.8 pt per row against 9.8 pt of type,
    # which happens to fit, but the pitch was implicit in the row count. Labels are
    # wrapped to a measured width so no name can run into its neighbour or clip.
    row_h = 0.30
    n_panel = len(cohorts_present)
    label_pt = max(text_width_pt(l, FS_TICK)
                   for t in taxa for l in [compact_label(t)])
    inner_pt = text_width_pt("0.0", FS_ANNOT) + 6.0     # widest x tick + gap
    avail_pt = FULL_W * 72 - inner_pt - (n_panel - 1) * 26.0
    per_panel_pt = avail_pt / n_panel
    if label_pt > per_panel_pt:
        raise SystemExit(
            f"Figure 4: labels need {label_pt:.1f}pt but each panel has "
            f"{per_panel_pt:.1f}pt; shorten taxa or widen the column.")
    fig, axes = plt.subplots(1, n_panel, figsize=(FULL_W, row_h * len(taxa) + 1.5),
                             sharey=True)
    axes = np.atleast_1d(axes)
    # Build one label per taxon up front so the same string is used by both the
    # bar order and the tick labels; deriving them separately let the two disagree.
    label_map = dict(zip(taxa, unique_labels([compact_label(t) for t in taxa])))
    order = [label_map[t] for t in taxa]
    melted = rel.melt(id_vars=['cohort', 'status'], var_name='Taxon', value_name='RelAbund')
    melted['Label'] = melted['Taxon'].map(label_map)
    # One shared x-scale across all three panels. The bars show
    # MEANS, so the axis must be sized to the means (plus their
    # standard error), not to the largest individual sample. Sizing
    # it to the raw samples let the axis run to ~33% while the
    # tallest bar reached only ~2.3%, leaving the bars about 2 pt
    # wide -- effectively invisible. Sizing to mean+SE fills the
    # panel and keeps the cohorts comparable on one scale.
    agg = melted.groupby(['cohort', 'status', 'Taxon'])['RelAbund'].agg(['mean', 'sem'])
    xmax = float((agg['mean'] + agg['sem']).max()) * 1.15
    for i, (ax, coh) in enumerate(zip(axes, cohorts_present)):
        sub = melted[melted['cohort'] == coh]
        sns.barplot(data=sub, y='Label', x='RelAbund', hue='status', order=order,
                    palette={'PD': COLOR_CASE, 'control': COLOR_CTRL},
                    errorbar='se', ax=ax)
        ax.set_title(f"{coh}  (n={int((rel['cohort'] == coh).sum())})",
                     fontsize=FS_PANEL, fontweight='bold')
        ax.set_ylabel('')
        # barplot(x='RelAbund') labels the axis with the raw column
        # name, so each panel otherwise reads "RelAbund". Clear it;
        # one shared label below the row carries the real description.
        ax.set_xlabel('')
        ax.set_xlim(0, xmax)
        # sharey=True means all three panels draw from ONE y-axis, so
        # set_yticklabels([]) on panels b/c clears the shared labels and
        # leaves panel a blank. tick_params(labelleft=False) hides the
        # labels without touching the shared tick positions.
        if i == 0:
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels(order, fontsize=FS_TICK)
            ax.tick_params(axis='y', length=3, labelleft=True)
        else:
            ax.set_yticks(range(len(order)))
            ax.tick_params(axis='y', length=0, labelleft=False)
        ax.tick_params(axis='x', labelsize=FS_TICK)
        # Row pitch must clear the tallest label; assert rather than assume.
        lo, hi = ax.get_ylim()
        pitch = abs(hi - lo) / max(len(order), 1) * (fig.get_size_inches()[1] * 72)
        if pitch < FS_TICK * 1.15:
            raise SystemExit(
                f"Figure 4: row pitch {pitch:.1f}pt < label height "
                f"{FS_TICK*1.15:.1f}pt; labels would overlap vertically.")
        # Panel labels sat at dy=1.03, the same height as the panel
        # title, so "(b)" and "(c)" overlapped "Ijaz2017" and "HMP2".
        # Lift them above the titles and widen the top margin to fit
        # both. Panel a's label stays in the y-label gutter (dx=-0.32);
        # b and c sit at their left edge (dx=0.0).
        panel_label(ax, chr(ord('a') + i),
                    dx=-0.32 if i == 0 else 0.0, dy=1.115)
        sns.despine()
    # One legend only. seaborn draws a legend on every panel, so panels
    # b and c carried a second and third copy at their top right, sitting
    # on the bars, and a figure-level legend at the bottom collided with
    # the shared x-axis title. Panel a now uses a common x-scale, so its
    # bars (which top out near 6.5 against a 37.9 axis) leave the right
    # ~80% of that panel empty -- the legend goes there, in clear space.
    for ax in axes[1:]:
        if ax.get_legend():
            ax.get_legend().remove()
    leg = axes[0].get_legend()
    if leg is not None:
        leg.set_loc('center right')
        leg.get_frame().set_alpha(0.95)
        leg.get_frame().set_edgecolor('#CCCCCC')
        leg.get_title().set_fontsize(FS_ANNOT)
        for txt in leg.get_texts():
            txt.set_fontsize(FS_ANNOT)
    # One x-axis title for the whole row, centred below the three
    # panels. supxlabel anchors to the bottom of the subplot grid,
    # so the bottom margin must clear the tick labels; 0.145 of
    # the figure height leaves ~41 pt, enough for the ticks plus
    # this label without touching them.
    fig.supxlabel('Mean relative abundance (%)', fontsize=FS_AXIS,
                  y=0.02)
    fig.subplots_adjust(left=0.215, right=0.995, top=0.875, bottom=0.155, wspace=0.08)
    plt.savefig(f"{OUT_DIR}/figure_1_composition.png", bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 4 saved."); report_coverage(f"{OUT_DIR}/figure_1_composition.png")

    # Figure 5: PCA of relative-abundance profiles. Previously this ran on raw
    # counts, so PC1 largely encoded library size, and the in-figure title said
    # "PCA" while the manuscript caption said "PCoA". It is a linear projection of
    # the relative-abundance matrix, so it is PCA; the caption is aligned to match.
    X = rel[taxa]
    pca = PCA(n_components=2)
    comp = pca.fit_transform(X.fillna(0))
    fig, ax = plt.subplots(figsize=(FULL_W, 3.9))
    sns.scatterplot(x=comp[:, 0], y=comp[:, 1], hue=rel['status'], style=rel['cohort'],
                    palette={'PD': COLOR_CASE, 'control': COLOR_CTRL},
                    s=26, edgecolor='white', linewidth=0.4, alpha=0.85,
                    legend='full', ax=ax)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} of variance)",
                  fontsize=FS_AXIS)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} of variance)",
                  fontsize=FS_AXIS)
    ax.legend(fontsize=FS_ANNOT, title='Status / cohort', title_fontsize=FS_ANNOT,
              frameon=True, framealpha=0.95, edgecolor='#CCCCCC',
              loc='center left', bbox_to_anchor=(1.01, 0.5), borderpad=0.4,
              labelspacing=0.3, handletextpad=0.4)
    ax.axhline(0, color='#CCCCCC', linewidth=0.5)
    ax.axvline(0, color='#CCCCCC', linewidth=0.5)
    sns.despine()
    fig.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_2_pcoa.png", bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 5 saved."); report_coverage(f"{OUT_DIR}/figure_2_pcoa.png")

    # Figure 6: correlation heatmap. Spearman rather than Pearson, because
    # relative-abundance data are compositional and heavily zero-inflated, where
    # Pearson is dominated by a few large values.
    corr = rel[taxa].corr(method='spearman')
    lbl = unique_labels(list(corr.index))
    corr.index = lbl; corr.columns = lbl
    # Cell-by-cell numbers on an 8x8 grid drove coverage to ~16% and read as a wall of
    # digits. Above HEATMAP_MAX_ANNOT cells the numbers go and the colour scale
    # carries the magnitude, which is what the eye reads fastest anyway.
    show_annot = corr.size <= HEATMAP_MAX_ANNOT
    fig, ax = plt.subplots(figsize=(FULL_W, FULL_W * 0.76))
    sns.heatmap(corr, annot=show_annot, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                fmt=".2f", linewidths=0.4,
                annot_kws={"fontsize": FS_ANNOT}, ax=ax)
    # seaborn's heatmap returns the Axes, so reach the mappable via the collection
    # it drew in order to label the colourbar.
    cbar = ax.figure.colorbar(ax.collections[0], ax=ax, shrink=0.8)
    cbar.set_label('Spearman $\\rho$', fontsize=FS_ANNOT)
    cbar.ax.tick_params(labelsize=FS_ANNOT)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=FS_ANNOT)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=FS_ANNOT)
    fig.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_3_heatmap.png", bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 6 saved."); report_coverage(f"{OUT_DIR}/figure_3_heatmap.png")

    # Figure 7: alpha diversity per cohort (the caption says within-cohort, but the
    # old version pooled all three cohorts into one pair of boxes).
    taxa_cols = taxa
    with np.errstate(divide='ignore', invalid='ignore'):
        p = rel[taxa_cols].to_numpy(dtype=float) / 100.0
        shannon = -np.nansum(np.where(p > 0, p * np.log(p, where=p > 0), 0), axis=1)
    div = pd.DataFrame({'status': rel['status'].values, 'cohort': rel['cohort'].values,
                        'shannon': shannon})
    fig, axes = plt.subplots(1, len(cohorts_present), figsize=(FULL_W, 2.9), sharey=True)
    axes = np.atleast_1d(axes)
    for i, (ax, coh) in enumerate(zip(axes, cohorts_present)):
        sub = div[div['cohort'] == coh]
        sns.boxplot(x='status', y='shannon', data=sub, order=['control', 'PD'],
                    palette={'PD': COLOR_CASE, 'control': COLOR_CTRL},
                    width=0.55, fliersize=0, linewidth=0.8, ax=ax)
        sns.stripplot(x='status', y='shannon', data=sub, order=['control', 'PD'],
                      color='#333333', alpha=0.35, size=1.6, ax=ax)
        ax.set_title(f"{coh}\n(n={len(sub)})", fontsize=FS_PANEL, fontweight='bold')
        ax.set_xlabel('')
        ax.set_xticklabels(['Control', 'IBD'], fontsize=FS_TICK)
        panel_label(ax, chr(ord('a') + i), dx=-0.16, dy=1.14)
        if i:
            ax.set_ylabel('')
        else:
            ax.set_ylabel('Shannon diversity (natural log)', fontsize=FS_AXIS)
        sns.despine()
    fig.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_4_diversity.png", bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 7 saved."); report_coverage(f"{OUT_DIR}/figure_4_diversity.png")

    # Figure 8: feature importance. This deliberately mirrors the reported
    # classifier in 05_ml/ml_crosscohort.py (raw counts, >=20% prevalence filter in
    # the pooled frame) so the panel describes the model whose AUC is quoted in the
    # text, rather than a different model. sklearn's feature_importances_ is mean
    # decrease in impurity, not "Gini importance", so it is labelled accordingly.
    X_all = df_all[[c for c in df_all.columns if c not in ('status', 'cohort')]]
    prev = (X_all > 0).mean(axis=0)
    keep = prev[prev >= 0.20].index
    if len(keep) < 5:
        print("Figure 8 skipped: fewer than 5 taxa at >=20% prevalence.")
    else:
        # Fit on the same >=20%-prevalence feature set the reported classifier uses,
        # then display only the leading features. Restricting the fit to the 8
        # consensus taxa would not reproduce that model at all, and would leave too
        # few taxa surviving the prevalence filter to draw a ranking.
        y_ml = df_all['status'].map({'control': 0, 'PD': 1})
        rf = RandomForestClassifier(n_estimators=100, random_state=42,
                                    class_weight='balanced').fit(X_all[keep], y_ml)
        imp = pd.Series(rf.feature_importances_, index=keep).sort_values()
        top = imp.tail(min(IMPORTANCE_TOP_N, len(imp)))
        # Long binomial labels on every bar, plus 15 of them, pushed coverage to
        # ~22%. Abbreviating to 'F. saccharivorans' and trimming the tail cuts the
        # ink without losing the ranking the panel exists to show.
        bar_labels = unique_labels([compact_label(t) for t in top.index])
        fig, ax = plt.subplots(figsize=(FULL_W, 0.26 * len(top) + 1.1))
        ax.barh(range(len(top)), top.values, color=COLOR_ACCENT, edgecolor='none', height=0.68)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(bar_labels, fontsize=FS_TICK)
        ax.set_xlabel('Mean decrease in impurity (MDI)', fontsize=FS_AXIS)
        ax.grid(axis='x', alpha=0.25, linestyle=':', linewidth=0.5)
        ax.set_axisbelow(True)
        sns.despine()
        fig.tight_layout()
        plt.savefig(f"{OUT_DIR}/figure_5_importance.png", bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"Figure 8 saved (top {len(top)} of {len(keep)} taxa at >=20% prevalence).")
        report_coverage(f"{OUT_DIR}/figure_5_importance.png")

    # Figure 9: co-occurrence networks, Spearman on relative abundance.
    thr = 0.3
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, FULL_W * 0.68))
    for i, (status, title) in enumerate([('control', "Controls"), ('PD', "IBD cases")]):
        sub = rel[rel['status'] == status]
        adj = (sub[taxa].corr(method='spearman').abs() > thr).astype(int)
        mat = np.array(adj.values, dtype=int, copy=True)
        np.fill_diagonal(mat, 0)
        G = nx.from_pandas_adjacency(pd.DataFrame(mat, index=adj.index, columns=adj.columns))
        G.remove_nodes_from(list(nx.isolates(G)))
        # Keep every node on a circle. A spring layout on a 4-6 node hub-and-spoke
        # graph collapses into a line or a clump, and the label positions become
        # arbitrary; a fixed ring is reproducible and leaves an even gap between
        # every label, which is what a 2-inch panel needs.
        nodes = sorted(G.nodes(), key=lambda t: (-G.degree(t), compact_label(t)))
        pos = nx.circular_layout(G, scale=1.0)
        labels = {t: compact_label(t) for t in nodes}
        deg = dict(G.degree())
        # Node area was 280 + 150*degree pt^2, i.e. up to 730 pt^2 => 54 pt across,
        # more than a quarter of a 210 pt panel. Labels then sat on top of nodes.
        # Sizes are now scaled to the largest degree present so the marker spans a
        # fixed, modest fraction of the panel no matter how the degrees fall.
        dmax = max(deg.values()) or 1
        sizes = [180 + 260 * (deg[t] / dmax) for t in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=sizes,
                               node_color=[COLOR_CTRL if s == 'control' else COLOR_CASE
                                           for s in [status] * len(G)],
                               edgecolors='white', linewidths=0.9, ax=axes[i])
        nx.draw_networkx_edges(G, pos, width=0.9, edge_color='#9AA5B1', alpha=0.65, ax=axes[i])
        # Labels are pushed outside the ring so they never overlap a marker.
        for t in G.nodes():
            x, y = pos[t]
            axes[i].annotate(labels[t], xy=(x, y), xytext=(x * 1.32, y * 1.32),
                             ha='center', va='center', fontsize=FS_ANNOT,
                             fontweight='bold',
                             bbox=dict(boxstyle='round,pad=0.15', fc='white',
                                       ec='none', alpha=0.85))
        axes[i].set_title(f"{title}\nn={len(sub)}, {G.number_of_nodes()} taxa, "
                          f"{G.number_of_edges()} pairs", fontsize=FS_PANEL, fontweight='bold')
        axes[i].set_xlim(-1.75, 1.75)
        axes[i].set_ylim(-1.75, 1.75)
        axes[i].set_aspect('equal')
        axes[i].axis('off')
        panel_label(axes[i], chr(ord('a') + i), dx=-0.02, dy=0.99)
    fig.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_6_network.png", bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 9 saved."); report_coverage(f"{OUT_DIR}/figure_6_network.png")


if __name__ == "__main__":
    print("Generating publication figures from analysis outputs...")
    scorecard = load_scorecard()
    taxa = consensus_taxa(scorecard)
    print(f"Consensus taxa (replicated in >=2 cohorts): {len(taxa)}")
    for t in taxa:
        print(f"   - {full_label(t)}")
    generate_figure_1()
    generate_figure_2(scorecard, taxa)
    generate_figure_3()
    generate_complex_analysis(taxa)
    print(f"\nAll figures saved to: {OUT_DIR}")