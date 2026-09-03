import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# ============================================================
# Configuration
# ============================================================

# ------------------------------------------------------------
# Load PCM
# ------------------------------------------------------------

bg_vn = 12
n_simul = 78
bg_cn = 4
Zc = 6

code = np.load(
    "Codes/TCOM_aSCED/5G_zc=6/max_rank_5G_zc=6.npz"
)

H = code["h"].astype(int)


# ------------------------------------------------------------
# Determine transmitted PCM dimensions
# ------------------------------------------------------------

number_vn_simul = n_simul + 2 * Zc

number_vns_start = (
    bg_vn * Zc
    + 2 * Zc
)

harq_part = (
    number_vn_simul - number_vns_start
) // Zc

print(f"harq_part = {harq_part}")

m = Zc * bg_cn + harq_part * Zc +1

number_vns = (
    number_vns_start
    + harq_part * Zc
)

print(f"number_vns = {number_vns}")
print(f"H.shape   = {H.shape}")


# ============================================================
# Parameters
# ============================================================

# ------------------------------------------------------------
# Starting point
# ------------------------------------------------------------

# i is the LAST original / green row.
#
# The first blue row is therefore i + 1.

i = m - 2


# ------------------------------------------------------------
# Last non-grey column
# ------------------------------------------------------------

# j is the LAST non-grey column.

j = number_vns - 1


# ------------------------------------------------------------
# Delta / L sweep
# ------------------------------------------------------------

# Each tuple is:
#
#     (Delta, L)
#
# Delta = number of consecutive rows per batch
# L     = number of Delta-sized batches
#
# Example:
#
#     (2, 3)
#
# gives:
#
#     i+1, i+2       batch 1
#     i+3, i+4       batch 2
#     i+5, i+6       batch 3


delta_L_values = [
    (1, 1),
    (1, 2),
    (1, 6),
    (2, 1),
    (2, 3),
    (2, 6),
    (3,1),
    (3,2),
    (3,4),
    (6,1),
    (6,3),
]



# ------------------------------------------------------------
# Zoom window
# ------------------------------------------------------------

# Additional grey rows after the LARGEST selected region.

zoom_rows_grey = 10

# Additional grey columns after j.

zoom_cols_grey = 10


# ============================================================
# Determine maximum zoom size
# ============================================================

# Total number of blue rows:
#
#     Delta * L
#
# Find the largest Delta * L.

max_delta, max_L = max(
    delta_L_values,
    key=lambda x: x[0] * x[1],
)

max_selected_rows = (
    max_delta * max_L
)

print()
print("Delta / L sweep:")

for delta, L in delta_L_values:

    print(
        f"  Delta = {delta:2d}, "
        f"L = {L:2d}, "
        f"blue rows = {delta * L:3d}"
    )

print()
print(
    f"Maximum selection: "
    f"Delta = {max_delta}, "
    f"L = {max_L}, "
    f"{max_selected_rows} blue rows"
)


# ============================================================
# Output
# ============================================================

output_plot1 = "pcm_spy.pdf"

output_plot2 = (
    f"pcm_spy_i{i}_j{j}.pdf"
)


# ============================================================
# KIT colors
# ============================================================

KIT_GREEN = (
    0 / 255,
    150 / 255,
    130 / 255,
)

KIT_BLUE = (
    70 / 255,
    100 / 255,
    170 / 255,
)

GREY = "0.5"
GREY_ALPHA = 0.20


# ============================================================
# Appearance
# ============================================================

ENTRY_SIZE = 1.25

QC_LINE_WIDTH = 0.7
QC_LINE_ALPHA = 0.30

DELTA_LINE_WIDTH = 1.0
DELTA_LINE_ALPHA = 0.75


# ============================================================
# Axis setup
# ============================================================

def setup_axis(
    ax,
    H,
    row_min=0,
    row_max=None,
    col_min=0,
    col_max=None,
):

    m, n = H.shape

    if row_max is None:
        row_max = m - 1

    if col_max is None:
        col_max = n - 1

    ax.set_xlim(
        col_min - 0.5,
        col_max + 0.5,
    )

    ax.set_ylim(
        row_max + 0.5,
        row_min - 0.5,
    )

    ax.set_aspect("equal")

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)


# ============================================================
# Draw PCM entries
# ============================================================

def draw_entries(
    ax,
    H,
    row_min,
    row_max,
    col_min,
    col_max,
    i=None,
    j=None,
    delta=None,
    L=None,
):

    rows, cols = np.nonzero(H)

    # --------------------------------------------------------
    # Only draw entries inside displayed window
    # --------------------------------------------------------

    mask = (
        (rows >= row_min)
        & (rows <= row_max)
        & (cols >= col_min)
        & (cols <= col_max)
    )

    rows = rows[mask]
    cols = cols[mask]

    # --------------------------------------------------------
    # Determine selected region
    # --------------------------------------------------------

    if delta is not None and L is not None:

        # i = LAST green/original row
        selected_start = i + 1

        # selected_end is EXCLUSIVE
        #
        # Therefore:
        #
        # selected_start ... selected_end-1
        #
        # contains exactly Delta * L rows.

        selected_end = min(
            i + delta * L + 1,
            H.shape[0],
        )

    else:

        selected_start = None
        selected_end = None

    # --------------------------------------------------------
    # Draw entries
    # --------------------------------------------------------

    for row, col in zip(rows, cols):

        # ====================================================
        # Plot 1: complete PCM
        # ====================================================

        if i is None:

            facecolor = KIT_GREEN
            alpha = 1.0

        # ====================================================
        # Plot 2
        # ====================================================

        elif delta is None:

            # Everything below i OR right of j is grey.

            if row > i or col > j:

                facecolor = GREY
                alpha = GREY_ALPHA

            else:

                facecolor = KIT_GREEN
                alpha = 1.0

        # ====================================================
        # Plot 3
        # ====================================================

        else:

            # ------------------------------------------------
            # Column cutoff ALWAYS has priority
            # ------------------------------------------------

            if col > j:

                facecolor = GREY
                alpha = GREY_ALPHA

            # ------------------------------------------------
            # Rows after selected region
            # ------------------------------------------------

            elif row >= selected_end:

                facecolor = GREY
                alpha = GREY_ALPHA

            # ------------------------------------------------
            # Selected rows
            # ------------------------------------------------

            elif row >= selected_start:

                facecolor = KIT_BLUE
                alpha = 1.0

            # ------------------------------------------------
            # Original rows
            # ------------------------------------------------

            else:

                facecolor = KIT_GREEN
                alpha = 1.0

        # ----------------------------------------------------
        # Draw square
        # ----------------------------------------------------

        square = Rectangle(
            (
                col - ENTRY_SIZE / 2,
                row - ENTRY_SIZE / 2,
            ),
            ENTRY_SIZE,
            ENTRY_SIZE,
            facecolor=facecolor,
            edgecolor="none",
            alpha=alpha,
            zorder=3,
        )

        ax.add_patch(square)


# ============================================================
# QC block boundaries
# ============================================================

def add_qc_lines(
    ax,
    row_min,
    row_max,
    Zc,
):

    # Draw boundaries every Zc rows.

    first_boundary = (
        (row_min // Zc) + 1
    ) * Zc

    for row in range(
        first_boundary,
        row_max + 1,
        Zc,
    ):

        ax.axhline(
            row - 0.5,
            linestyle="--",
            linewidth=QC_LINE_WIDTH,
            alpha=QC_LINE_ALPHA,
            zorder=1,
        )


# ============================================================
# Delta batch boundaries
# ============================================================

def add_delta_lines(
    ax,
    i,
    delta,
    L,
):

    """
    i is the LAST original row.

    Blue rows are:

        i+1, ..., i+Delta
        i+Delta+1, ..., i+2*Delta
        ...

    Therefore the dotted lines are placed between batches.

    Example:

        i = 10
        Delta = 2
        L = 3

        row 10       green
        row 11       blue ┐
        row 12       blue ┘ batch 1
        ------------------- boundary
        row 13       blue ┐
        row 14       blue ┘ batch 2
        ------------------- boundary
        row 15       blue ┐
        row 16       blue ┘ batch 3
    """

    for batch in range(1, L):

        # Number of selected rows before this boundary:
        #
        #     batch * delta

        boundary_row = (
            i + 1 + batch * delta
        )

        ax.axhline(
            boundary_row - 0.5,
            linestyle=":",
            linewidth=DELTA_LINE_WIDTH,
            alpha=DELTA_LINE_ALPHA,
            zorder=4,
        )


# ============================================================
# Plot 1
# ============================================================

def create_plot1():

    m, n = H.shape

    fig, ax = plt.subplots(
        figsize=(8, 8),
    )

    setup_axis(
        ax,
        H,
        row_min=0,
        row_max=m - 1,
        col_min=0,
        col_max=n - 1,
    )

    draw_entries(
        ax,
        H,
        0,
        m - 1,
        0,
        n - 1,
    )

    plt.savefig(
        output_plot1,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.05,
    )

    plt.close(fig)

    print(f"Saved: {output_plot1}")


# ============================================================
# Plot 2
# ============================================================

def create_plot2():

    m, n = H.shape

    fig, ax = plt.subplots(
        figsize=(8, 8),
    )

    setup_axis(
        ax,
        H,
        row_min=0,
        row_max=m - 1,
        col_min=0,
        col_max=n - 1,
    )

    draw_entries(
        ax,
        H,
        0,
        m - 1,
        0,
        n - 1,
        i=i,
        j=j,
    )

    plt.savefig(
        output_plot2,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.05,
    )

    plt.close(fig)

    print(f"Saved: {output_plot2}")


# ============================================================
# Plot 3
# ============================================================

def create_plot3(
    delta,
    L,
):

    m, n = H.shape

    # --------------------------------------------------------
    # Selected region for THIS Delta / L
    # --------------------------------------------------------

    selected_start = i + 1

    # EXCLUSIVE upper bound

    selected_end = min(
        i + delta * L + 1,
        m,
    )

    # --------------------------------------------------------
    # Common zoom determined by the largest Delta * L
    # --------------------------------------------------------

    max_selected_end = min(
        i + max_selected_rows + 1,
        m,
    )

    # --------------------------------------------------------
    # Row window
    #
    # Keep complete upper part of PCM.
    #
    # Then show zoom_rows_grey additional rows after the
    # LARGEST selected region.
    # --------------------------------------------------------

    row_min = 0

    row_max = min(
        max_selected_end + zoom_rows_grey - 1,
        m - 1,
    )

    # --------------------------------------------------------
    # Column window
    # --------------------------------------------------------

    col_min = 0

    col_max = min(
        j + zoom_cols_grey,
        n - 1,
    )

    # --------------------------------------------------------
    # Output filename
    # --------------------------------------------------------

    output_plot3 = (
        f"pcm_spy_zoom_i{i}_j{j}"
        f"_delta{delta}_L{L}.pdf"
    )

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(10, 8),
    )

    setup_axis(
        ax,
        H,
        row_min=row_min,
        row_max=row_max,
        col_min=col_min,
        col_max=col_max,
    )

    # --------------------------------------------------------
    # QC block boundaries
    # --------------------------------------------------------

    add_qc_lines(
        ax,
        row_min,
        row_max,
        Zc,
    )

    # --------------------------------------------------------
    # PCM
    # --------------------------------------------------------

    draw_entries(
        ax,
        H,
        row_min,
        row_max,
        col_min,
        col_max,
        i=i,
        j=j,
        delta=delta,
        L=L,
    )

    # --------------------------------------------------------
    # Delta batch boundaries
    # --------------------------------------------------------

    add_delta_lines(
        ax,
        i,
        delta,
        L,
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    plt.savefig(
        output_plot3,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.08,
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------

    print(f"Saved: {output_plot3}")

    print(
        f"  Delta = {delta}, "
        f"L = {L}"
    )

    print(
        f"  selected rows: "
        f"{selected_start} ... {selected_end - 1}"
    )

    print(
        f"  number of blue rows: "
        f"{selected_end - selected_start}"
    )

    print(
        f"  common zoom rows: "
        f"{row_min} ... {row_max}"
    )

    print(
        f"  common zoom cols: "
        f"{col_min} ... {col_max}"
    )


# ============================================================
# Generate all plots
# ============================================================

create_plot1()

create_plot2()

for delta, L in delta_L_values:

    create_plot3(
        delta,
        L,
    )