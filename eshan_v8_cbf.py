#!/usr/bin/env python3
import os
import sys
import time
import multiprocessing as mp
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt
from matplotlib import use as mpl_use
from matplotlib.patches import Rectangle
from matplotlib import font_manager
from Bio.PDB import PDBParser, DSSP, PPBuilder
try:
    # Biopython 1.80 or later
    from Bio.Data.PDBData import protein_letters_3to1
except ImportError:
    # Older Biopython versions (prior to 1.80, which still had three_to_one in Polypeptide, but for forward compatibility)
    from Bio.SeqUtils.IUPACData import protein_letters_3to1

# Define a function to use the dictionary for consistency with the original code's function call
def three_to_one(three_letter_code):
    """Converts a three-letter amino acid code to a one-letter code."""
    return protein_letters_3to1[three_letter_code.upper()]
from Bio.PDB import PDBParser
mpl_use("Agg")  # Safe for multiprocessing (no GUI)

# CONFIGURATION
models_dir = "./models"
x_offset = 80

# GFF PARSING HELPERS
def parse_attributes(attr_str):
    return dict(item.split("=", 1) for item in str(attr_str).split(";") if "=" in item)

def load_gff(gff_file):
    colnames = ["seqid","source","type","start","end","score","strand","phase","attributes"]
    gff = pd.read_csv(gff_file, sep="\t", comment="#", names=colnames, dtype=str)
    gff[["start","end"]] = gff[["start","end"]].astype(int)
    return gff

# PARSE GFF FOR A PARENT (full function from v5)
def parse_gff_for_parent(gff, parent_id):
    gene_rows = gff[gff["type"].isin(["gene","protein_coding_gene"])]
    gene_row = None
    gene_id = None
    for _, row in gene_rows.iterrows():
        attrs = parse_attributes(row["attributes"])
        if parent_id == attrs.get("ID") or parent_id == attrs.get("Name"):
            gene_row = row
            gene_id = attrs.get("ID")
            break
    if gene_row is None:
        for _, row in gene_rows.iterrows():
            attrs = parse_attributes(row["attributes"])
            name = attrs.get("Name","")
            if name.startswith(parent_id):
                gene_row = row
                gene_id = attrs.get("ID")
                break
    if gene_row is None:
        for _, row in gff[gff["type"].isin(["mRNA","transcript"])].iterrows():
            attrs = parse_attributes(row["attributes"])
            if parent_id == attrs.get("ID") or (attrs.get("Name","").startswith(parent_id)):
                gene_row = row
                gene_id = attrs.get("Parent")
                break
    if gene_row is None:
        raise ValueError(f"Parent ID {parent_id} not found in GFF.")
    strand = gene_row["strand"]
    gene_start, gene_end = int(gene_row["start"]), int(gene_row["end"])
    transcripts = []
    tx_rows = gff[gff["type"].isin(["mRNA","transcript"])]
    for _, row in tx_rows.iterrows():
        attrs = parse_attributes(row["attributes"])
        tx_id = attrs.get("ID")
        tx_parent = attrs.get("Parent")
        tx_name = attrs.get("Name", "")
        if tx_parent == gene_id or tx_id == parent_id or tx_name.startswith(parent_id):
            transcripts.append({"id": tx_id, "exons": [], "cds": []})
    if not transcripts:
        for _, row in tx_rows.iterrows():
            tx_start, tx_end = int(row["start"]), int(row["end"])
            if tx_start <= gene_end and tx_end >= gene_start:
                attrs = parse_attributes(row["attributes"])
                tx_id = attrs.get("ID")
                transcripts.append({"id": tx_id, "exons": [], "cds": []})
    tx_map = {tx["id"]: tx for tx in transcripts if tx["id"] is not None}
    exon_rows = gff[gff["type"] == "exon"]
    cds_rows = gff[gff["type"] == "CDS"]
    for _, row in exon_rows.iterrows():
        attrs = parse_attributes(row["attributes"])
        parent = attrs.get("Parent")
        if parent in tx_map:
            tx_map[parent]["exons"].append((int(row["start"]), int(row["end"])))
    for _, row in cds_rows.iterrows():
        attrs = parse_attributes(row["attributes"])
        parent = attrs.get("Parent")
        if parent in tx_map:
            source = "GFF"
            phase = row["phase"] if (not pd.isna(row["phase"]) and row["phase"] != "nan") else None
            tx_map[parent]["cds"].append((int(row["start"]), int(row["end"]), phase, source))
    for tx in transcripts:
        cds_list = tx["cds"]
        if not cds_list:
            continue
        cds_sorted = sorted(cds_list, key=lambda x: x[0], reverse=(strand == "-"))
        previous_frame = 0
        previous_length = 0
        enriched = []
        for start, end, ph, src in cds_sorted:
            length = end - start + 1
            if ph is None:
                ph_used = previous_frame
                src_used = "Inferred"
            else:
                ph_used = int(ph) if str(ph).isdigit() else previous_frame
                src_used = src
            frame = (previous_frame + previous_length) % 3
            enriched.append((start, end, ph_used, src_used, frame, previous_frame, previous_length))
            previous_frame = frame
            previous_length = length
        tx["cds"] = enriched
    ref_tx = None
    for tx in transcripts:
        if tx["cds"]:
            ref_tx = tx
            break
    if ref_tx is None:
        return gene_start, gene_end, strand, transcripts
    def get_first_cds_start(tx, strand):
        if not tx["cds"]:
            return None
        starts = []
        for s, e, *rest in tx["cds"]:
            starts.append(min(s, e))
        if strand == "+":
            return min(starts)
        else:
            return max(starts)
    ref_start_cds = get_first_cds_start(ref_tx, strand)
    if ref_start_cds is None:
        return gene_start, gene_end, strand, transcripts
    for tx in transcripts:
        tx_start = get_first_cds_start(tx, strand)
        if tx_start is None:
            continue
        offset = (ref_start_cds - tx_start) % 3
        if offset != 0:
            adjusted = []
            for (s, e, ph, src, frame, prev_frame, prev_len) in tx["cds"]:
                new_frame = (frame + offset) % 3
                adjusted.append((s, e, ph, src, new_frame, prev_frame, prev_len))
            tx["cds"] = adjusted
    return gene_start, gene_end, strand, transcripts

# GENE SCHEMATIC
# GENE SCHEMATIC PLOTTING
def plot_gene(ax, gene_start, gene_end, strand, transcripts, label=None):
    y_spacing = 0.5
    color_map = {0: "navy", 1: "orange", 2: "green"}

    for i, tx in enumerate(transcripts):
        y = - (i * y_spacing)
        cds_blocks = tx.get("cds", [])

        # Determine color: red if no protein model exists
        tx_id = tx.get("id", f"tx{i+1}")
        pdb_path = os.path.join(models_dir, f"{tx_id}.pdb")
        label_color = "red" if not os.path.exists(pdb_path) else "black"

        # Draw transcript label
        ax.text(gene_start - 300, y, tx_id, fontsize=10, ha="right", va="center", color=label_color)

        # draw exons if no CDS
        if not cds_blocks:
            for (s, e) in tx.get("exons", []):
                ax.add_patch(Rectangle((s, y - 0.08), e - s, 0.16,
                                       facecolor="skyblue", edgecolor="black", lw=0.4))
            continue

        # draw CDS blocks
        cds_blocks_sorted = sorted(cds_blocks, key=lambda x: x[0], reverse=(strand == "-"))
        for j in range(len(cds_blocks_sorted) - 1):
            end_j = cds_blocks_sorted[j][1]
            start_next = cds_blocks_sorted[j + 1][0]
            ax.plot([end_j, start_next], [y, y], color="black", lw=0.8)

        for (start, end, phase, source, frame, prev_frame, prev_len) in cds_blocks_sorted:
            ax.add_patch(Rectangle((start, y - 0.08), end - start, 0.16,
                                   facecolor=color_map.get(frame, "gray"), edgecolor="black", lw=0.4))

        # draw arrow for strand direction
        if strand == "+":
            arrow_start = cds_blocks_sorted[-1][1]
            ax.arrow(arrow_start, y, 40, 0, head_width=0.18, head_length=50,
                     fc="black", ec="black", lw=0, length_includes_head=True)
        else:
            arrow_start = cds_blocks_sorted[-1][0]
            ax.arrow(arrow_start, y, -40, 0, head_width=0.18, head_length=50,
                     fc="black", ec="black", lw=0, length_includes_head=True)

    ax.set_xlim(gene_start - 400, gene_end + 200)  # more space for labels
    ax.set_ylim(-len(transcripts) * y_spacing - 0.2, 0.5)
    ax.set_yticks([])

    # Keep only the x-axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_visible(False)

    if label:
        ax.set_title(label, fontsize=14)

# SIMPLE PDB CARTOON RENDERING WITH SECONDARY STRUCTURE + pLDDT
def render_pdb_cartoon(pdb_path, offset_x=0, width=3000, height=1500):
   # Render a PDB file as a cartoon (CA trace) with colouring by pLDDT.
   # Returns the temporary PNG path.

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("model", pdb_path)
    except Exception as e:
        print(f"Failed to parse PDB {pdb_path}: {e}")
        return None

    coords = []
    colors = []

    def plddt_color(atom):
        try:
            b = float(atom.get_bfactor())
        except Exception:
            b = 0.0
        if b >= 90:
            return (0/255, 76/255, 202/255)
        elif b >= 70:
            return (73/255, 196/255, 238/255)
        elif b >= 50:
            return (255/255, 213/255, 57/255)
        else:
            return (255/255, 113/255, 67/255)

    for model in structure:
        for chain in model:
            chain_coords = []
            chain_colors = []
            for res in chain:
                if "CA" in res:
                    atom = res["CA"]
                    chain_coords.append(atom.coord[:2])
                    chain_colors.append(plddt_color(atom))
            if chain_coords:
                coords.append(chain_coords)
                colors.append(chain_colors)

    if not coords:
        print(f"No CA atoms found in {pdb_path}")
        return None

    # compute bounds
    all_x = [x for chain in coords for x, y in chain]
    all_y = [y for chain in coords for x, y in chain]
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    margin_x = (x_max - x_min) * 0.05
    margin_y = (y_max - y_min) * 0.05

    fig_height = max(4, height/300 * ((y_max - y_min + 2*margin_y)/(x_max - x_min + 2*margin_x)))
    fig, ax = plt.subplots(figsize=(width/300, fig_height))
    ax.axis("off")
    ax.set_aspect("equal")
    ax.set_xlim(x_min - margin_x + offset_x, x_max + margin_x + offset_x)
    ax.set_ylim(y_min - margin_y, y_max + margin_y)

    for chain_coords, chain_colors in zip(coords, colors):
        xs, ys = zip(*chain_coords)
        for i in range(len(xs)-1):
            ax.plot([xs[i]+offset_x, xs[i+1]+offset_x],
                    [ys[i], ys[i+1]],
                    color=chain_colors[i], lw=2)
        ax.scatter([x+offset_x for x in xs], ys, color=chain_colors, s=12, zorder=3)

    tmp_png = pdb_path.replace(".pdb", "_tmp.png")
    plt.tight_layout()
    plt.savefig(tmp_png, dpi=300, bbox_inches="tight")
    plt.close()
    return tmp_png

# MAKE GENE PNG
def make_gene_png_for_parent(parent_id, gff, out_png):
    gene_start, gene_end, strand, transcripts = parse_gff_for_parent(gff, parent_id)
    fig, ax = plt.subplots(figsize=(16, max(2, 0.6 * max(1, len(transcripts)))))
    plot_gene(ax, gene_start, gene_end, strand, transcripts, label=parent_id)
    plt.tight_layout()
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close()
    return [tx["id"] for tx in transcripts if tx.get("id")]



# WORKER
def process_parent(args):
    parent_id, gff = args
    try:
        # Generate gene schematic PNG
        gene_png = f"{parent_id}_schematic.png"
        transcript_ids = make_gene_png_for_parent(parent_id, gff, gene_png)

        # Render structures (or placeholders)
        structure_pngs = []
        width_per_model = 800  # width of each model image
        height_per_model = 150  # default height of placeholder

        for tx in transcript_ids:
            pdb_path = os.path.join(models_dir, f"{tx}.pdb")
            if os.path.exists(pdb_path):
                tmp_png = render_pdb_cartoon(pdb_path)
            else:
                # create blank placeholder
                tmp_png = f"{tx}_placeholder.png"
                img = Image.new("RGB", (width_per_model, height_per_model), (255,255,255))
                img.save(tmp_png)
            structure_pngs.append(tmp_png)

        # Combine all structure images horizontally with aligned bottom labels
        imgs = [Image.open(p) for p in structure_pngs]
        if not imgs:
            imgs = [Image.new("RGB", (width_per_model, height_per_model), (255,255,255))]

        max_model_h = max(im.height for im in imgs)

        # Font settings
        label_font_size = 60  # moderate font size for clarity
        font_path = font_manager.findfont("DejaVu Sans:style=normal")
        try:
            font = ImageFont.truetype(font_path, label_font_size)
        except:
            font = ImageFont.load_default()

        # Add extra space for labels below the models
        total_w = sum(im.width for im in imgs)
        total_h = max_model_h + label_font_size + 40  # margin for labels
        combined = Image.new("RGB", (total_w, total_h), (255,255,255))
        draw = ImageDraw.Draw(combined)

        x_offset = 0
        y_bottom = max_model_h + 20  # vertical position for labels
        for i, im in enumerate(imgs):
            # Paste model / placeholder
            combined.paste(im, (x_offset, 0))

            # Transcript ID label
            tx_id = transcript_ids[i]
            color = "red" if not os.path.exists(os.path.join(models_dir, f"{tx_id}.pdb")) else "black"

            # Compute horizontal center
            x_pos = x_offset + im.width // 2
            draw.text((x_pos, y_bottom), tx_id, fill=color, font=font, anchor="ms")  # "ms" = middle, baseline

            x_offset += im.width

        structure_png = f"{parent_id}_structures.png"
        combined.save(structure_png)

        # Combine gene schematic + structures vertically
        img1 = Image.open(gene_png)
        img2 = Image.open(structure_png)
        w, h = max(img1.width, img2.width), img1.height + img2.height
        final = Image.new("RGB", (w, h), (255,255,255))
        final.paste(img1, (0,0))
        final.paste(img2, (0,img1.height))
        combined_png = f"{parent_id}_combined.png"
        final.save(combined_png)

        return f"{parent_id} done."

    except Exception as e:
        return f"{parent_id} error: {e}"

###########################################################################
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python eshan_v6.py <GFF_FILE> <PARENT_LIST_FILE>")
        sys.exit(1)
    gff_file, parent_list_file = sys.argv[1], sys.argv[2]
    gff = load_gff(gff_file)
    with open(parent_list_file) as f:
        parent_ids = [line.strip() for line in f if line.strip()]
    start_time = time.time()
    ncpu = min(mp.cpu_count(), len(parent_ids))
    print(f"Using {ncpu} CPU cores...\n")
    with mp.Pool(ncpu, maxtasksperchild=1) as pool:
        for result in pool.imap_unordered(process_parent, [(pid, gff) for pid in parent_ids]):
            print(result)
    print(f"\n⏱ Total runtime: {time.time()-start_time:.1f} s")

