# Motif-mark
- BGMP class project.
  
## Goal: Visualize motifs on sequences

## Requirements:
- Python 3 program that uses:
    - argparse
    - pycairo
- Program should handle
    - Ambiguous motifs
    - Multiple motifs
    - Multiple sequences
- Input
    - Fasta file
        - sequence lowercase alpha, except exons (uppercase alpha)
    - Motifs txt file (1 motif per line)
- Output
    - SVG created using pycairo
    - Graphic should include labeling
