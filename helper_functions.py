from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")
GAP_PENALTY = - 8

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    # Create a scoring matrix
    rows = len(seq1)
    cols = len(seq2)

    matrix = [[0] * (cols + 1) for _ in range(rows + 1)]

    # Initialise first column
    for row in range(1, rows + 1):
        matrix[row][0] = matrix[row - 1][0] + scoring_function(seq1[row - 1], "-")

    # Initialise first row
    for col in range(1, cols + 1):
        matrix[0][col] = matrix[0][col - 1] + scoring_function("-", seq2[col - 1])

    ## Fill in the matrix
    for row in range(1, rows + 1):
        for col in range(1, cols + 1):
            # Align current character between seq1 and seq2
            diagonal = matrix[row - 1][col - 1] + scoring_function(seq1[row - 1], seq2[col - 1])
            # Align current character from seq1 with a gap in seq2
            up = matrix[row - 1][col] + scoring_function(seq1[row - 1], "-")
            # Align current character from seq2 with a gap in seq1
            left = matrix[row][col - 1] + scoring_function("-", seq2[col - 1])
            # Choose the highest-scoring option
            matrix[row][col] = max(diagonal,up,left)

    score = float(matrix[rows][cols])

    # Traceback
    aligned1 = []
    aligned2 = []
    row = len(seq1)
    col = len(seq2)

    while row > 0 or col > 0:
        current = matrix[row][col]
        # Diagonal move
        if row > 0 and col > 0 and current == matrix[row - 1][col - 1] + scoring_function(seq1[row - 1], seq2[col - 1]):
            aligned1.append(seq1[row-1])
            aligned2.append(seq2[col-1])
            row -= 1
            col -= 1
        # Up move
        elif row > 0 and current == matrix[row - 1][col] + scoring_function(seq1[row - 1], "-"):
            aligned1.append(seq1[row - 1])
            aligned2.append("-")
            row -= 1
        # Left move
        elif col > 0 and current == matrix[row][col - 1] + scoring_function("-", seq2[col - 1]):
            aligned2.append(seq2[col-1])
            aligned1.append("-")
            col -= 1
        else:
            raise RuntimeError(f"Traceback failed)")

    # Traceback constructs the sequences backwards so reverse them
    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), score


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """

    # Create a scoring matrix
    rows = len(seq1) + 1 
    cols = len(seq2) + 1

    matrix = [[0] * cols for _ in range(rows)]

    # Keep track of the highest score and its location
    final_row = 0
    final_col = 0
    final_score = 0

    # Fill in the matrix
    for row in range(1, rows):
        for col in range(1, cols):
            # Align current character between seq1 and seq2
            diagonal = matrix[row - 1][col - 1] + scoring_function(seq1[row - 1], seq2[col - 1])
            # Align current character from seq1 with a gap in seq2
            up = matrix[row - 1][col] + scoring_function(seq1[row - 1], "-")
            # Align current character from seq2 with a gap in seq1
            left = matrix[row][col - 1] + scoring_function("-", seq2[col - 1])
            matrix[row][col] = max(diagonal, up, left)

            # If number is negative reset to 0
            if matrix[row][col] < 0:
                matrix[row][col] = 0

            # Remember the highest-scoring cell anywhere in the matrix
            if matrix[row][col] > final_score:
                final_score = matrix[row][col]
                final_row = row
                final_col = col

    # Traceback: start at the best cell, stop when hits 0
    aligned1 = []
    aligned2 = []
    row, col = final_row, final_col

    # Stop traceback when a cell with score 0 is reached or at the top left
    while row > 0 and col > 0 and matrix[row][col] > 0:
        current = matrix[row][col]
        # Diagonal move
        if current == matrix[row-1][col-1] + scoring_function(seq1[row-1], seq2[col-1]):
            aligned1.append(seq1[row-1])
            aligned2.append(seq2[col-1])
            row -= 1
            col -= 1
        # Up move
        elif current == matrix[row-1][col] + scoring_function(seq1[row-1], "-"):
            aligned1.append(seq1[row-1])
            aligned2.append("-")
            row -= 1
        # Left move
        elif current == matrix[row][col - 1] + scoring_function("-", seq2[col - 1]):
            aligned2.append(seq2[col-1])
            aligned1.append("-")
            col -= 1
        else:
            raise RuntimeError(f"Traceback failed)")

    # Build sequences forward
    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2), float(final_score)


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

def scoring_function_blosum(aa_i, aa_j):
    # Handle gaps separately
    if aa_i == "-" or aa_j == "-":
        return GAP_PENALTY
    # Return the BLOSUM62 substitution score
    return blosum62[aa_i, aa_j]

