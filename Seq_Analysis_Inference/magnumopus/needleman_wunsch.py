
                            ###############################
                            ##### alignment functions #####
                            ###############################


def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    """
    aligns two sequences and returns alignments, score, and pairwise distance
    """
    rows, cols = len(seq_a) + 1, len(seq_b) + 1
    empty_matrix = construct_matrix(rows, cols)
    initialized_matrix = initialize_matrix(empty_matrix, rows, cols, gap)
    filled_matrix, score = fill_matrix(initialized_matrix, seq_a, seq_b, match, mismatch, gap)
    alignments, score, pairwise_distance = backtrace_matrix(filled_matrix, score, seq_a, seq_b, match, mismatch, gap)
    return alignments, score, pairwise_distance


def construct_matrix(rows:int, cols:int) -> list[list[int]]:
    """
    create an empty matrix 
    """
    empty_matrix = []
    for r in range(rows):
        row = []
        for c in range(cols):
            row.append(0)
        empty_matrix.append(row)
    return empty_matrix


def initialize_matrix(matrix: list[list[int]], rows:int, cols:int, gap:int) -> list[list[int]]:
    """
    initialize matrix by filling first row and column using gap scores 
    """
    #initialize first column 
    for r in range(1, rows):
        matrix[r][0] = r * gap
    #initialize first row 
    for c in range(1, cols):
        matrix[0][c] = c * gap 
    return matrix 


def fill_matrix(matrix: list[list[int]], seq_a: str, seq_b: str, match: int, mismatch: int, gap:int)  -> tuple[list[list[int]], int]:
    """
    fill matrix based on match, mismatch, and gap scores
    """
    score = 0
    rows, cols = len(seq_a) + 1, len(seq_b) + 1
    for i in range(1, rows):
        for j in range(1, cols):

            if seq_a[i-1] == seq_b[j-1]:
                diagonal = matrix[i-1][j-1] + match
            else:
                diagonal = matrix[i-1][j-1] + mismatch
            
            delete = matrix[i-1][j] + gap
            insert = matrix[i][j-1] + gap
            matrix[i][j] = max(diagonal, delete, insert)

    score = matrix[-1][-1]
    return matrix, score


def backtrace_matrix(filled_matrix, score: tuple[list[list[int]], int], seq_a: str, seq_b: str, match: int, mismatch: int, gap:int):
    """
    find best alignment, score, and pairwise_distance
    """
    a_aligned = []
    b_aligned = []
    rows, cols = len(seq_a), len(seq_b)
    score = score
    pairwise_distance = 0 
    
    while rows > 0 or cols > 0:
        current = filled_matrix[rows][cols]
        if seq_a[rows-1] == seq_b[cols-1]:
            m_score = match 
        else:
            m_score = mismatch 

        #diagonal (match and mismatch)
        if rows > 0 and cols > 0 and current == filled_matrix[rows-1][cols-1] + m_score:
            a_aligned.append(seq_a[rows-1])
            b_aligned.append(seq_b[cols-1])
            if seq_a[rows-1] != seq_b[cols-1]:
                pairwise_distance += 1
            rows -= 1
            cols -= 1
        #up (gap in seq_b)
        elif rows > 0 and current == filled_matrix[rows-1][cols] + gap:
            a_aligned.append(seq_a[rows-1])
            b_aligned.append('-')
            pairwise_distance += 1
            rows -= 1
        #left (gap in seq_a)
        else:
            a_aligned.append('-')
            b_aligned.append(seq_b[cols-1])
            pairwise_distance += 1
            cols -= 1
            
    aligned_1 = ''.join(reversed(a_aligned))
    aligned_2 = ''.join(reversed(b_aligned))
    return (aligned_1, aligned_2), score, pairwise_distance

