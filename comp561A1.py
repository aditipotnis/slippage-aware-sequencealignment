def read_fasta(file_path):
    """Reads a FASTA file and returns two sequences."""
    sequences = []
    with open(file_path, 'r') as file:
        seq = ""
        for line in file:
            if line.startswith(">"):
                if seq:  # Save previous sequence
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        sequences.append(seq)  # Append the last sequence
    return sequences[0], sequences[1]

def initialize_substitution_matrix(match_score, mismatch_score):
    """Creates a 5x5 substitution matrix for nucleotides."""
    N = ['-','A', 'G', 'C', 'T']
    substitution_cost_matrix = {}

    for i in N:
        substitution_cost_matrix[i] = {}
        for j in N:
            if i == j:
                substitution_cost_matrix[i][j] = match_score
            else:
                substitution_cost_matrix[i][j] = mismatch_score 
    return substitution_cost_matrix

def slippage_aware_alignment(file_path, match_score, mismatch_score, cs, cn):
    S, T = read_fasta(file_path)
    m = len(S)
    n = len(T)

    M = initialize_substitution_matrix(match_score, mismatch_score)
    
    # Initialize the DP and direction matrices
    DP = [[0] * (n + 1) for _ in range(m + 1)]
    TB = [[None] * (n + 1) for _ in range(m + 1)]
    
    # Initialize the first column with proper gap penalties
    DP[0][0] = 0  # The starting point
    
    for i in range(1, m + 1):
        if i == 1:
            GAP_PENS = cn
        elif S[i - 1] == S[i - 2]:
            GAP_PENS = cs
        else:
            GAP_PENS = cn
        
        DP[i][0] = DP[i - 1][0] + GAP_PENS
        TB[i][0] = 'U'  # Up
    
    # Initialize the first row with proper gap penalties
    for j in range(1, n + 1):
        if j == 1:
            GAP_PENT = cn
        elif T[j - 1] == T[j - 2]:
            GAP_PENT = cs
        else:
            GAP_PENT = cn
        DP[0][j] = DP[0][j - 1] + GAP_PENT
        TB[0][j] = 'L'  # Left
    
    # Fill the DP and traceback matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate match/mismatch score
            D = DP[i - 1][j - 1] + M[S[i - 1]][T[j - 1]]

            # Calculate the gap penalty for moving up (insertion in T, deletion in S)
            if i > 1 and S[i - 1] == S[i - 2]:
                U = DP[i - 1][j] + cs  # Slippage condition
            else:
                U = DP[i - 1][j] + cn  # Non-slippage condition
            
            # Calculate the gap penalty for moving left (insertion in S, deletion in T)
            if j > 1 and T[j - 1] == T[j - 2]:
                L = DP[i][j - 1] + cs  # Slippage condition
            else:
                L = DP[i][j - 1] + cn  # Non-slippage condition
            
            # Choose the maximum score among match/mismatch, insertion, and deletion
            DP[i][j] = max(D, U, L)
            
            # Track the path
            if DP[i][j] == D:
                TB[i][j] = 'D'  # Diagonal
            elif DP[i][j] == U:
                TB[i][j] = 'U'  # Up
            else:
                TB[i][j] = 'L'  # Left
    
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    i, j = m, n
    while i > 0 or j > 0:
        if TB[i][j] == 'D':
            aligned_seq1 = S[i - 1] + aligned_seq1
            aligned_seq2 = T[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif TB[i][j] == 'U':
            aligned_seq1 = S[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        elif TB[i][j] == 'L':
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = T[j - 1] + aligned_seq2
            j -= 1
    
    return DP[m][n], aligned_seq1, aligned_seq2

# Example usage
match_score = 1   
mismatch_score = -1  
cs = -1  
cn = -2  

print("Alignment 1:/n")

# Run the alignment algorithm
alignment_score, aligned_seq1, aligned_seq2 = slippage_aware_alignment('/Users/aditipotnis/Downloads/sequence1.fa.txt', match_score, mismatch_score, cs, cn)

# Print the results
print("Optimal Alignment Score:", alignment_score)
print("Aligned Sequence 1:", aligned_seq1)
print("Aligned Sequence 2:", aligned_seq2)
print()


alignment_score2, aligned_seq12, aligned_seq22 = slippage_aware_alignment('/Users/aditipotnis/Downloads/sequence2.fa.txt', match_score, mismatch_score, cs, cn)

# Print the results
print()
print("Optimal Alignment Score:", alignment_score)
print("Aligned Sequence 1:", aligned_seq1)
print("Aligned Sequence 2:", aligned_seq2)
print()
