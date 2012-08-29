import sys

# =============================================================================
def align(array1, array2, gap = -2, match = 1, mismatch = -1):
    """Performs Needleman-Wunsch alignment of string1 and string2.
    Prints out the alignment and returns the array of scores and pointers(arrows).

    Example usage from an interactive shell:
        from NeedlemanWunsch import NW
        Scores, Pointers = NW('PELICAN','COELACANTH')

    This is modified from a Perl implementation in the book BLAST by Korf, et al.
    """
    # initialize scoring and 'arrow' matrices to 0
    Scores = [[0 for x in range(len(array2)+1)] for y in range(len(array1)+1)]
    Pointers = [[0 for x in range(len(array2)+1)] for y in range(len(array1)+1)]

    # initialize borders
    # for pointers (arrows), use 2 for diagonal, -1 for horizontal, and 1 for vertical moves (an arbitrary system).
    # I have tried to consistently use i for rows (vertical positions) in the score and pointer tables, and j for columns (horizontal positions).
    for i in range(len(array1)+1):
        Scores[i][0] = gap*i
        Pointers[i][0] = 1 
    for j in range(len(array2)+1):
        Scores[0][j] = gap*j
        Pointers[0][j] = -1

    # fill with scores
    for i in range(1,len(array1)+1):
        for j in range(1,len(array2)+1):
            letter1 = array1[i-1]
            letter2 = array2[j-1]
            if letter1 == letter2: 
              DiagonalScore = Scores[i-1][j-1] + match
            else: DiagonalScore = Scores[i-1][j-1] + mismatch
            HorizontalScore = Scores[i][j-1] + gap 
            UpScore = Scores[i-1][j] + gap
            # TempScores is list of the three scores and their pointers
            TempScores = [[DiagonalScore,2],[HorizontalScore,-1],[UpScore,1]]
            # Now we keep the highest score, and the associated direction (pointer)
            Scores[i][j], Pointers[i][j] = max(TempScores)

    # backtrace from the last entry.  
    [i,j] = [len(array1),len(array2)]
    align1 = []
    align2 = []
    while [i,j] != [0,0]:
        if Pointers[i][j] == 2:
            align1.append(array1[i-1])
            align2.append(array2[j-1])
            i = i - 1
            j = j - 1
        elif Pointers[i][j] == -1:
            align1.append('-')
            align2.append(array2[j-1])
            j = j - 1
        else:
            align1.append(array1[i-1])
            align2.append('-')
            i = i - 1

    # the alignments have been created backwards, so we need to reverse them:
    align1 = align1[::-1]
    align2 = align2[::-1]

    # print out alignment
    #print align1
    #print align2

    # in case you want to look at the scores and pointers, the function returns them
    return [Scores,Pointers, [align1, align2]]

