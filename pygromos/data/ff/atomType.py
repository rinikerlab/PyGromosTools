class AtomType(object):
    def __init__(self, atomNum, atomName, sqrtC06, sqrtC12_1, sqrtC12_2, sqrtC12_3, sqrtC06NB, sqrtC12NB, matrix):
        """

        Parameters
        ----------
        atomNum
        atomName
        sqrtC06
        sqrtC12_1
        sqrtC12_2
        sqrtC12_3
        sqrtC06NB
        sqrtC12NB
        matrix
        """
        self.atomNum = atomNum
        self.atomName = atomName
        self.sqrtC06 = sqrtC06
        self.sqrtC12_1 = sqrtC12_1
        self.sqrtC12_2 = sqrtC12_2
        self.sqrtC12_3 = sqrtC12_3
        self.sqrtC06NB = sqrtC06NB
        self.sqrtC12NB = sqrtC12NB
        self.matrix = matrix

    def __str__(self):
        return "[" + self.atomNum + "  " + self.atomName + "  " + self.sqrtC06 + "  " + self.sqrtC12_1 + "]"

    def addToMatrix(self, n):
        """

        Parameters
        ----------
        n

        Returns
        -------

        """
        self.matrix.append(n)
