"""
File: calculation of error estimates
Warnings: this class is WIP!

Description:
    Implementation of functions estimating the Error for calculated values using Gromos

Author: Paul Katzberger
"""

import numpy as np
from pygromos.utils.typing import List


class error_estimator:
    """
    Class to calculate the Error Estimate as implemented in ene_ana
    """

    def __init__(self, values: List[float]):
        """
        Initialize calculation
        Parameters
        ----------
        values : List[float]
            list of ordered values for which the error should be estimated
        """
        # Prepare blocks
        self.values = values
        self.blksz = 50
        self.old = 2
        self.d_counter = len(self.values)
        self.d_blocksize = []

    def calculate_rmsd(self):
        """
        Calculate rmsd
        Returns
        -------
        rmsd : float
            rmsd of values
        """
        sum = 0
        ssum = 0
        for v in self.values:
            sum += v
            ssum += v * v
        sum /= len(self.values)
        ssum /= len(self.values)
        msd = ssum - sum * sum

        return np.sqrt(msd)

    def calculate_error_estimate(self):
        """
        Calculation of the Error Estimate as in ene_ana
        Returns
        -------
        d_ee : float
            Error Estimate for provided list
        """
        # Setup Blocks
        while 4 * self.blksz < self.d_counter:
            self.d_blocksize.append(int(self.blksz))
            self.old = int(self.blksz)
            while self.old == int(self.blksz):
                self.blksz = self.blksz * 1.07177

        # Set start Variables
        Nblocks = len(self.d_blocksize)
        rmsd2 = 0
        ave = 0
        runave = np.mean(self.values)
        runrmsd = self.calculate_rmsd()

        fit = np.zeros(Nblocks)
        x = np.zeros(Nblocks)

        for j in range(Nblocks):
            Nblcki = self.d_counter // self.d_blocksize[j]

            rmsd2 = 0
            for i in range(Nblcki):
                start = i * self.d_blocksize[j]
                end = (i + 1) * self.d_blocksize[j]
                ave = np.mean(self.values[start:end])
                rmsd2 += (ave - runave) * (ave - runave)

            rmsd2 = rmsd2 / Nblcki

            fit[j] = (self.d_blocksize[j] * rmsd2) / runrmsd / runrmsd
            x[j] = 1 / self.d_blocksize[j]

        # Start addup
        sx, sf, sfx, sxx = 0, 0, 0, 0

        for i in range(Nblocks):
            sx += x[i]
            sf += fit[i]
            sfx += x[i] * fit[i]
            sxx += x[i] * x[i]

        a = (sf * sx / Nblocks - sfx) / (sx * sx / Nblocks - sxx)
        b = (sf - a * sx) / Nblocks

        d_ee = np.sqrt(b / self.d_counter) * runrmsd

        return d_ee
