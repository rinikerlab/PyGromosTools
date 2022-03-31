import sys
from pygromos.data.ff.atomType import AtomType


def main():
    fileName = "54a7.ifp"
    outName = "NEW_FILE.ifp"
    row1, row2, block = readFile(fileName)

    atomTypes = readBlock(block)
    addAtomType(atomTypes, 1)

    check(atomTypes)

    writeOut(row1, row2, atomTypes, fileName, outName)


def readFile(fileName):
    """

    Parameters
    ----------
    fileName

    Returns
    -------

    """
    lines = []
    row1 = 0
    row2 = 0

    with open(fileName, "r") as f:
        lines = f.readlines()

    lines = [row.strip().split() for row in lines]

    for i in range(len(lines)):
        if lines[i][0] == "SINGLEATOMLJPAIR":
            row1 = i
            row2, block = getBlock(lines, i + 1)
            return row1, row2, block


def getBlock(lines, row):
    """

    Parameters
    ----------
    lines
    row

    Returns
    -------

    """
    block = []

    for i in range(row, len(lines)):
        if lines[i][0] == "END":
            return i, block
        else:
            if lines[i][0] != "#":
                block.append(lines[i])


def readBlock(block):
    """

    Parameters
    ----------
    block

    Returns
    -------

    """
    nAtoms = 0
    atomTypes = []

    for i in range(len(block)):
        if block[i][0] == "#number":
            nAtoms = block[i + 1][0]  # noqa: F841

        elif block[i][0] == "#CS6":
            atm = getAtomType(block, i - 1)
            atomTypes.append(atm)

    return atomTypes


def getAtomType(block, row):
    """

    Parameters
    ----------
    block
    row

    Returns
    -------

    """
    atomBlock = []

    for i in range(row, len(block)):
        if block[i][0] == "#---":
            break
        else:
            if block[i][0] != "#CS6":
                atomBlock.append(block[i])

    atomNum = atomBlock[0][0]
    atomName = atomBlock[0][1]
    sqrtC06 = atomBlock[0][2]
    sqrtC12_1 = atomBlock[0][3]
    sqrtC12_2 = atomBlock[0][4]
    sqrtC12_3 = atomBlock[0][5]

    sqrtC06NB = atomBlock[1][0]
    sqrtC12NB = atomBlock[1][1]

    matrix = []
    for i in range(2, len(atomBlock)):
        for j in range(len(atomBlock[i])):
            matrix.append(atomBlock[i][j])

    atm = AtomType(atomNum, atomName, sqrtC06, sqrtC12_1, sqrtC12_2, sqrtC12_3, sqrtC06NB, sqrtC12NB, matrix)

    return atm


def addAtomType(atomTypes, N):
    """

    Parameters
    ----------
    atomTypes
    N

    Returns
    -------

    """
    for i in range(1, N + 1):
        # extend matrix of all existing atom type by 1
        for a in atomTypes:
            a.addToMatrix(1)

        # create dummy atom and add it to list of atomTypes
        numAtoms = len(atomTypes)
        num = str(numAtoms + 1)
        name = "MAP" + num
        val = "0.000000e-00"
        matrix = [1 for j in range(numAtoms)]
        matrix.append(1)
        atomTypes.append(AtomType(num, name, val, val, val, val, val, val, matrix))


def check(atomTypes):
    """

    Parameters
    ----------
    atomTypes

    Returns
    -------

    """
    N = len(atomTypes)

    for a in atomTypes:
        if len(a.matrix) != N:
            print(a)
            sys.exit("Matrix with wrong size\n")


def writeOut(row1, row2, atomTypes, fileName, outName):
    """

    Parameters
    ----------
    row1
    row2
    atomTypes
    fileName
    outName

    Returns
    -------

    """
    rawLines = []
    N = len(atomTypes)

    with open(fileName, "r") as f:
        rawLines = f.readlines()

    lines = [raw.strip().split() for raw in rawLines]

    out = open(outName, "w")
    for i in range(0, row2):
        if i <= row1:
            writeLine(rawLines[i], out)
        else:
            if lines[i][0] == "#" or lines[i][0] == "#number":
                if len(lines[i]) == 2:
                    if lines[i][1] == "NRATT":
                        writeLine(rawLines[i], out)
                        out.write("{0:>10}\n".format(N))
                else:
                    writeLine(rawLines[i], out)

    writeAtomTypes(atomTypes, out)

    for i in range(row2, len(lines)):
        writeLine(rawLines[i], out)

    out.close()


def writeLine(line, file):
    """

    Parameters
    ----------
    line
    file

    Returns
    -------

    """
    for item in line:
        file.write_out_energy_trajs("%s" % item)


def writeAtomTypes(atomTypes, out):
    """

    Parameters
    ----------
    atomTypes
    out

    Returns
    -------

    """
    for a in atomTypes:
        out.write_out_energy_trajs(
            "{0:>9} {1:>7} {2:>14} {3:>14} {4:>13} {5:>13}\n".format(
                a.atomNum, a.atomName, a.sqrtC06, a.sqrtC12_1, a.sqrtC12_2, a.sqrtC12_3
            )
        )

        out.write_out_energy_trajs("#CS6 CS12 parameters LJ14PAIR\n")
        out.write_out_energy_trajs("   {0:<12} {1:>14}\n".format((a.sqrtC06NB).strip(), a.sqrtC12NB))

        matrix = a.matrix

        for i in range(1, len(matrix) + 1):
            if i % 20 == 0 and i != 0:
                out.write_out_energy_trajs("   %s\n" % matrix[i - 1])
            else:
                out.write_out_energy_trajs("   %s" % matrix[i - 1])

        out.write_out_energy_trajs("\n#---\n")


if __name__ == "__main__":
    main()
