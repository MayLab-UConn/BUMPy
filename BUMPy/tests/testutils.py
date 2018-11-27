

class FileComp:

    @staticmethod
    def filesMatchExactly(file1, file2):
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            for line1, line2 in zip(f1, f2):
                if line1 != line2:
                    return False
        return True

    @staticmethod
    def filesMatchLinesStartingWith(file1, file2, lineStart):
        pass

    @staticmethod
    def filesMatchLinesContaining(file1, file2, lineStart):
        pass

    @staticmethod
    def compareFirstNChars(line1, line2, nLines):
        if len(line1) < nLines or len(line2) < nLines:
            return False
        return line1[:nLines] == line2[:nLines]


class PDBComp:

    @staticmethod
    def compareAtomFields(file1, file2):
        numCharsToCompare = 54  # pdb file, up to z coordinate
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            lines1 = [line[:numCharsToCompare] in line for line in f1 if line[:4] == "ATOM" ]
            lines2 = [line[:numCharsToCompare] in line for line in f2 if line[:4] == "ATOM" ]

            for line1 , line2 in zip(lines1, lines2):
                if line1 != line2:
                    return False
        return True


class GROComp:

    @staticmethod
    def compareAtomFields(file1, file2):
        numCharsToCompare = 44  # pdb file, up to z coordinate
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            # comment line
            next(f1)
            next(f2)
            nLinesF1 = int(next(f1))
            nLinesF2 = int(next(f2))
            if (nLinesF1 != nLinesF2):
                print("Line number mismatch, %d vs %d\n", nLinesF1, nLinesF2)
                return False

            lines1 = f1.readlines()
            lines2 = f2.readlines()

            for line1 , line2 in zip(lines1[:nLinesF1][:numCharsToCompare], lines2[:nLinesF1][:numCharsToCompare]):
                if line1 != line2:
                    return False
        return True


# class that stdout is redirected to. Modified from
# http://pragmaticpython.com/2017/03/23/unittesting-print-statements/
class stdout_checker:
    def __init__(self):
        self.data = []

    def write(self, string):
        self.data.append(string)

    def clear(self):
        self.data = []

    def __str__(self):
        return "".join(self.data)

    def flush(self):
        pass
