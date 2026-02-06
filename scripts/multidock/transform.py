""" Transform class """
from pathlib import Path
from typing import List


class Transform:
    """ Transform class, contains transformation files """
    def __init__(self) -> None:
        self.trans_file: List[str] = []
        self.trans_reverse: List[bool] = []
        self.size = 0

    def __len__(self) -> int:
        return self.size

    def __iter__(self):
        return iter(zip(self.trans_file, self.trans_reverse))

    def append(self, file: str, inverse: bool):
        """
        Append a new transformation
        Args:
            file (str): hdock output file
            inverse (bool): whether to fix the ligand

        Returns:
            Self: return self to chain call
        """
        self.trans_file.append(file)
        self.trans_reverse.append(inverse)
        self.size += 1
        return self

    def check(self):
        """ Check if all files are existed """
        for i in self.trans_file:
            if not Path(i).exists:
                raise FileNotFoundError
