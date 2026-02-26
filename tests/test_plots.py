import unittest

from pathlib import Path


class Test(unittest.TestCase):

    def setUp(self):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.testDir = Path(f"{os.getcwd()}/tests")
        self.cwd = Path(os.getcwd())

        self.baf = Baf()
    
    def cleanUp(self, in_dir: Path):
        """
            Remove directory and contents
        """
        for f in in_dir.glob("*"):
            f.unlink(missing_ok=True)
        in_dir.rmdir()