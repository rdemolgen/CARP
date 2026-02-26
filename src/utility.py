import logging, os, re, time
from pathlib import Path

class Utility:

    def __init__(self, logPath, prefix, proband_id, now):
        # self.catch_the_carp(Path("resources") / "carp")
        self.logger = self.setup_logging(logPath, prefix, proband_id, now)

    def catch_the_carp(self, file_path: Path, delay: float = 0.001):
        """
            Show the operator the carp!
        """
        with file_path.open("r") as f:
            while True:
                char = f.read(1)
                if not char:
                    break
                print(char, end="", flush=True)
                time.sleep(delay)

    def setup_logging(self, logPath: Path, prefix: str, proband_id: str, now: str) -> logging.Logger:
        """Create log file in specificied folder"""
        if not logPath.is_dir():
            logPath.mkdir(parents=True, exist_ok=True)
            logPath.chmod(0o775) 
        
        logName = f"{prefix}{proband_id}_{now}"
        logger = logging.getLogger(logName)
        logger.setLevel(logging.INFO)
        logger.propagate = False

        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

        file_handler = logging.FileHandler(logPath / f"{logName}.log")
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(name)s: %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        logger.info(f"✅ Log file created: '{logName}.log'")

        return logger   

    def set_capture(self, proband_id: str) -> str:
        """
            Set capture based on sample id prefix
        """
        captures = {
            "TwEx": "exome",
            "WGS": "genome",
        }

        for prefix, capture in captures.items():
            if proband_id.startswith(prefix):
                self.logger.info(f"Capture: {capture}")
                return capture

        msg = f"Unknown capture type for sample: {proband_id}"
        self.logger.error(msg)
        raise ValueError(msg)
    
    def verify_dir(self, dirStr: str) -> Path:
        """Set input/output direct"""
        try:
            dir = Path(dirStr)
        except TypeError:
            self.logger.error(f"Directory: '{dirStr}' not provided or invalid")
            dir = Path(os.getcwd())

        if not dir.is_dir():
            self.logger.info(f"{str(dir)} folder created")
            dir.mkdir(mode=0o777, parents=True, exist_ok=False)

        return dir

    def find_file(self, in_dir: Path, regex_pattern: str) -> Path:
        """
            Find a file based on the provided regex_pattern.
            Error if no file is found or if duplicates exist.
        """
        self.logger.info(f"Searching for '{regex_pattern}' in '{str(in_dir)}'")
        pattern = re.compile(regex_pattern)

        reg_files = [p for p in in_dir.iterdir() if p.is_file() and pattern.match(p.name)]

        if not reg_files:
            msg = f"No file found in {in_dir} matching {regex_pattern!r}"
            self.logger.error(msg)
            raise FileNotFoundError(msg)

        if len(reg_files) > 1:
            msg = f"Multiple files found in {in_dir} matching {regex_pattern!r}: " + ", ".join(p.name for p in reg_files)
            self.logger.error(msg)
            raise RuntimeError(msg)

        self.logger.info(f"File found:  '{reg_files[0]}'")
        return reg_files[0]