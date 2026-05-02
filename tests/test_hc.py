import pathlib
import subprocess
import shutil

import pytest

BASE_DIR = pathlib.Path(__file__).parent.parent
INPUT_DIR = BASE_DIR / "sample_input"
SAMPLE_OUTPUT_DIR = BASE_DIR / "sample_output"
STUDENT_OUTPUT_DIR = BASE_DIR / "student_output"
BINARY_DIR = BASE_DIR / "build"
