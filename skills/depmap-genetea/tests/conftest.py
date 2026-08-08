"""Put the skill root on sys.path so `import kernel` works from tests/."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
