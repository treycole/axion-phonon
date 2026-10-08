"""Put the repository root on sys.path for pytest.

Lets ``validation/unit`` import ``modules.*`` and ``validation.*`` without any
per-file path manipulation.
"""

import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parent.parent
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))
