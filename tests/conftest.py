"""Put `notebooks/` on the import path so `import labcompass` works.

There is no package to install.  `labcompass.py` sits beside the loop folders,
where the notebooks import it from, and the tests reach it the same way.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "notebooks"))
