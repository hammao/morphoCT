import numpy as np
import pytest


@pytest.fixture(scope="session")
def blobs():
    import porespy as ps

    return ps.generators.blobs(shape=[48, 48, 48], porosity=0.3, blobiness=0.5, seed=3).astype(np.uint8)
