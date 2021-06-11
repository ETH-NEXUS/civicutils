from civicpy import civic

# TODO: can we just load the civic database once when the package is first loaded?
success = civic.load_cache(on_stale='ignore')
if not success:
    raise ValueError("Could not load cache file from `civicpy`!")
