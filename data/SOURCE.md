# Data source

`spontaneous_swim.npz` is a conversion of `spontaneous_swim.mat` from

> J. Fernandez de Cossio Diaz and S. Cocco, *ZebrafishHMM* — Cogmaster tutorial
> on hidden Markov models applied to zebrafish trajectories.
> https://github.com/SCocco/ZebrafishHMM

523 recordings of spontaneous swimming from 39 larval zebrafish; per bout:
time, position, heading, displacement and turn angle.

The .mat file is MATLAB v7.3 (HDF5), so scipy.io cannot read it.  Conversion:

    python3 -c "
    import h5py, numpy as np
    f = h5py.File('spontaneous_swim.mat','r'); g = f['Es']
    d = {k: np.array(g[k]) for k in ('Angle','R','T','TimeBout','FishN')}
    d['x'] = np.array(g['Coordinates/x']); d['y'] = np.array(g['Coordinates/y'])
    np.savez_compressed('spontaneous_swim.npz', **d)"

Held here for reproducibility only; the repository carrying it is private, and
the public page redistributes no data.  Cite the authors above.
