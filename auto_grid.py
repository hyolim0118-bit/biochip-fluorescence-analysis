"""Geometry-first detection; intensity preprocessing is never used for measurement."""
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.optimize import least_squares


class DetectionError(ValueError):
    pass


def detect_candidates(raw):
    image = np.log1p(np.maximum(raw.astype(np.float32), 0))
    contrast = ndi.gaussian_filter(image, 1.2) - ndi.gaussian_filter(image, 7)
    noise = np.median(np.abs(contrast - np.median(contrast))) * 1.4826
    peaks = (contrast == ndi.maximum_filter(contrast, size=9)) & (contrast > max(0.025, 3*noise))
    peaks[:6] = peaks[-6:] = False
    peaks[:, :6] = peaks[:, -6:] = False
    ys, xs = np.nonzero(peaks)
    if len(xs) < 120:
        raise DetectionError('Too few spot candidates for a supported 35×35 lattice.')
    points = []
    for x, y in zip(xs, ys):
        patch = contrast[y-3:y+4, x-3:x+4]
        weights = np.maximum(patch - patch.min(), 0)
        yy, xx = np.mgrid[-3:4, -3:4]
        total = weights.sum()
        points.append([x+(weights*xx).sum()/total, y+(weights*yy).sum()/total])
    return np.asarray(points)


def estimate_basis(points):
    tree = cKDTree(points)
    _, neighbors = tree.query(points, k=min(9, len(points)))
    vectors = (points[neighbors[:, 1:]]-points[:, None]).reshape(-1, 2)
    bases = []
    for axis in (0, 1):
        other = 1-axis
        good = (vectors[:, axis] > 11) & (np.abs(vectors[:, other]) < vectors[:, axis]*0.5)
        samples = vectors[good]
        if len(samples) < 50:
            raise DetectionError('No stable horizontal/vertical lattice periodicity (rotation limit ~25 degrees).')
        # Every candidate votes equally; a saturated dust object has no extra weight.
        bins = np.floor(samples / 1.5).astype(int)
        cells, counts = np.unique(bins, axis=0, return_counts=True)
        center = (cells[np.argmax(counts)] + 0.5)*1.5
        chosen = samples[np.linalg.norm(samples-center, axis=1) < max(2.5, np.linalg.norm(center)*0.12)]
        bases.append(np.median(chosen, axis=0))
    basis = np.stack(bases, axis=1)
    if np.linalg.det(basis) <= 0 or np.linalg.cond(basis) > 2:
        raise DetectionError('Unstable lattice pitch/skew.')
    return basis


def candidate_components(points, basis):
    uv = points @ np.linalg.inv(basis).T
    tree = cKDTree(uv)
    pairs = tree.query_pairs(1.65, output_type='ndarray')
    parent = np.arange(len(points))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in pairs:
        delta = uv[b]-uv[a]
        if np.linalg.norm(delta-np.round(delta)) < 0.28:
            parent[root(a)] = root(b)
    roots = np.array([root(i) for i in range(len(points))])
    return [points[roots == label] for label, count in zip(*np.unique(roots, return_counts=True)) if count >= 120]


def design(coords, smooth=False):
    u, v = coords[:, 0]/34, coords[:, 1]/34
    terms = [np.ones(len(u)), u, v, u*v]
    if smooth:
        terms += [u*u, v*v, u*u*v, u*v*v, u**3, v**3]
    return np.stack(terms, axis=1)


def robust_fit(matrix, targets, initial=None):
    if initial is None:
        initial = np.linalg.lstsq(matrix, targets, rcond=None)[0]
    return least_squares(lambda p: (matrix @ p.reshape(initial.shape)-targets).ravel(),
                         initial.ravel(), loss='soft_l1', f_scale=1.0).x.reshape(initial.shape)


def fit_block(points, basis):
    # Refine spacing per block: a small pitch difference accumulates over 34 intervals.
    basis = estimate_basis(points)
    pitch = min(np.linalg.norm(basis, axis=0))
    uv = points @ np.linalg.inv(basis).T
    phase = np.angle(np.mean(np.exp(2j*np.pi*uv), axis=0))/(2*np.pi)
    indices = np.floor(uv-phase+0.5).astype(int)
    phase_error = np.linalg.norm(uv-phase-indices, axis=1)
    keep = phase_error < 0.4
    indices, points = indices[keep], points[keep]
    # Find the most-supported 35×35 index window, including missing interior cells.
    lower = np.percentile(indices, 2, axis=0).astype(int)
    upper = np.percentile(indices, 98, axis=0).astype(int)
    starts = []
    for axis in range(2):
        options = range(upper[axis]-34, lower[axis]+1)
        if not options:
            options = range(lower[axis]-2, lower[axis]+3)
        starts.append(list(options))
    best = None
    for x in starts[0]:
        for y in starts[1]:
            local = indices-[x, y]
            mask = ((local >= 0) & (local <= 34)).all(axis=1)
            count = len(np.unique(local[mask], axis=0))
            if best is None or count > best[0]:
                best = count, local, mask
    _, local, mask = best
    local, targets = local[mask], points[mask]
    # One geometric observation per lattice node.
    unique, first = np.unique(local, axis=0, return_index=True)
    local, targets = unique, targets[first]
    if len(local) < 120:
        raise DetectionError('Insufficient unique lattice support.')
    global_coeff = robust_fit(design(local), targets)
    row, col = np.mgrid[:35, :35]
    nodes = np.column_stack([col.ravel(), row.ravel()])
    global_grid = design(nodes) @ global_coeff
    tree = cKDTree(points)
    # Reassociate to global geometry, never relocate a node to a local brightest pixel.
    for _ in range(5):
        distance, matches = tree.query(global_grid)
        good = distance < pitch*0.23
        if good.sum() < 120:
            raise DetectionError('Global fit lost lattice support.')
        global_coeff = robust_fit(design(nodes[good]), points[matches[good]])
        global_grid = design(nodes) @ global_coeff
    distance, matches = tree.query(global_grid)
    good = distance < pitch*0.2
    if good.sum() < 120:
        raise DetectionError('Too few inliers after robust fitting.')
    coverage = good.reshape(35, 35)
    supported_rows = int(np.count_nonzero(coverage.sum(axis=1) >= 3))
    supported_cols = int(np.count_nonzero(coverage.sum(axis=0) >= 3))
    grid = global_grid.copy()
    correction = np.zeros_like(grid)
    if good.sum() >= 300 and supported_rows >= 30 and supported_cols >= 30:
        # Low-dimensional shared cubic residual field; damped, bounded and smooth.
        residual = points[matches[good]]-global_grid[good]
        coeff = robust_fit(design(nodes[good], True), residual)
        correction = design(nodes, True) @ coeff
        if np.max(np.linalg.norm(correction, axis=1)) <= pitch*0.15:
            grid += correction
        else:
            correction[:] = 0
    final_distance, _ = tree.query(grid)
    fraction = float(np.mean(final_distance < pitch*0.2))
    rms = float(np.sqrt(np.mean(final_distance[good]**2)))
    status = 'OK' if fraction >= 0.6 and supported_rows == 35 and supported_cols == 35 and rms < pitch*0.12 else 'LOW_CONFIDENCE'
    return dict(grid=grid.reshape(35, 35, 2), status=status, matched_fraction=fraction,
                rms_pixels=rms, supported_rows=supported_rows, supported_columns=supported_cols,
                pitch_x=float(np.linalg.norm(basis[:, 0])), pitch_y=float(np.linalg.norm(basis[:, 1])),
                rotation_degrees=float(np.degrees(np.arctan2(basis[1, 0], basis[0, 0]))),
                skew_degrees=float(np.degrees(np.arccos(np.clip(np.dot(basis[:, 0], basis[:, 1])/
                    np.prod(np.linalg.norm(basis, axis=0)), -1, 1)))-90),
                center=grid.mean(axis=0).tolist(),
                local_correction_max=float(np.max(np.linalg.norm(correction, axis=1))),
                node_distance=final_distance.reshape(35, 35))


def detect_grids(raw):
    points = detect_candidates(raw)
    basis = estimate_basis(points)
    blocks = []
    for component in candidate_components(points, basis):
        try:
            blocks.append(fit_block(component, basis))
        except (DetectionError, ValueError):
            continue
    if not blocks:
        raise DetectionError('No consistent 35×35 block found.')
    # Group by vertical block bands, then left to right.
    blocks.sort(key=lambda b: b['center'][1])
    bands = []
    band = min(np.linalg.norm(basis, axis=0))*17
    for block in blocks:
        if not bands or abs(block['center'][1]-bands[-1][0]['center'][1]) > band:
            bands.append([])
        bands[-1].append(block)
    blocks = [block for group in bands for block in sorted(group, key=lambda b: b['center'][0])]
    for number, block in enumerate(blocks, 1):
        block['block'] = number
    return blocks
