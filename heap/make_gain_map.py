"""make_gain_map

Functions for calculating a pixel-pixel relative gain map from pedestal variances
"""
import numpy as np
import matplotlib.pyplot as plt

def gain_from_pedvars(
        pedvar_preflip: np.ndarray,
        pedvar_postflip: np.ndarray,
        plotting: bool = False
    ):
    """
    Calculate a pixel-pixel relative gain map using pedestal variances.

    Minimizes the impact of bright stars by using two star fields (i.e. pre and postmeridian flip).
    If a pixel is bright (i.e. has an inflated pedestal variance) in one field but not the other,
    the lower of the two variances is used for that pixel instead of the average.

    Parameters:
        pedvar_preflip: pedvar map before meridian flip. numpy array with shape (32, 32) or (1024,)
        pedvar_postflip: pedvar map after meridian flip. numpy array with shape (32, 32) or (1024,)
        plotting: if True, plot the pedvar comparison, the ratio distribution and the resulting gain map

    Returns:
        relgain: relative gain map, in variance units, normalized to a mean of 1.
            numpy array with shape (32, 32)
    """

    pedvar_preflip = np.asarray(pedvar_preflip, dtype=float).reshape(32, 32)
    pedvar_postflip = np.asarray(pedvar_postflip, dtype=float).reshape(32, 32)

    # ratio of postflip to preflip pedvar, used to flag pixels with a star in one field
    ratio = np.ones_like(pedvar_preflip)
    np.divide(pedvar_postflip, pedvar_preflip, out=ratio, where=(pedvar_preflip > 0) & (pedvar_postflip > 0))

    # bright pixel thresholds: 3 sigma away from the mean of the ratio distribution
    ratio_mean = np.mean(ratio)
    ratio_std = np.std(ratio)
    brightpix_upper = ratio_mean + 3 * ratio_std
    brightpix_lower = ratio_mean - 3 * ratio_std

    # convert to variance units, so that noise is combined linearly (Var = G^2<N> + eps^2)
    var_preflip = pedvar_preflip ** 2
    var_postflip = pedvar_postflip ** 2

    # default: no star in either field - use the average of the two variances
    relgain = (var_preflip + var_postflip) / 2.0
    relgain = np.where(ratio > brightpix_upper, var_preflip, relgain)  # star in postflip
    relgain = np.where(ratio < brightpix_lower, var_postflip, relgain)  # star in preflip

    relgain = relgain / np.mean(relgain)

    if plotting:
        fig, axs = plt.subplots(2, 2, figsize=(10, 10))
        ax_scatter, ax_ratio, ax_gain, ax_map = axs.flat

        ax_scatter.scatter(pedvar_preflip, pedvar_postflip, marker='.', alpha=0.2)
        lims = (0, max(ax_scatter.get_xlim()[1], ax_scatter.get_ylim()[1]))
        x = np.array(lims)
        ax_scatter.plot(x, x * brightpix_upper, color='red')
        ax_scatter.plot(x, x * max(brightpix_lower, 0), color='red')
        ax_scatter.set_xlim(lims)
        ax_scatter.set_ylim(lims)
        ax_scatter.set_aspect('equal')
        ax_scatter.set_xlabel("Pedestal variance before meridian flip")
        ax_scatter.set_ylabel("Pedestal variance after meridian flip")
        ax_scatter.set_title("Pedestal variance comparison")

        ax_ratio.hist(ratio.ravel(), bins=40, range=(0.5, 2.0))
        ax_ratio.axvline(brightpix_upper, color='red')
        ax_ratio.axvline(brightpix_lower, color='red')
        ax_ratio.set_xlabel("Pedestal variance ratio (postflip/preflip)")
        ax_ratio.set_ylabel("Counts")
        ax_ratio.set_title("Pedestal variance ratio")

        ax_gain.hist(relgain.ravel(), bins=50, range=(0.0, 2.0))
        ax_gain.set_xlabel("Relative gain")
        ax_gain.set_ylabel("Counts")
        ax_gain.set_title("Relative gain")

        im = ax_map.imshow(relgain, vmin=0.5, vmax=1.6, origin='lower')
        ax_map.set_title("Relative gain map")
        fig.colorbar(im, ax=ax_map)

        fig.tight_layout()

        plt.show()

    return relgain
