import numpy as np


def perpendicular(a):
    b = np.empty_like(a)
    b[0] = a[1]
    b[1] = -a[0]
    return b


def normalize(a):
    a = np.array(a)
    return a / np.linalg.norm(a)


def create_subs_df(vu0, vu1, vu2, i_point, Vs, dist, from_i_div, to_i_div, centr_dist):
    vdiff = vu0 - vu2
    vdiffnorm = vdiff / (np.linalg.norm(vdiff))
    vnorm = perpendicular(vdiffnorm) / np.linalg.norm(vdiffnorm)
    inflex_zero = 2 + 5 * i_point

    disp_max = []
    x1 = []
    x2 = []

    lneg_1 = np.linspace(inflex_zero * -1, i_point * -1, from_i_div).tolist()
    lneg_0 = np.linspace(i_point * -1, 0, to_i_div).tolist()[1:]
    lcentr = np.linspace(0, dist, centr_dist).tolist()[1:-1]
    lpos_0 = np.linspace(dist, i_point + dist, to_i_div).tolist()[:-1]
    lpos_1 = np.linspace(i_point + dist, inflex_zero +
                         dist, from_i_div).tolist()
    l_y = lneg_1 + lneg_0 + lcentr + lpos_0 + lpos_1

    for i in l_y:
        si = -Vs / (np.sqrt(2 * np.pi) * i_point) * np.exp(-1 / 2 * i**2 / i_point**2) - \
            Vs / (np.sqrt(2 * np.pi) * i_point) * \
            np.exp(-1 / 2 * (i - dist)**2 / i_point**2)

        v_out = (vnorm * i * -1) + vu1

        disp_max.append(si)
        x1.append(v_out[0])
        x2.append(v_out[1])

    df = pd.DataFrame({'xloc': l_y, 'z': disp_max,
                       'easting_perp': x1, 'northing_perp': x2})

    return df
