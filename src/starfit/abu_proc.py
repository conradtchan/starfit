import isotope
import numpy as np


def ap(combine, star, db, z_max, upper_lim, z_exclude, z_lolim):
    # Mask star, remove elements Z>z_max
    mask1 = np.array([ion.Z <= z_max for ion in star.element_abundances.element])
    eval_data = star.element_abundances[mask1]

    if not upper_lim:
        mask2 = np.array([error > 0 for error in eval_data.error])
        eval_data = eval_data[mask2]

    full_abudata = np.copy(db.data.transpose())
    full_ions = np.copy(db.ions)

    if combine == 1:
        ab = []
        ab += [
            10 ** eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("C")])]
        ]  # 0
        ab += [
            10 ** eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("N")])]
        ]  # 1

        eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("C")])] = np.log10(
            ab[0] + ab[1]
        )

        er = []
        er += [eval_data.error[np.in1d(eval_data.element, [isotope.Ion("C")])]]  # 0
        er += [eval_data.error[np.in1d(eval_data.element, [isotope.Ion("N")])]]  # 1

        eval_data.error[np.in1d(eval_data.element, [isotope.Ion("C")])] = np.log10(
            (
                (ab[0] * (10 ** er[0] - 1) ** 2 + ab[1] * (10 ** er[1] - 1) ** 2)
                / (ab[0] + ab[1])
            )
            + 1
        )

        eval_data = eval_data[np.invert(np.in1d(eval_data.element, [isotope.Ion("N")]))]

        full_abudata[np.in1d(db.ions, [isotope.Ion("C")])] += full_abudata[
            np.in1d(db.ions, [isotope.Ion("N")])
        ]
        cut_abudata = full_abudata[np.invert(np.in1d(db.ions, [isotope.Ion("N")]))]
        cut_ions = full_ions[np.invert(np.in1d(db.ions, [isotope.Ion("N")]))]

    if combine == 2:
        ab = []
        ab += [
            10 ** eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("C")])]
        ]  # 0
        ab += [
            10 ** eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("N")])]
        ]  # 1
        ab += [
            10 ** eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("O")])]
        ]  # 2

        eval_data.abundance[np.in1d(eval_data.element, [isotope.Ion("C")])] = np.log10(
            ab[0] + ab[1] + ab[2]
        )

        er = []
        er += [eval_data.error[np.in1d(eval_data.element, [isotope.Ion("C")])]]  # 0
        er += [eval_data.error[np.in1d(eval_data.element, [isotope.Ion("N")])]]  # 1
        er += [eval_data.error[np.in1d(eval_data.element, [isotope.Ion("O")])]]  # 2

        eval_data.error[np.in1d(eval_data.element, [isotope.Ion("C")])] = np.log10(
            (
                (
                    ab[0] * (10 ** er[0] - 1) ** 2
                    + ab[1] * (10 ** er[1] - 1) ** 2
                    + ab[2] * (10 ** er[2] - 1) ** 2
                )
                / (ab[0] + ab[1] + ab[2])
            )
            + 1
        )

        eval_data = eval_data[
            np.invert(np.in1d(eval_data.element, [isotope.Ion("N"), isotope.Ion("O")]))
        ]

        full_abudata[np.in1d(db.ions, [isotope.Ion("C")])] += (
            full_abudata[np.in1d(db.ions, [isotope.Ion("N")])]
            + full_abudata[np.in1d(db.ions, [isotope.Ion("O")])]
        )
        cut_abudata = full_abudata[
            np.invert(np.in1d(db.ions, [isotope.Ion("N"), isotope.Ion("O")]))
        ]
        cut_ions = full_ions[
            np.invert(np.in1d(db.ions, [isotope.Ion("N"), isotope.Ion("O")]))
        ]
    else:
        cut_abudata = full_abudata
        cut_ions = full_ions

    # Mask db, remove elements that are not in the star file
    mask3 = np.in1d(cut_ions, eval_data.element, assume_unique=True)
    trimmed_db = cut_abudata[mask3]

    # z_max
    # With every element
    z_list1 = np.array([i.Z for i in db.ions])
    # Less than z_max
    z_list2 = z_list1[np.where(z_list1 <= z_max)]

    if combine == 1:
        # Without the combined elements (for example, if we are doing CN, then we get rid of N)
        z_list3 = z_list2[np.invert(np.in1d(z_list2, [isotope.Ion("N")]))]
    elif combine == 2:
        # Without the combined elements (for example, if we are doing CN, then we get rid of N)
        z_list3 = z_list2[
            np.invert(np.in1d(z_list2, [isotope.Ion("N"), isotope.Ion("O")]))
        ]
    else:
        z_list3 = z_list2
    # Excluded z, in the star's index
    z_exclude_index = np.in1d(
        [i.Z for i in eval_data.element], z_exclude, assume_unique=True
    )

    # Lower lim index for all elements
    lolim_index1 = np.in1d(z_list3, z_lolim, assume_unique=True)
    # Lower lim Z, in the star's index
    lolim_index2 = np.in1d(
        [i.Z for i in eval_data.element], z_lolim, assume_unique=True
    )
    eval_data.error[lolim_index2] = -np.abs(eval_data.error[lolim_index2])

    return (
        eval_data,
        trimmed_db,
        full_abudata,
        z_list1,
        z_list2,
        z_list3,
        z_exclude_index,
        lolim_index1,
        lolim_index2,
    )
