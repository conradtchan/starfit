import starfit

# starfit.Single(
#         "HE1327-2326.dat", "znuc.S4.star.el.y.stardb.gz", silent=True, z_max=30)

if __name__ == "__main__":
    starfit.Multi(
        "HE1327-2326.dat",
        "znuc.S4.star.el.y.stardb.gz",
        z_max=30,
        sol_size=2,
        threads=8,
    )
