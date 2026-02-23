import pytest
from sevnpy.sevn import SEVNmanager, Star
import numpy as np


def check_and_initialise():
    if not SEVNmanager.check_initiliased(): SEVNmanager.init()


def compare_df_last(df1, df2, exclude=[]):
    for prop in df1.columns:
        if prop in exclude:
            continue
        elif df1[prop].values[-1] == pytest.approx(df2[prop].values[-1], abs(0.001 * df2[prop].values[-1])):
            continue
        else:
            print(f"Property {prop} is not the same in the two star (s1: {df1[prop].values[-1]},"
                  f" s2: {df2[prop].values[-1]})")
            return False

    return True


def compare_properties_last(star1, star2, exclude=[]):
    return compare_df_last(star1.getp(mode="last"), star2.getp(mode="last"), exclude=exclude)


class Test_star:

    def __del__(self):
        if SEVNmanager.check_initiliased(): SEVNmanager.close()

    def test_star_init(self):

        check_and_initialise()
        tables_info = None
        try:
            s1 = Star(Mzams=30, Z=0.001, tini="zams")
            tables_info = s1.tables_info
        except Exception as err:
            pytest.fail(f"Star initialisation failed with message {err}")

        mzams, z = self.draw_mzams_z()
        s1 = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                  star_flag="H", rseed=10, ID=5, Mass=max(mzams - 1, tables_info["min_zams"]))

        assert s1.Mzams == pytest.approx(mzams, 0.0001)
        assert s1.Zmet == pytest.approx(z, 0.0001)
        assert s1.snmodel == "rapid"
        assert s1.star_flag == "H"
        assert s1.rseed == 10
        assert s1.ID == 5
        assert s1.getp_array("Mass", mode="last")[0] == max(mzams - 1, tables_info["min_zams"])
        assert s1.used_sevnParams == SEVNmanager.get_sevnParams()
        assert s1.SEVNmanager_ID == SEVNmanager.get_ID()

    def test_star_init_HE(self):

        check_and_initialise()
        try:
            s1 = Star(Mzams=30, Z=0.001, tini="cheb", star_flag="HE")
        except Exception as err:
            pytest.fail(f"Star initialisation failed with message {err}")

        mzams, z = self.draw_mzams_z_he()
        s1 = Star(Mzams=mzams, Z=z, tini="cheb", snmodel="rapid",
                  star_flag="HE", rseed=10, ID=5)

        assert s1.Mzams == pytest.approx(mzams, 0.0001)
        assert s1.Zmet == pytest.approx(z, 0.0001)
        assert s1.snmodel == "rapid"
        assert s1.star_flag == "HE"
        assert s1.rseed == 10
        assert s1.ID == 5
        assert s1.getp_array("Mass", mode="last") == s1.getp_array("MHE", mode="last")

    def test_masses_flag_star(self):

        check_and_initialise()
        tables_info = SEVNmanager.tables_info()

        for i in range(10):
            mzams, z = self.draw_mzams_z()
            with pytest.raises(ValueError):
                s1 = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                          star_flag="HE", rseed=10, ID=5, MHE=0.9 * mzams)

            mzams, z = self.draw_mzams_z_he()
            s1 = Star(Mzams=mzams, Z=z, tini="cheb", snmodel="rapid",
                      star_flag="HE", rseed=10, ID=5, Mass=0.9 * mzams)
            assert s1.getp_array("Mass", mode="last") == pytest.approx(0.9 * mzams, 0.001)

            # Check automatically correct star flag
            mzams, z = self.draw_mzams_z_he()
            mzams = min(mzams,tables_info["max_zams"])
            mzams = max(mzams,tables_info["min_zams"])
            z = min(z,tables_info["max_z"])
            z = max(z,tables_info["min_z"])
            s1 = Star(Mzams=mzams, Z=z, tini="cheb", snmodel="rapid",
                      star_flag="H", rseed=10, ID=5, Mass=0.9 * mzams, MHE=0.9 * mzams)
            assert s1.star_flag == "HE"

            s1 = Star(Mzams=mzams, Z=z, tini="cheb", snmodel="rapid",
                      star_flag="HE", rseed=10, ID=5, Mass=1 * mzams, MHE=0.9 * mzams)

            assert s1.star_flag == "H"

    def test_star_basic_evolve(self):
        check_and_initialise()
        s1 = Star(Mzams=50, Z=0.0001, tini="zams", snmodel="rapid",
                  star_flag="H", rseed=10, ID=5)

        ret, extf = s1._evolve_basic(tstart="sheb", just_init=True,
                                     Mass=20, MHE=10, MCO=2)

        assert ret["Mass"] == pytest.approx(20, 0.001)
        assert ret["MHE"] == pytest.approx(10, 0.001)
        assert ret["MCO"] == pytest.approx(2, 0.001)

    def test_evolve(self):
        check_and_initialise()
        s1 = Star(Mzams=50, Z=0.0001, tini="zams", snmodel="rapid",
                  star_flag="H", rseed=10, ID=5)
        s1.evolve(tend="end")
        assert s1.getp_array("Phase", mode="last") == 7

        s1 = Star(Mzams=50, Z=0.0001, tini="zams", snmodel="rapid",
                  star_flag="H", rseed=10, ID=5)
        s1.evolve(tend=10)
        assert s1.getp_array("Worldtime", mode="last") == 10

    def test_evolve_for(self):

        check_and_initialise()
        for i in range(10):
            mzams, z = self.draw_mzams_z()
            s1 = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                      star_flag="H", rseed=10, ID=5)
            Phasec=100
            tevolve = 0.9 * s1.tlife
            i=0
            # Do not consider remnant
            while Phasec>=7:
                tevolve = (0.9-0.05*i) * s1.tlife
                s1.evolve(tend=tevolve)
                Phasec=s1.getp("PhaseBSE",mode="last").values
                i+=1

            s2 = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                      star_flag="H", rseed=10, ID=5)
            s2.evolve_for(0.2 * tevolve)
            s2.evolve_for(0.5 * tevolve)
            s2.evolve_for(0.3 * tevolve)

            print(mzams,z,tevolve)
            assert compare_properties_last(s1, s2)

    def test_evolve_for_he(self):

        check_and_initialise()
        mzams, z = self.draw_mzams_z_he()
        s1 = Star(Mzams=mzams, Z=z, tini="cheb", snmodel="rapid",
                  star_flag="HE", rseed=10, ID=5)
        s1.evolve(tend=0.1)

        s2 = Star(Mzams=mzams, Z=z, tini="cheb", snmodel="rapid",
                  star_flag="HE", rseed=10, ID=5)
        s2.evolve_for(0.05)
        s2.evolve_for(0.05)

        assert compare_properties_last(s1, s2)

    def test_get_remnant(self):
        check_and_initialise()
        s1 = Star(Mzams=50, Z=0.0001, tini="zams", snmodel="rapid",
                  star_flag="H", rseed=10, ID=5)
        Phase = s1.get_remnant()["Phase"].values
        assert len(Phase) == 1
        assert Phase[0] == 7

    def test_look_at_track(self):
        check_and_initialise()

        for i in range(10):
            mzams, z = self.draw_mzams_z()
            s1 = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                      star_flag="H", rseed=10, ID=5)
            Tevolve = 0.5 * s1.tlife
            Tevolves = Tevolve - s1.getp_array("Localtime",mode="last")[0]
            s1.evolve(Tevolves)

            s2 = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                      star_flag="H", rseed=10, ID=5, Mass=0.5 * mzams)
            s2.evolve(Tevolves)

            st = Star(Mzams=mzams, Z=z, tini="zams", snmodel="rapid",
                      star_flag="H", rseed=10, ID=5)
            dl = st.look_at_track(Tevolve)

            if s1.getp_array("Phase",mode="last") < 4:
                print(s1.getp(mode="last"))
                print(s2.getp(mode="last"))
                print(dl)
                assert compare_df_last(s1.getp(mode="last"), dl, exclude=["Worldtime", ])
                assert compare_df_last(s2.getp(mode="last"), dl, exclude=["Worldtime", ]) == False
    def draw_mzams_z(self):
        check_and_initialise()
        tables_info = SEVNmanager.tables_info()
        mzams = np.random.uniform(tables_info["min_zams"], tables_info["max_zams"])
        z = np.random.uniform(tables_info["min_z"], tables_info["max_z"])


        return mzams, z

    def draw_mzams_z_he(self):
        check_and_initialise()
        tables_info = SEVNmanager.tables_info()
        mzams = np.random.uniform(tables_info["min_zams_he"], tables_info["max_zams_he"])
        z = np.random.uniform(tables_info["min_z_he"], tables_info["max_z_he"])

        return mzams, z
