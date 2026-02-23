import pytest
import numpy as np
import re
from sevnpy.io.logschema import LogSchema
from sevnpy.io.logreader import LogReader, WDLogReader


class TestLogSchema:

    def test_initialise(self):
        """
        Test if the LogSchema is initiliased

        """

        try:
            LogSchema({"name": (str, ""), "id": ("id", ""), "event": ("SN", "")})
        except Exception as err:
            pytest.fail(f"LogSchema initilisation failed with message {err}")

    def test_failed_initialisation(self):
        """
        Test if an error is raised

        """
        with pytest.raises(ValueError):
            LogSchema({"name": (str, ""), "id": ("id", ""), "event": (id, "")})

    def test_failed_insertion(self):
        """
        Test if an error is raised

        """
        ls = LogSchema({"name": (str, ""), "id": ("id", ""), "event": ("SN", "")})
        with pytest.raises(ValueError):
            ls.update({"pippo": "pluto"})


class TestLogReader:

    def test_initialisation_one(self):
        lr = LogReader(log_event="CE",
                       log_type="B",
                       body_schema=LogSchema({}))

        lr2 = LogReader(log_event="CE",
                        log_type="B",
                        body_schema={})

    def test_initialisation_two(self):
        names_list = [10320, 3245]
        ids_list = [1, 2, 3]

        lr = LogReader(log_event="CE",
                       log_type="B",
                       body_schema=LogSchema({}),
                       names=[10320, 3245],
                       IDs=[1, 2, 3])

        print(lr.summary())
        assert lr.summary()["name"]["pattern"] == "|".join(map(str, names_list))
        assert lr.summary()["ID"]["pattern"] == "|".join(map(str, ids_list))

    def test_multi_event(self):
        lr = LogReader(log_event=["CE", "BSN"],
                       log_type=["B", "S"],
                       body_schema=LogSchema({}))

    def test_multi_event_fail(self):
        with pytest.raises(TypeError):
            lr = LogReader(log_event=["CE", "BSN", 3],
                           log_type=["B", "S"],
                           body_schema=LogSchema({}))

    def test_multi_event_fail_two(self):
        with pytest.raises(ValueError):
            lr = LogReader(log_event="CE",
                           log_type=["B", "S", "D"],
                           body_schema=LogSchema({}))


class TestWDLogReader:

    def test_initilisation(self):
        from sevnpy.io.logreader import WDLogReader
        t = WDLogReader()

        assert len(t._body_schema) == 5
        assert t.log_type == "S"

    def test_one(self):
        log_entry = "S;857175750378006;0;WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3"

        tlog = WDLogReader()
        ma = re.findall(tlog.matching_pattern(capturing_names="default")[0], log_entry)
        assert ma == [('857175750378006', '0', '2.108597e+01', '2.39881', '2.39881', '1.28715', '1.28715', '3')]

    def test_reading_str(self):

        string = f"S;857175750378006;0;WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n"
        string += f"B;857175750378006;0;WDR;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n"

        df = WDLogReader().readstring(string)

        assert df["WDMass"].values[0] == pytest.approx(1.28715, 0.0001)

    def test_reading(self, tmp_path):

        fo = tmp_path / "logfile_0_tmp.dat"
        fo.write_text("S;857175750378006;0;WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n" * 10 +
                      "B;857175750378006;0;WDR;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n" * 5)

        so = WDLogReader().readfiles(str(tmp_path / "logfile_0_tmp.dat"), capturing_cols="default", nproc=1)

        df = so.data
        assert len(df) == 1
        assert np.sum(df.duplicated()) == 0
        assert df["WDMass"].values[0] == pytest.approx(1.28715, 0.0001)

        df = WDLogReader().readfiles(str(tmp_path / "logfile_0_tmp.dat"), capturing_cols=("name", "time", "Mtot_preWD"),
                                     nproc=1).data
        print(list(df.columns))
        assert list(df.columns) == ["name", "time", "Mtot_preWD", 'IDfile']

        df = WDLogReader().readfiles(str(tmp_path / "logfile_0_tmp.dat"),
                                     capturing_cols=("name", "time", 0, 1, "Mtot_preWD", 4), nproc=1).data
        assert len(df.columns) == 5  # Including the extra IDfile

    def test_reading_multi(self, tmp_path):

        nmulti = 3

        for i in range(nmulti):
            fo = tmp_path / f"logfile_{i}_tmp.dat"
            fo.write_text(f"S;857175750378006;{i};WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n" * 10)

        df = WDLogReader().readfiles(str(tmp_path / "logfile_*_tmp.dat"),
                                     capturing_cols=("name", "time", "ID", "Mtot_preWD"),
                                     nproc=1).data

        assert len(df.drop_duplicates()) == 3

    def test_reading_multi_parallel(self, tmp_path):

        nmulti = 3

        for i in range(nmulti):
            fo = tmp_path / f"logfile_{i}_tmp.dat"
            fo.write_text(f"S;857175750378006;{i};WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n" * 10)

        sevno = WDLogReader().readfiles(str(tmp_path / "logfile_*_tmp.dat"),
                                        capturing_cols=("name", "time", "ID", "Mtot_preWD"),
                                        nproc=3)

        print(sevno.data)
        print(sevno.sevn_files)
        print(sevno.data_info)
        print(sevno.columns_info)
        assert len(sevno.data.drop_duplicates()) == 3

    def test_reading_names(self, tmp_path):

        name = "857175750378006"
        ID = "0"

        fo = tmp_path / "logfile_0_tmp.dat"
        fo.write_text(f"S;{name};{ID};WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n")

        sevno = WDLogReader(names="3", IDs=ID).readfiles(str(tmp_path / "logfile_*_tmp.dat"))

        assert len(sevno.data) == 0

        sevno = WDLogReader(names=name).readfiles(str(tmp_path / "logfile_*_tmp.dat"), capturing_cols="all")

        assert len(sevno.data) == 1

        sevno = WDLogReader(names=name).readfiles(str(tmp_path / "logfile_*_tmp.dat"), capturing_cols="header")

        assert list(sevno.data.columns) == ['logtype', 'name', 'ID', 'event', 'time', 'IDfile']

        ll = WDLogReader(names="3", IDs=ID)
        sevno = ll.readfiles(str(tmp_path / "logfile_*_tmp.dat"))

        assert len(sevno.data) == 0


class TestSNLogReader:
    def test_initilisation(self):
        from sevnpy.io.logreader import SNLogReader
        t = SNLogReader()

        assert len(t._body_schema) == 11


class TestNSLogReader:
    def test_initilisation(self):
        from sevnpy.io.logreader import NSLogReader
        t = NSLogReader()

        assert len(t._body_schema) == 5
        assert t.log_type == "S"

    def test_duplicates(self):
        logfile = "S;951702780946457;20;SN;2.997599e+00;93.1927:39.1415:37.8287:35.044:6:3:206.209:0:0:0:0\n"
        logfile += "S;671196681173444;21;SN;2.902868e+00;102.005:44.7439:44.328:39.6789:6:3:460.293:0:0:0:-0\n"
        logfile += "S;871990168366889;22;SN;2.827043e+00;110.141:46.0504:45.8528:41.3679:6:3:315.714:0:0:0:0\n"
        logfile += "S;903266641389177;23;HENAKED;3.348935e+01;2.54616:2.47845:1.26626:295.107:6.98231:6:8:2.40867:1\n"
        logfile += "S;903266641389177;23;WD;3.348935e+01;2.52802:2.52802:1.26632:1.26632:3\n"
        logfile += "S;603781501433866;24;NS;2.298488e+01;5:1.19322:5.57449e+11:141.903:0.99243\n"
        logfile += "S;603781501433866;24;NS;5.598488e+01;5:1.19322:5.57449e+11:141.903:0.99243\n"

        from sevnpy.io.logreader import NSLogReader
        t = NSLogReader()

        df = t.readstring(logfile)

        assert len(df) == 1
        assert df["time"].values[0] == pytest.approx(5.598488e+01, 0.0001)

class TestQHELogReader:
    def test_initilisation(self):
        from sevnpy.io.logreader import QHELogReader
        t = QHELogReader()

        assert len(t._body_schema) == 3
        assert t.log_type == "S"

class TestHENAKEDLogReader:
    def test_initilisation(self):
        from sevnpy.io.logreader import HENAKEDLogReader
        t = HENAKEDLogReader()

        assert len(t._body_schema) == 9
        assert t.log_type == "S"


class TestCONAKEDLogReader:
    def test_initilisation(self):
        from sevnpy.io.logreader import CONAKEDLogReader
        t = CONAKEDLogReader()

        assert len(t._body_schema) == 8
        assert t.log_type == "S"


class TestBSNLogReader:

    def test_initilisation(self):
        from sevnpy.io.logreader import BSNLogReader
        t = BSNLogReader()

        assert len(t._body_schema) == 18
        assert t.log_type == "B"


class TestRLO_BEGINLogReader:

    def test_initilisation(self):
        from sevnpy.io.logreader import RLO_BEGINLogReader
        t = RLO_BEGINLogReader()

        assert len(t._body_schema) == 20
        assert t.log_type == "B"


class TestRLO_ENDLogReader:

    def test_initilisation(self):
        from sevnpy.io.logreader import RLO_ENDLogReader
        t = RLO_ENDLogReader()

        assert len(t._body_schema) == 20
        assert t.log_type == "B"


class TestCIRCLogReader:

    def test_initilisation(self):
        from sevnpy.io.logreader import CIRCLogReader
        t = CIRCLogReader()

        assert len(t._body_schema) == 4
        assert t.log_type == "B"


class TestCELogReader:

    def test_initalisation(self):
        from sevnpy.io.logreader import CELogReader
        t = CELogReader()

        assert len(t._body_schema) == 15
        assert t.log_type == "B"


class TestMERGERLogReader:

    def test_initalisation(self):
        from sevnpy.io.logreader import MERGERLogReader
        t = MERGERLogReader()

        assert len(t._body_schema) == 17
        assert t.log_type == "B"


class TestCOLLISIONLogReader:

    def test_initalisation(self):
        from sevnpy.io.logreader import COLLISIONLogReader
        t = COLLISIONLogReader()

        assert len(t._body_schema) == 12
        assert t.log_type == "B"

class TestSWALLOWEDLogReader:

    def test_initalisation(self):
        from sevnpy.io.logreader import SWALLOWEDLogReader
        t = SWALLOWEDLogReader()

        assert len(t._body_schema) == 13
        assert t.log_type == "B"


class Test_readlogfiles:

    def test_reading(self, tmp_path):
        from sevnpy.io.logreader import readlogfiles

        name = "857175750378006"
        ID = "0"

        fo = tmp_path / "logfile_0_tmp.dat"
        fo.write_text(f"S;{name};{ID};WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n")

        so = readlogfiles(str(tmp_path / "logfile_*_tmp.dat"), events="WD")

        assert len(so.data) == 1

    def test_reading_multi_one(self, tmp_path):
        from sevnpy.io.logreader import readlogfiles

        name = "857175750378006"
        ID = "0"

        fo = tmp_path / "logfile_0_tmp.dat"
        fo.write_text(f"S;{name};{ID};WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n")

        so = readlogfiles(str(tmp_path / "logfile_*_tmp.dat"), events=["WD", "CE"])

        assert len(so["CE"].data) == 0
        assert len(so["WD"].data) == 1

    def test_reading_multi_two(self, tmp_path):
        from sevnpy.io.logreader import readlogfiles

        name = "857175750378006"
        ID = "0"

        fo = tmp_path / "logfile_0_tmp.dat"
        fo.write_text(f"S;{name};{ID};WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n" +
                      f"S;{name};{ID};WD;2.108597e+01;2.39881:3.39881:1.28715:1.28715:3\n" +
                      f"S;1000;3;WD;2.108597e+01;2.39881:2.39881:1.28715:1.28715:3\n")

        so = readlogfiles(str(tmp_path / "logfile_*_tmp.dat"), events=["WD", "CE"])

        assert len(so["CE"].data) == 0
        assert len(so["WD"].data) == 2

        so = readlogfiles(str(tmp_path / "logfile_*_tmp.dat"), events=["WD", "CE"], names=name)

        assert len(so["CE"].data) == 0
        assert len(so["WD"].data) == 1


class Test_readlogstring:

    def test_reading(self, tmp_path):
        from sevnpy.io.logreader import readlogstring

        logfile = "S;951702780946457;20;SN;2.997599e+00;93.1927:39.1415:37.8287:35.044:6:3:206.209:0:0:0:0\n"
        logfile += "S;671196681173444;21;SN;2.902868e+00;102.005:44.7439:44.328:39.6789:6:3:460.293:0:0:0:-0\n"
        logfile += "S;871990168366889;22;SN;2.827043e+00;110.141:46.0504:45.8528:41.3679:6:3:315.714:0:0:0:0\n"
        logfile += "S;903266641389177;23;HENAKED;3.348935e+01;2.54616:2.47845:1.26626:295.107:6.98231:6:8:2.40867:1\n"
        logfile += "S;903266641389177;23;WD;3.348935e+01;2.52802:2.52802:1.26632:1.26632:3\n"
        logfile += "S;603781501433866;24;NS;2.298488e+01;5:1.19322:5.57449e+11:141.903:0.99243\n"
        logfile += "S;603781501433866;24;NS;5.598488e+01;5:1.19322:5.57449e+11:141.903:0.99243\n"

        dfo = readlogstring(string=logfile, events=("WD", "NS", "SN", "HENAKED", "BSN"))

        assert "WD" in dfo
        assert "NS" in dfo
        assert "SN" in dfo
        assert "HENAKED" in dfo
        assert "BSN" in dfo
        assert len(dfo["SN"]) == 3
        assert len(dfo["NS"]) == 1 # because we expected that the replicated lines are removed
        assert len(dfo["HENAKED"]) == 1
        assert len(dfo["WD"]) == 1
        assert dfo["BSN"].empty

