import math
import unittest
import ase.io
import numpy as np
import os
import tempfile

def read_allxyz(filename):
    with open(filename, "r") as f:
        file = f.read()
    xyz = file.replace(">\n", "")
    tmp_dir = tempfile.gettempdir()
    file = os.path.join(tmp_dir, "structure.xyz")
    f = open(file, "w")
    f.write(xyz)
    f.close()
    read_atoms = ase.io.read(file, index=':')
    return read_atoms


def check_dist(file, atom1, atom2):
    epsilon = 0.02
    atoms = read_allxyz(file)
    dist = math.inf
    print("Distances: ", end="")
    for frame in atoms:
        p1 = frame.positions[atom1]
        p2 = frame.positions[atom2]
        distance = np.linalg.norm(p2-p1)
        print(f"{distance} ", end="")
        if distance - dist > epsilon:
            return False
        else:
            dist = distance
    print("Checked")
    return True


class TestOrcatools(unittest.TestCase):

    def setUp(self):
        # Write test molecules
        with open("1.xyz", "w") as f:
            f.write("""14
                        Coordinates
                          C   -5.36529061559212     -0.86608651653211      0.02012559880300
                          C   -4.18873312236111     -0.21172863592802      0.01469640492368
                          O   -2.98771175733077     -0.85231795398421      0.01175393524835
                          C   -1.82886454940400     -0.01993565007077      0.00779403230037
                          C   -0.60798051030226     -0.90611930402959      0.00404676994501
                          C   0.66401787330200     -0.47436127219837     -0.00339500166889
                          H   -6.30380392490365     -0.29662450234992      0.02207186188845
                          H   -5.40862134172539     -1.96561103959184      0.02272283110365
                          H   -4.13966423956783      0.89502750626295      0.01209721206499
                          H   -1.82958515910466      0.64616306985257      0.90576602800680
                          H   -1.83540179463090      0.64580025871473     -0.89040028658903
                          H   -0.83252620785424     -1.98742712430642      0.00819536376870
                          H   0.92054039701130      0.59877377097230     -0.00783384847264
                          H   1.50556525246363     -1.18349160681131     -0.00545567132246""")
        with open("2.xyz", "w") as f:
            f.write("""14
                        Coordinates from ORCA-job linear
                          C   -3.63626535166847     -2.07423067883330     -0.77915827685653
                          C   -4.68187496059893     -0.97835746933313     -0.67262871454517
                          O   -5.58551648346604     -0.95165236328446      0.13768973670054
                          C   -1.20915474452647      0.52116794928100      0.02901963684231
                          C   -1.64302471774267     -0.75006530972356      0.09799850540314
                          C   -2.21099200393643     -1.54467154107331     -1.05069028500252
                          H   -3.95168795386129     -2.72338764579245     -1.63038545886542
                          H   -3.68894617849574     -2.69694256628652      0.13849478275028
                          H   -4.57310467256605     -0.15762315967687     -1.45211946593362
                          H   -0.80494348132481      1.03993998652113      0.91235528041416
                          H   -1.24028813722746      1.09078623887523     -0.91568405048280
                          H   -1.59544517398356     -1.27781237079394      1.06968851348365
                          H   -2.20725955611129     -0.92618785099849     -1.97541705013850
                          H   -1.55390396449078     -2.41832129888132     -1.25886261376950""")

    def tearDown(self) -> None:
        import os
        # Delete all files created during testing
        files = [
            "traj.xyz",
            "traj.allxyz",
            "traj2.xyz",
            "1.xyz",
            "2.xyz"
        ]
        for file in files:
            if os.path.isfile(file):
                os.remove(file)
            else:
                print(f"[WARN] file {file} not present after test")

    def test_dummy(self):
        pass

    def test_interpolation(self):
        from orcatools.reinterpolate import reinterpolate
        ### Interpolation of structures to allxyz
        reinterpolate(images=12, structures=["1.xyz", "2.xyz"], output="traj.allxyz", verbose=True)
        # check if the distance between C0 and C5 decreases
        assert check_dist("traj.allxyz", 0, 5)
        ### Reinterpolation of allxyz to xyz
        reinterpolate(images=16, trajectory="traj.allxyz", output="traj.xyz", verbose=True)
        # check if the distance between C0 and C5 decreases
        assert check_dist("traj.xyz", 0, 5)
        ### Interpolation of structures to xyz using distance based sampling
        reinterpolate(images=12, structures=["1.xyz", "2.xyz"], output="traj2.xyz", distance=True,
                      verbose=True)
        # check if the distance between C0 and C5 decreases
        assert check_dist("traj2.xyz", 0, 5)




