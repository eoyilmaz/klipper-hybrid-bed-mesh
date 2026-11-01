"""Tests for the hybrid eddy rapid scan + tap bed mesh."""

# Standard Library Imports
import math
import sys
import types

# Third-Party Imports
import pytest


# The real module imports Klipper's ``extras.bed_mesh`` and ``extras.probe``;
# stub them so the module under test can be imported standalone.
class BedMeshError(Exception):
    pass


class FakeZMesh:
    def __init__(self, params, name):
        self.params = dict(params)
        self.name = name
        self.matrix = None
        self.zero_ref = None

    def build_mesh(self, matrix):
        self.matrix = [list(row) for row in matrix]

    def set_zero_reference(self, x, y):
        self.zero_ref = (x, y)

    def get_mesh_params(self):
        return self.params

    def get_probed_matrix(self):
        return self.matrix


class FakeProbePointsHelper:
    instances = []

    def __init__(self, config, finalize_callback, default_points=None):
        self.printer = config.get_printer()
        self.finalize_callback = finalize_callback
        self.points = default_points
        self.gcmds = []
        self.surface = None
        FakeProbePointsHelper.instances.append(self)

    def update_probe_points(self, points, min_points):
        self.points = points

    def _move_next(self, probe_num):
        self.printer.events.append(("move", probe_num))

    def start_probe(self, gcmd):
        self.gcmds.append(gcmd)
        self.finalize_callback(
            [
                types.SimpleNamespace(bed_x=x, bed_y=y, bed_z=self.surface(x, y))
                for x, y in self.points
            ]
        )


sys.modules.setdefault(
    "extras.bed_mesh",
    types.SimpleNamespace(BedMeshError=BedMeshError, ZMesh=FakeZMesh),
)
sys.modules.setdefault(
    "extras.probe",
    types.SimpleNamespace(ProbePointsHelper=FakeProbePointsHelper),
)

# Local Imports
from extras import hybrid_bed_mesh  # noqa: E402
from extras.hybrid_bed_mesh import (  # noqa: E402
    HybridBedMesh,
    HybridPlan,
    check_algorithm,
    merge_rows,
    natural_cubic_spline,
    plan_rows,
    tap_columns,
)

# the stubs above only take effect if nothing imported the real modules first
BedMeshError = hybrid_bed_mesh.BedMeshError


# ---------------------------------------------------------------------------
# plan_rows
# ---------------------------------------------------------------------------
def test_plan_rows_returns_none_when_eddy_reaches_the_whole_mesh() -> None:
    assert plan_rows(24.3, 225, 21, 24.3) is None


def test_plan_rows_e3ng_geometry() -> None:
    plan = plan_rows(5, 225, 21, 24.3)
    assert plan.rows[:3] == pytest.approx([5, 16, 27])
    assert plan.rows[-1] == pytest.approx(225)
    assert plan.split == 2
    assert plan.step == 1
    assert plan.scan_row_count == 19
    assert plan.scan_y_dist == pytest.approx(11)


def test_plan_rows_floors_the_row_distance_like_bed_mesh() -> None:
    plan = plan_rows(5, 225, 8, 24.3)
    # 220 / 7 = 31.428... -> 31.42
    assert plan.rows[1] - plan.rows[0] == pytest.approx(31.42)


def test_plan_rows_densifies_the_scan_to_at_least_four_rows() -> None:
    # rows at 5, 18.33, 31.66, 45 -> the eddy reaches 31.66 and 45 only
    plan = plan_rows(5, 45, 4, 24.3)
    assert plan.split == 2
    assert plan.step == 3
    assert plan.scan_row_count == 4
    # the row distance is a multiple of step hundredths
    dist = round((plan.rows[1] - plan.rows[0]) * 100)
    assert dist % plan.step == 0
    assert plan.scan_y_dist * 100 == pytest.approx(round(plan.scan_y_dist * 100))


def test_plan_rows_raises_when_eddy_reaches_less_than_two_rows() -> None:
    with pytest.raises(BedMeshError, match="less than two mesh rows"):
        plan_rows(5, 30, 3, 24.3)


def test_plan_rows_raises_when_rows_are_too_close() -> None:
    with pytest.raises(BedMeshError, match="too close"):
        plan_rows(23, 26, 5, 24.3)


def test_hybrid_plan_properties() -> None:
    plan = HybridPlan(rows=[0, 10, 20, 30], split=1, step=2)
    assert plan.scan_row_count == 5
    assert plan.scan_y_dist == 5


# ---------------------------------------------------------------------------
# tap_columns
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("x_count", "tap_x_count", "expected"),
    [
        (21, 5, [0, 5, 10, 15, 20]),
        (21, 7, [0, 3, 7, 10, 13, 17, 20]),
        (5, 7, [0, 1, 2, 3, 4]),
        (21, 2, [0, 20]),
    ],
)
def test_tap_columns(x_count, tap_x_count, expected) -> None:
    assert tap_columns(x_count, tap_x_count) == expected


# ---------------------------------------------------------------------------
# natural_cubic_spline
# ---------------------------------------------------------------------------
def test_natural_cubic_spline_passes_through_the_nodes() -> None:
    xs = [0, 3, 7, 10]
    ys = [0.1, -0.2, 0.05, 0.3]
    assert natural_cubic_spline(xs, ys, xs) == pytest.approx(ys)


def test_natural_cubic_spline_reproduces_a_line() -> None:
    xs = [0, 3, 7, 10, 20]
    ys = [2 * x - 1 for x in xs]
    x_new = [i * 0.5 for i in range(41)]
    assert natural_cubic_spline(xs, ys, x_new) == pytest.approx(
        [2 * x - 1 for x in x_new]
    )


def test_natural_cubic_spline_two_points_is_linear() -> None:
    assert natural_cubic_spline([0, 10], [0, 1], [0, 5, 10]) == pytest.approx(
        [0, 0.5, 1]
    )


def test_natural_cubic_spline_follows_a_smooth_curve() -> None:
    xs = [0, 5, 10, 15, 20]
    ys = [math.sin(x / 10) for x in xs]
    x_new = list(range(21))
    result = natural_cubic_spline(xs, ys, x_new)
    # the natural end condition (no curvature) costs a little at the ends
    for x, z in zip(x_new, result):
        assert z == pytest.approx(math.sin(x / 10), abs=0.02)


# ---------------------------------------------------------------------------
# merge_rows
# ---------------------------------------------------------------------------
def test_merge_rows_shifts_the_tap_rows_onto_the_scan() -> None:
    columns = [0, 2, 4]
    scan_rows = [[0.0, 0.1, 0.2, 0.3, 0.4], [1.0] * 5]
    # taps read 0.05 higher than the scan
    tap_rows = [[0.55, 0.65, 0.75], [0.05, 0.25, 0.45]]
    matrix, corrections, deviation = merge_rows(tap_rows, columns, scan_rows)
    assert corrections == pytest.approx([0.05, 0.05, 0.05])
    assert deviation == pytest.approx(0)
    assert len(matrix) == 3
    assert matrix[0] == pytest.approx([0.5, 0.55, 0.6, 0.65, 0.7])
    assert matrix[1:] == scan_rows


def test_merge_rows_reports_the_overlap_deviation() -> None:
    columns = [0, 4]
    scan_rows = [[0.0] * 5]
    tap_rows = [[0.0, 0.0], [0.1, -0.1]]
    _, corrections, deviation = merge_rows(tap_rows, columns, scan_rows)
    assert corrections == pytest.approx([0, 0])
    assert deviation == pytest.approx(0.1)


def test_merge_rows_linear_takes_out_an_x_slope() -> None:
    columns = [0, 2, 4]
    scan_rows = [[0.0, 0.1, 0.2, 0.3, 0.4], [1.0] * 5]
    # taps read 0.05 higher than the scan plus 0.1 lower per 4 columns
    tap_rows = [[0.55, 0.6, 0.65], [0.05, 0.2, 0.35]]
    matrix, corrections, deviation = merge_rows(
        tap_rows, columns, scan_rows, linear=True
    )
    assert corrections == pytest.approx([0.05, 0.0, -0.05])
    assert deviation == pytest.approx(0)
    assert matrix[0] == pytest.approx([0.5, 0.55, 0.6, 0.65, 0.7])
    assert matrix[1:] == scan_rows


def test_merge_rows_offset_keeps_the_x_slope() -> None:
    columns = [0, 2, 4]
    scan_rows = [[0.0] * 5]
    tap_rows = [[0.0, 0.0, 0.0], [0.1, 0.0, -0.1]]
    _, corrections, deviation = merge_rows(tap_rows, columns, scan_rows)
    assert corrections == pytest.approx([0, 0, 0])
    assert deviation == pytest.approx(0.1)
    _, corrections, deviation = merge_rows(
        tap_rows, columns, scan_rows, linear=True
    )
    assert corrections == pytest.approx([0.1, 0, -0.1])
    assert deviation == pytest.approx(0)


# ---------------------------------------------------------------------------
# check_algorithm
# ---------------------------------------------------------------------------
def test_check_algorithm_keeps_a_valid_algorithm() -> None:
    assert check_algorithm("bicubic", 21, 21, 2, 2) == "bicubic"
    assert check_algorithm("lagrange", 5, 5, 2, 2) == "lagrange"


def test_check_algorithm_no_interpolation_is_direct() -> None:
    assert check_algorithm("bicubic", 21, 21, 0, 0) == "direct"


def test_check_algorithm_bicubic_with_three_points_is_lagrange() -> None:
    assert check_algorithm("bicubic", 3, 5, 2, 2) == "lagrange"


@pytest.mark.parametrize(
    ("algo", "x_count", "y_count"),
    [("lagrange", 7, 5), ("bicubic", 3, 7), ("spline", 5, 5)],
)
def test_check_algorithm_raises(algo, x_count, y_count) -> None:
    with pytest.raises(BedMeshError):
        check_algorithm(algo, x_count, y_count, 2, 2)


# ---------------------------------------------------------------------------
# HybridBedMesh
# ---------------------------------------------------------------------------
def surface(x, y):
    """A bed that curves down towards the front."""
    return 0.002 * (x - 117.5) ** 2 / 100 - 0.4 * math.exp(-y / 30)


class FakeGcmd:
    def __init__(self, params):
        self.params = params
        self.messages = []

    def error(self, msg):
        return BedMeshError(msg)

    def get(self, name, default=None):
        return self.params.get(name, default)

    def get_int(self, name, default=None, minval=None):
        value = self.params.get(name)
        return default if value is None else int(value)

    def get_float(self, name, default=None, above=None, minval=None):
        value = self.params.get(name)
        return default if value is None else float(value)

    def respond_info(self, msg):
        self.messages.append(msg)


class FakeBedMesh:
    def __init__(self):
        self.bmc = types.SimpleNamespace(
            orig_config={
                "radius": None,
                "mesh_min": (10.0, 24.3),
                "mesh_max": (220.0, 225.0),
                "x_count": 21,
                "y_count": 21,
                "algo": "bicubic",
            },
            probe_mgr=types.SimpleNamespace(
                get_zero_ref_pos=lambda: (117.5, 117.5)
            ),
        )
        self.mesh = None
        self.mesh_history = []
        self.saved = []

    def get_mesh(self):
        return self.mesh

    def set_mesh(self, mesh):
        self.mesh = mesh
        self.mesh_history.append(mesh)

    def save_profile(self, name):
        self.saved.append(name)


class FakeGcode:
    def __init__(self, bed_mesh, scan_bias):
        self.bed_mesh = bed_mesh
        self.scan_bias = scan_bias
        self.scripts = []
        self.commands = {}

    def register_command(self, name, func, desc=None):
        self.commands[name] = func

    def create_gcode_command(self, cmd, line, params):
        return FakeGcmd(params)

    def run_script_from_command(self, script):
        # rapid scan the fake surface like bed_mesh would
        self.scripts.append(script)
        params = dict(p.split("=") for p in script.split()[1:])
        min_x, min_y = (float(v) for v in params["MESH_MIN"].split(","))
        max_x, max_y = (float(v) for v in params["MESH_MAX"].split(","))
        x_cnt, y_cnt = (int(v) for v in params["PROBE_COUNT"].split(","))
        x_dist = math.floor((max_x - min_x) / (x_cnt - 1) * 100) / 100
        y_dist = math.floor((max_y - min_y) / (y_cnt - 1) * 100) / 100
        xs = [min_x + i * x_dist for i in range(x_cnt)]
        ys = [min_y + i * y_dist for i in range(y_cnt)]
        mesh = FakeZMesh(
            {
                "min_x": xs[0], "max_x": xs[-1], "min_y": ys[0],
                "max_y": ys[-1], "x_count": x_cnt, "y_count": y_cnt,
                "mesh_x_pps": 2, "mesh_y_pps": 2, "algo": "bicubic",
                "tension": 0.2,
            },
            params["PROFILE"],
        )
        mesh.build_mesh(
            [[surface(x, y) + self.scan_bias for x in xs] for y in ys]
        )
        self.bed_mesh.set_mesh(mesh)


class FakePrinter:
    def __init__(self, bed_mesh, gcode):
        self.objects = {
            "gcode": gcode,
            "bed_mesh": bed_mesh,
            "probe": types.SimpleNamespace(get_offsets=lambda: (0, 22.245, 0)),
            "toolhead": types.SimpleNamespace(
                dwell=lambda delay: self.events.append(("dwell", delay))
            ),
        }
        self.handlers = {}
        self.events = []

    def lookup_object(self, name, default=None):
        return self.objects.get(name, default)

    def load_object(self, config, name):
        return self.objects[name]

    def register_event_handler(self, event, callback):
        self.handlers[event] = callback


class FakeConfig:
    def __init__(self, printer, values):
        self.printer = printer
        self.values = values

    def get_printer(self):
        return self.printer

    def error(self, msg):
        return BedMeshError(msg)

    def _get(self, name, default):
        return self.values.get(name, default)

    def getfloatlist(self, name, default=None, count=None):
        return self._get(name, default)

    def getintlist(self, name, default=None, count=None):
        return self._get(name, default)

    def getfloat(self, name, default=None, **kwargs):
        return self._get(name, default)

    def getint(self, name, default=None, **kwargs):
        return self._get(name, default)

    def getchoice(self, name, choices, default=None):
        return choices[self._get(name, default)]

    def getsection(self, name):
        return FakeConfig(self.printer, {"position_min": 0.0})


@pytest.fixture
def hybrid_factory():
    def factory(scan_bias=0.0, **values):
        bed_mesh = FakeBedMesh()
        gcode = FakeGcode(bed_mesh, scan_bias)
        printer = FakePrinter(bed_mesh, gcode)
        hybrid = HybridBedMesh(FakeConfig(printer, values))
        hybrid.tap_helper.surface = surface
        printer.handlers["klippy:connect"]()
        return hybrid, bed_mesh, gcode

    return factory


def test_scan_min_y_defaults_to_the_eddy_reach(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory()
    assert hybrid.scan_min_y == pytest.approx(22.245)


def test_scan_min_y_from_the_config_is_kept(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory(scan_min_y=24.3)
    assert hybrid.scan_min_y == 24.3


def test_calibrate_without_front_strip_only_scans(hybrid_factory) -> None:
    hybrid, bed_mesh, gcode = hybrid_factory(scan_min_y=24.3)
    gcmd = FakeGcmd({})
    hybrid.cmd_bed_mesh_hybrid_calibrate(gcmd)
    assert len(gcode.scripts) == 1
    assert "MESH_MIN=10.000,24.300" in gcode.scripts[0]
    assert hybrid.tap_helper.gcmds == []


def test_calibrate_merges_scan_and_tap(hybrid_factory) -> None:
    hybrid, bed_mesh, gcode = hybrid_factory(
        scan_bias=0.03, scan_min_y=24.3, mesh_min=(10.0, 5.0)
    )
    gcmd = FakeGcmd({"ALGORITHM": "bicubic", "PROFILE": "default"})
    hybrid.cmd_bed_mesh_hybrid_calibrate(gcmd)

    # the scan covers Y 27 - 225 with 19 rows
    assert len(gcode.scripts) == 1
    script = gcode.scripts[0]
    assert "METHOD=rapid_scan" in script
    assert "MESH_MIN=10.000,27.000" in script
    assert "PROBE_COUNT=21,19" in script
    assert "ALGORITHM=bicubic" in script

    # 3 rows x 7 points tapped in a zig-zag, at the nozzle position
    tap_gcmd = hybrid.tap_helper.gcmds[0]
    assert tap_gcmd.params["METHOD"] == "tap"
    points = hybrid.tap_helper.points
    assert len(points) == 21
    assert [p[1] for p in points[::7]] == pytest.approx([5, 16, 27])
    assert points[0][0] < points[6][0]
    assert points[7][0] > points[13][0]

    # the partial scan mesh is cleared before tapping
    assert bed_mesh.mesh_history[1] is None

    mesh = bed_mesh.mesh
    assert mesh.params["min_y"] == pytest.approx(5)
    assert mesh.params["max_y"] == pytest.approx(225)
    assert mesh.params["y_count"] == 21
    assert mesh.params["algo"] == "bicubic"
    assert len(mesh.matrix) == 21
    assert all(len(row) == 21 for row in mesh.matrix)
    assert mesh.zero_ref == (117.5, 117.5)
    assert bed_mesh.saved == ["default"]

    # the tap rows are shifted onto the (biased) scan, so the whole mesh
    # follows the surface + scan bias, the spline is exact on tapped columns
    for i, y in enumerate([5, 16]):
        for c in (0, 3, 7, 10, 13, 17, 20):
            x = 10 + c * 10.5
            assert mesh.matrix[i][c] == pytest.approx(
                surface(x, y) + 0.03, abs=1e-6
            )
    assert any("offset: -0.0300" in m for m in gcmd.messages)


def test_calibrate_raises_on_large_overlap_offset(hybrid_factory) -> None:
    hybrid, bed_mesh, _ = hybrid_factory(
        scan_bias=0.8, scan_min_y=24.3, mesh_min=(10.0, 5.0),
        max_overlap_offset=0.5,
    )
    with pytest.raises(BedMeshError, match="max_overlap_offset"):
        hybrid.cmd_bed_mesh_hybrid_calibrate(FakeGcmd({}))


def tilted_surface(x, y):
    """The surface as the nozzle taps it, 0.2mm lower at the right edge."""
    return surface(x, y) - 0.2 * (x - 10) / 210


@pytest.mark.parametrize(
    ("params", "expected_deviation"),
    [({}, 0.0), ({"OVERLAP_ALIGNMENT": "offset"}, 0.1)],
)
def test_calibrate_overlap_alignment(
    hybrid_factory, params, expected_deviation
) -> None:
    hybrid, bed_mesh, _ = hybrid_factory(scan_min_y=24.3, mesh_min=(10.0, 5.0))
    hybrid.tap_helper.surface = tilted_surface
    gcmd = FakeGcmd(params)
    hybrid.cmd_bed_mesh_hybrid_calibrate(gcmd)
    summary = next(m for m in gcmd.messages if "X slope" in m)
    assert f"largest remaining difference: {expected_deviation:.4f}" in summary
    front = bed_mesh.mesh.matrix[0]
    if not params:
        # the slope is taken out, the front row follows the scanned surface
        for c in (0, 10, 20):
            assert front[c] == pytest.approx(surface(10 + c * 10.5, 5), abs=1e-6)
    else:
        assert front[20] - front[0] == pytest.approx(
            tilted_surface(220, 5) - tilted_surface(10, 5), abs=1e-6
        )


def test_calibrate_rejects_invalid_overlap_alignment(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory(scan_min_y=24.3, mesh_min=(10.0, 5.0))
    with pytest.raises(BedMeshError, match="OVERLAP_ALIGNMENT"):
        hybrid.cmd_bed_mesh_hybrid_calibrate(
            FakeGcmd({"OVERLAP_ALIGNMENT": "spline"})
        )


def test_calibrate_rejects_adaptive(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory()
    with pytest.raises(BedMeshError, match="ADAPTIVE"):
        hybrid.cmd_bed_mesh_hybrid_calibrate(FakeGcmd({"ADAPTIVE": "1"}))


def test_calibrate_passes_tap_parameters(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory(scan_min_y=24.3, mesh_min=(10.0, 5.0))
    hybrid.cmd_bed_mesh_hybrid_calibrate(
        FakeGcmd(
            {
                "TAP_SAMPLES": "3",
                "TAP_SAMPLE_RETRACT_DIST": "4",
                "TAP_PROBE_SPEED": "3",
                "TAP_HORIZONTAL_MOVE_Z": "6",
            }
        )
    )
    params = hybrid.tap_helper.gcmds[0].params
    assert params == {
        "METHOD": "tap",
        "SAMPLES": "3",
        "SAMPLE_RETRACT_DIST": "4.0",
        "PROBE_SPEED": "3.0",
        "HORIZONTAL_MOVE_Z": "6.0",
    }


def test_taps_wait_for_the_toolhead_to_settle(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory(scan_min_y=24.3, mesh_min=(10.0, 5.0))
    hybrid.cmd_bed_mesh_hybrid_calibrate(FakeGcmd({}))
    assert hybrid.tap_helper.settle_time == 0.5
    hybrid.tap_helper._move_next(0)
    assert hybrid.printer.events == [("move", 0), ("dwell", 0.5)]


def test_tap_settle_time_can_be_turned_off(hybrid_factory) -> None:
    hybrid, _, _ = hybrid_factory(scan_min_y=24.3, mesh_min=(10.0, 5.0))
    hybrid.cmd_bed_mesh_hybrid_calibrate(FakeGcmd({"TAP_SETTLE_TIME": "0"}))
    hybrid.tap_helper._move_next(0)
    assert hybrid.printer.events == [("move", 0)]
